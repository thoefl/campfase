! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_polycrystal
  use mod_base_functions
  use mod_regular_grid
  use mod_particle
  use mod_normal_dist
  use mod_crystal_indices
  use mod_globals, only: verbosity_level
  use mod_user_input
  use mod_random
  implicit none
  private

  type voronoi_t
    integer :: attempts
    real(wp) :: free_fill_ratio
  end type voronoi_t

  type distribution_t
    real(wp) :: median, sigma
  end type distribution_t

  type, public :: polycrystal_t
    type(voronoi_t) :: voronoi
    type(distribution_t) :: radius_distribution
    real(wp) :: rel_density
  contains
    procedure, public :: construct
  end type polycrystal_t

  real(wp), dimension(:), allocatable :: dist_map

contains

  subroutine construct(this, grid, idx, eta_labels)
  class(polycrystal_t), intent(in) :: this
    type(regular_grid), intent(in) :: grid
    type(crystal_indices), intent(in) :: idx
    integer, dimension(:), intent(out) :: eta_labels
    type(spheric_particle), dimension(:), allocatable :: seeds, best_seeds
    integer, dimension(:), allocatable :: grains, best_grains, intersection_map, seeds_id, seed_eta_idx, eta_count, &
      best_seed_eta_idx, seed_n_neighbours, randseed_backup
    real(wp), dimension(:), allocatable :: best_dist_map, seed_min_dist, eta_min_dist, &
      best_eta_min_dist
    logical, dimension(:), allocatable :: fill_mask, eta_filled, possible_coordinates
    real(wp) :: pairs_per_img, get_radius(1), radius, seed_area, total_area, max_total_area, &
      target_area, sigma_area, mu_area, sigma_n_seeds, mu_n_seeds, dist, grow_radius, particle_relative_density, &
      global_eta_min_dist, best_global_eta_min_dist, neighbour_node_max_dist
    integer :: i, j, pair_lo, pair_hi, rand, node, seed_idx, eta_idx, best_img,&
      n_seeds, dist_try, n_seeds_max, max_contacts, n_seed_pairs, &
      pair_idx, seed_1, seed_2, n_eta, min_n_eta
    type(int_ragged_matrix), allocatable :: surf_p(:)
    type(real_vector), allocatable :: seed_dist(:)
    real(wp), parameter :: n_seeds_max_quantile = 0.99_wp
    character(*), parameter :: proc_name = 'construct'

    allocate(grains(grid%n_nodes_total))
    allocate(intersection_map(grid%n_nodes_total))
    allocate(dist_map(grid%n_nodes_total))
    allocate(possible_coordinates(grid%n_nodes_total))
    allocate(eta_filled(grid%n_nodes_total))
    allocate(fill_mask(grid%n_nodes_total))

    neighbour_node_max_dist = norm2(grid%axes%spacing)
    n_eta = idx%eta%last - idx%eta%first + 1
    eta_filled = .false.
    eta_labels = 0
    grains = 0
    n_seeds = 0
    target_area = this%rel_density*grid%n_nodes_total
    randseed_backup = get_random_seed()

    ! Calculate the maximum amount of required particles as the 99th percentile of the total area distribution
    ! This is assuming that all particles are of equal size as defined by the given percentile (and thus inaccurate)
    sigma_area = sqrt(2.0_wp)*this%radius_distribution%sigma
    mu_area = 2*log(this%radius_distribution%median) + log(pi)
    sigma_n_seeds = sigma_area
    mu_n_seeds = log(target_area) - mu_area
    n_seeds_max = nint(exp(mu_n_seeds + sqrt(2.0_wp)*sigma_n_seeds*erfinv(2*n_seeds_max_quantile - 1)))
    allocate(seeds(n_seeds_max))

    ! Generate unique random number streams for each image
    call init_random_seed(image_distinct = .true., repeatable = .true.)

    max_total_area = 0
    do dist_try = 1, this%voronoi%attempts
      if(verbosity_level >= 2) call print_message('voronoi: Attempt ' // convert_to_char(dist_try) // '/' //&
        convert_to_char(this%voronoi%attempts))
      seed_idx = 0
      total_area = 0
      dist_map = huge(1.0_wp)
      grains = 0
      intersection_map = 0
      do
        if(seed_idx >= n_seeds_max) then
          call print_message('WARNING: Number of maximum seeds exceeded during voronoi construction!')
          exit
        end if
        get_radius =&
          exp(rnorm_box_muller_vec(1, log(this%radius_distribution%median), this%radius_distribution%sigma))
        radius = get_radius(1)
        seed_area = radius**2*pi
        if(total_area + seed_area > target_area) exit

        ! Find possible nodes for the next seed, given its radius
        ! Rules:
        !  - Particles must not overlap
        !  - Particles up to a certain fraction of the total area can be placed freely
        !  - Other particles are randomly placed for maximum amount of contacts with other particles

        if(total_area > this%voronoi%free_fill_ratio*grid%n_nodes_total) then
          possible_coordinates = dist_map > radius .and. dist_map <= (radius + neighbour_node_max_dist)
          max_contacts = maxval(intersection_map, possible_coordinates)
          where(intersection_map < max_contacts)
            possible_coordinates = .false.
          end where
        else
          possible_coordinates = dist_map > radius
        end if

        ! If there are no possible locations, exit
        ! Else, randomly select one of the possible locations
        if(.not. any(possible_coordinates)) then
          exit
        else
          seed_idx = seed_idx + 1
          if(this_image() == 1 .and. verbosity_level >= 4) then
            call print_message('voronoi: Placing particle ' // convert_to_char(seed_idx))
          end if
          total_area = total_area + seed_area
          rand = random_int(1, count(possible_coordinates))
          node = nth_true(possible_coordinates, rand)
          if(node == 0) error stop 'ERROR: voronoi: Could not find nth true node!'
          call seeds(seed_idx)%init(grid, seed_idx, radius, grid%linear_to_dim_idx(node))
          do concurrent(node = 1:grid%n_nodes_total)
            dist = grid%point_distance(node, seeds(seed_idx)%coordinates) - seeds(seed_idx)%radius
            if(dist < dist_map(node)) then
              if(dist + neighbour_node_max_dist > dist_map(node)) then
                intersection_map(node) = intersection_map(node) + 1
              else
                intersection_map(node) = 0
              end if
              grains(node) = seed_idx
              dist_map(node) = dist
            else if(dist <= dist_map(node) + neighbour_node_max_dist) then
              intersection_map(node) = intersection_map(node) + 1
            end if
          end do
        end if
      end do
      if(total_area > max_total_area) then
        best_seeds = seeds
        best_grains = grains
        best_dist_map = dist_map
        n_seeds = seed_idx
        max_total_area = total_area
      end if
    end do
    call move_alloc(best_grains, grains)
    call move_alloc(best_dist_map, dist_map)
    best_img = co_maxloc(max_total_area)
    if(best_img <= 0) error stop 'ERROR: voronoi: Could not determine densest particle distribution across images!'
    if(verbosity_level >= 5) call print_message('voronoi: Best image.', co_img_opt = best_img)
    call co_broadcast(n_seeds, best_img)
    deallocate(seeds)
    allocate(seeds_id(n_seeds))
    if(this_image() == best_img) then
      seeds = best_seeds(1:n_seeds)
      seeds_id = seeds%id
    else
      allocate(seeds(n_seeds))
    end if
    deallocate(best_seeds)
    call co_broadcast(grains, best_img)
    call co_broadcast(dist_map, best_img)
    call co_broadcast(seeds_id, best_img)

    particle_relative_density = relative_density_at_radius(0.0_wp)
    if(verbosity_level >= 2) then
      call print_message('Generated particle density: ' // convert_to_char(particle_relative_density), co_img_opt = 1)
    end if

    if(this%rel_density < 1.0_wp) then
      grow_radius = findroot_bisect(relative_density_at_radius, 0.0_wp, maxval(dist_map),&
        spacing(maxval(dist_map)), y_target_opt = this%rel_density, y_tol_opt = 2.0E-5_wp)
    else
      grow_radius = huge(1.0_wp)
    end if

    if(this_image() == 1) call output_image(merge(grains, 0, dist_map <= 0.0_wp), 'grains_particle', grid%axes(1)%n_nodes)
    grains = merge(grains, 0, dist_map <= grow_radius)
    if(this_image() == 1) call output_image(grains, 'grains', grid%axes(1)%n_nodes)
    if(this_image() == 1) call output_image(grains > 0, 'grains_fill', grid%axes(1)%n_nodes)

    if(verbosity_level >= 2) call print_message('voronoi: Calculating grain pair distances...', co_img_opt = 1)
    if(n_seeds > n_eta) then
      allocate(surf_p(n_seeds))
      n_seed_pairs = n_seeds*(n_seeds - 1)//2
      allocate(seed_dist(n_seeds - 1))
      do seed_idx = 1, n_seeds
        fill_mask = grains == seed_idx
        surf_p(seed_idx)%row = grid%get_surface_points(fill_mask)
        if(seed_idx < n_seeds) then
          allocate(seed_dist(seed_idx)%vals(seed_idx+1:n_seeds))
          seed_dist(seed_idx)%vals = huge(0.0_wp)
        end if
      end do

      call img_partition(n_seed_pairs, pairs_per_img, pair_lo, pair_hi)
      if(verbosity_level >= 5) then
        call print_message('Handling seed pairs ' // convert_to_char(pair_lo) // ' to ' // convert_to_char(pair_hi) &
          // ' out of ' // convert_to_char(n_seed_pairs))
      end if

      ! Calculate the minimum distance between all grain pairs, distributed among the images
      seed_1 = 1
      seed_2 = 2
      do pair_idx = 1, n_seed_pairs
        if(pair_idx >= pair_lo .and. pair_idx <= pair_hi) then
          do concurrent(i = 1:size(surf_p(seed_1)%row), j = 1:size(surf_p(seed_2)%row))
            seed_dist(seed_1)%vals(seed_2) = min(seed_dist(seed_1)%vals(seed_2), &
              grid%point_distance(surf_p(seed_1)%row(i)%v, surf_p(seed_2)%row(j)%v))
          end do
        end if
        if(seed_1 == n_seeds - 1 .and. seed_2 == n_seeds) exit
        seed_2 = seed_2 + 1
        if(seed_2 > n_seeds) then
          seed_1 = seed_1 + 1
          seed_2 = seed_1 + 1
        end if
      end do

      ! Broadcast the locally calculated distances to all images, so that all images have the information for all grains
      do seed_idx = 1, n_seeds - 1
        call co_min(seed_dist(seed_idx)%vals)
        if(any(seed_dist(seed_idx)%vals > norm2(grid%axes%length))) then
          call error_msg(proc_name, 'Seed distance at index ' // convert_to_char(seed_idx) &
            // ' is larger than physically possible in the given grid (expected maximum: ' &
            // convert_to_char(norm2(grid%axes%length)) // ' < ' // convert_to_char(seed_dist(seed_idx)%vals) // ')')
        else if(any(seed_dist(seed_idx)%vals < 0)) then
          call error_msg(proc_name, 'Seed distance at index ' // convert_to_char(seed_idx) // ' is negative: ' &
            // convert_to_char(seed_dist(seed_idx)%vals))
        end if
      end do


      ! Determine the number of direct neighbours for all grains
      allocate(seed_n_neighbours(n_seeds))
      seed_n_neighbours = 0
      do concurrent(seed_1 = 1:n_seeds - 1)
        do concurrent(seed_2 = seed_1 + 1:n_seeds)
          if(seed_dist(seed_1)%vals(seed_2) <= neighbour_node_max_dist) then
            seed_n_neighbours(seed_1) = seed_n_neighbours(seed_1) + 1
            seed_n_neighbours(seed_2) = seed_n_neighbours(seed_2) + 1
          end if
        end do
      end do
      min_n_eta = maxval(seed_n_neighbours)
      if(n_eta < min_n_eta) then
        if(this_image() == 1) then
          call warning_msg(proc_name, 'Number of order parameter fields (' // convert_to_char(n_eta) &
            // ') is lower than the minimum amount required to prevent merging of grains (' &
            // convert_to_char(min_n_eta) // ')')
          ! This doesn't work for some reason, maybe because inside team?
        end if
        call millisleep(5000)
      end if

      allocate(eta_min_dist(n_eta))
      allocate(seed_min_dist(n_eta))
      allocate(seed_eta_idx(n_seeds))

      if(verbosity_level >= 2) then
        call print_message('voronoi: Distributing grains among order parameter fields...', co_img_opt = 1)
      end if
      best_global_eta_min_dist = 0.0_wp
      do dist_try = 1, 1000
        eta_min_dist = huge(0.0_wp)
        seed_eta_idx = 0
        eta_idx = 1
        do
          seed_idx = random_int(1, n_seeds)
          if(seed_eta_idx(seed_idx) > 0) cycle
          seed_eta_idx(seed_idx) = eta_idx
          if(eta_idx == n_eta) exit
          eta_idx = eta_idx + 1
        end do

        do seed_idx = 1, n_seeds
          if(seed_eta_idx(seed_idx) > 0) cycle
          seed_min_dist = maxval(grid%axes%length)
          do concurrent(seed_1 = 1:n_seeds, seed_1 /= seed_idx .and. seed_eta_idx(seed_1) > 0)
            eta_idx = seed_eta_idx(seed_1)
            seed_min_dist(eta_idx) = min(seed_min_dist(eta_idx), &
              seed_dist(min(seed_1,seed_idx))%vals(max(seed_1,seed_idx)))
          end do
          if(any(seed_min_dist >= eta_min_dist .and. eta_min_dist > neighbour_node_max_dist)) then
            do i = 1, n_eta
              eta_idx = maxloc(seed_min_dist, 1)
              if(seed_min_dist(eta_idx) >= eta_min_dist(eta_idx)) then
                exit
              else
                seed_min_dist(eta_idx) = -seed_min_dist(eta_idx)
              end if
              if(i == n_eta) call error_msg(proc_name)
            end do
          else
            eta_idx = maxloc(seed_min_dist, 1)
          end if
          seed_eta_idx(seed_idx) = eta_idx
          eta_min_dist(eta_idx) = min(eta_min_dist(eta_idx), seed_min_dist(eta_idx))
        end do

        global_eta_min_dist = minval(eta_min_dist, eta_min_dist > 0)
        if(global_eta_min_dist > best_global_eta_min_dist) then
          best_eta_min_dist = eta_min_dist
          best_global_eta_min_dist = global_eta_min_dist
          best_seed_eta_idx = seed_eta_idx
        end if
      end do

      best_img = co_maxloc(best_global_eta_min_dist)
      call co_broadcast(best_global_eta_min_dist, best_img)
      call co_broadcast(best_seed_eta_idx, best_img)
      call co_broadcast(best_eta_min_dist, best_img)
      deallocate(seed_eta_idx)
      deallocate(eta_min_dist)
      call move_alloc(best_seed_eta_idx, seed_eta_idx)
      call move_alloc(best_eta_min_dist, eta_min_dist)

      allocate(eta_count(n_eta))
      do eta_idx = 1, n_eta
        eta_count(eta_idx) = count(seed_eta_idx == eta_idx)
      end do

      do i = 1, n_eta
        if(count(eta_count > 0) == 0) exit
        eta_idx = maxloc(eta_count, 1)
        eta_count(eta_idx) = -1
        where(seed_eta_idx == eta_idx)
          seed_eta_idx = -i
        end where
      end do
      seed_eta_idx = abs(seed_eta_idx)

      do seed_idx = 1, n_seeds
        where(grains == seed_idx)
          eta_labels = seed_eta_idx(seed_idx)
        end where
      end do

      if(verbosity_level >= 3) then
        call print_message('Distributed ' // convert_to_char(n_seeds) // ' grains among ' &
          // convert_to_char(n_eta) // ' order parameter fields with minimum distance = ' &
          // convert_to_char(best_global_eta_min_dist), co_img_opt = 1)
      end if
      if(verbosity_level >= 4 .and. this_image() == 1) then
        call print_message('Data by order parameter index: ')
        write(output_unit, '(a)', advance = 'yes') 'eta_idx  n_grains  min_dist  fill_ratio'
        do eta_idx = 1, n_eta
          write(output_unit, '(i7)', advance = 'no') eta_idx
          write(output_unit, '(a)', advance = 'no') '  '
          write(output_unit, '(i8)', advance = 'no') count(seed_eta_idx == eta_idx)
          write(output_unit, '(a)', advance = 'no') '  '
          write(output_unit, '(f8.4)', advance = 'no') eta_min_dist(eta_idx)
          write(output_unit, '(a)', advance = 'no') '  '
          write(output_unit, '(f8.4)', advance = 'no') 100*count(eta_labels == eta_idx)/real(grid%n_nodes_total, wp)
          write(output_unit, '(a)', advance = 'yes') ' %'
        end do
      end if
    else
      eta_labels = grains
    end if

    where(eta_labels > 0)
      eta_labels = eta_labels - 1 + idx%eta%first
    end where
    if(this_image() == 1) call output_image(eta_labels, 'grains_eta_idx', grid%axes(1)%n_nodes)
    call restore_random_seed(randseed_backup)
    deallocate(dist_map)
  end subroutine construct

  pure real(wp) function relative_density_at_radius(r)
    real(wp), intent(in) :: r

    relative_density_at_radius = true_ratio(dist_map <= r)
  end function relative_density_at_radius
end module mod_polycrystal
