! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

!TODO: Cases where spacing /= 1 need checking

module mod_regular_grid
  use mod_base_functions
  use mod_boundary_condition
  use mod_axis
  implicit none

  private
  character(*), parameter :: module_name = 'mod_regular_grid'

  type, public :: regular_grid
    integer :: n_dims, n_nodes_total
    type(axis), allocatable, dimension(:) :: axes
  contains
    procedure, public :: init
    procedure, public :: is_equispaced
    procedure, public :: linear_to_dim_idx
    procedure, public :: dim_to_linear_idx
    procedure, private :: linear_idx_point_distance
    procedure, private :: linear_dim_idx_point_distance
    procedure, private :: dim_idx_point_distance
    procedure, private :: real_point_distance
    generic, public :: point_distance => linear_idx_point_distance, dim_idx_point_distance,&
      linear_dim_idx_point_distance, real_point_distance
    procedure, public :: get_circle_mask
    procedure, public :: get_circle_nodes
    procedure, public :: get_circle_sections
    procedure, public :: get_surface_points
    procedure, public :: offset_dim_idx
    procedure, private :: increment_idx_full_range
    procedure, private :: increment_idx_limits
    generic, public :: increment_idx => increment_idx_full_range, increment_idx_limits
  end type regular_grid

contains

  pure subroutine init(this, axes)
  class(regular_grid), intent(inout) :: this
    type(axis), intent(in) :: axes(:)
    character(*), parameter :: proc_name = 'init'

    this%n_dims = size(axes)
    this%axes = axes
    this%n_nodes_total = product(axes%n_nodes)
  end subroutine init

  pure logical function is_equispaced(this)
  class(regular_grid), intent(in) :: this

    is_equispaced = all_equal(this%axes%spacing)
  end function is_equispaced

  pure function linear_to_dim_idx(this, lin_idx) result(dim_idx)
  class(regular_grid), intent(in) :: this
    integer, intent(in) :: lin_idx
    integer, dimension(this%n_dims) :: dim_idx
    integer :: dim, remainder, quotient, divisor

    remainder = lin_idx - 1
    do dim = this%n_dims, 2, -1
      divisor = product(this%axes(1:(dim - 1))%n_nodes)
      quotient = remainder/divisor
      dim_idx(dim) = quotient + 1
      remainder = remainder - quotient*divisor
    end do
    dim_idx(1) = remainder + 1
  end function linear_to_dim_idx

  pure integer function dim_to_linear_idx(this, dim_idx) result(lin_idx)
  class(regular_grid), intent(in) :: this
    integer, intent(in), dimension(:) :: dim_idx
    integer :: dim

    lin_idx = dim_idx(1)
    do dim = 2, this%n_dims
      lin_idx = lin_idx + (dim_idx(dim) - 1)*product(this%axes(:dim - 1)%n_nodes)
    end do
  end function dim_to_linear_idx

  pure real(wp) function linear_idx_point_distance(this, lin_idx_1, lin_idx_2) result(dist)
  class(regular_grid), intent(in) :: this
    integer, intent(in) :: lin_idx_1, lin_idx_2
    integer, dimension(this%n_dims) :: p_1, p_2

    p_1 = this%linear_to_dim_idx(lin_idx_1)
    p_2 = this%linear_to_dim_idx(lin_idx_2)
    dist = this%dim_idx_point_distance(p_1, p_2)
  end function linear_idx_point_distance

  pure real(wp) function dim_idx_point_distance(this, p_1, p_2) result(dist)
  class(regular_grid), intent(in) :: this
    integer, intent(in), dimension(:) :: p_1, p_2
    integer, dimension(this%n_dims) :: dist_nodes_dim

    dist_nodes_dim = abs(p_2 - p_1)
    where(this%axes%bound_cond%is_periodic)
      where(dist_nodes_dim > this%axes%n_nodes/2)
        dist_nodes_dim = this%axes%n_nodes - dist_nodes_dim
      end where
    end where
    if(this%is_equispaced()) then
      dist = sqrt(real(sum(dist_nodes_dim**2), wp))*this%axes(1)%spacing
    else
      dist = norm2(dist_nodes_dim*this%axes%spacing)
    end if
  end function dim_idx_point_distance

  pure real(wp) function linear_dim_idx_point_distance(this, lin_idx_1, p_2) result(dist)
  class(regular_grid), intent(in) :: this
    integer, intent(in) :: lin_idx_1
    integer, intent(in), dimension(:) :: p_2
    integer, dimension(this%n_dims) :: p_1

    p_1 = this%linear_to_dim_idx(lin_idx_1)
    dist = this%dim_idx_point_distance(p_1, p_2)
  end function linear_dim_idx_point_distance

  pure function get_circle_sections(this, p_center, radius) result(sections)
  class(regular_grid), intent(in) :: this
    integer, intent(in), dimension(:) :: p_center
    real(wp), intent(in) :: radius
    real(wp) :: radius_nodes
    type(array_index), dimension(:), allocatable :: sections, sections_temp
    integer, dimension(this%n_dims) :: p_low, p_high, p
    integer :: n_sections_upper_lim, n_sections, first, last
    logical :: inside
    character(*), parameter :: proc_name = 'get_circle_sections'

    if(.not. this%is_equispaced()) then
      call error_msg(proc_name, 'Not defined for non-equispaced grids (spacings: ' &
        // convert_to_char(this%axes%spacing), module_name_opt = module_name)
    end if
    radius_nodes = radius/this%axes(1)%spacing
    if(radius_nodes < 0.0_wp) then
      allocate(sections(0))
      return
    else if(radius_nodes < 1.0_wp) then
      allocate(sections(1))
      first = this%dim_to_linear_idx(p_center)
      last = first
      sections(1) = array_index(first, last)
      return
    else if(all(radius_nodes > this%axes%n_nodes &
      .or. (this%axes%bound_cond%is_periodic .and. radius_nodes >= 0.5_wp*this%axes%n_nodes))) then
      allocate(sections(1))
      first = 1
      last = this%n_nodes_total
      sections(1) = array_index(first, last)
      return
    end if

    ! Factor of 2 comes from possible wrapping over first dimension
    select case(this%n_dims)
    case(1)
      n_sections_upper_lim = 2
    case(2)
      n_sections_upper_lim = 2*ceiling(2*radius_nodes)
    case(3)
      n_sections_upper_lim = 2*ceiling(radius_nodes**2*pi)
    case default
      n_sections_upper_lim = 0
      call error_msg(proc_name, 'Not defined for the given amount of dimensions (' // convert_to_char(this%n_dims) // ')')
    end select
    allocate(sections(n_sections_upper_lim))

    !FIXME: this implicitly assumes periodic boundary conditions
    where(radius_nodes < 0.5_wp*this%axes%n_nodes)
      p_low = mod_idx(p_center - ceiling(radius_nodes), this%axes%n_nodes)
      p_high = mod_idx(p_center + ceiling(radius_nodes), this%axes%n_nodes)
    elsewhere
      p_low = 1
      p_high = this%axes%n_nodes
    end where
    p = p_low
    n_sections = 0
    first = 0
    last = 0
    inside = .false.
    do
      if(.not. inside .and. this%point_distance(p, p_center) <= radius_nodes) then
        inside = .true.
        n_sections = n_sections + 1
        if(n_sections > n_sections_upper_lim) then
          call error_msg(proc_name,&
            'Number of sections (n_sections) is greater than the assumed upper limit (n_sections_upper_lim)')
        end if
        first = this%dim_to_linear_idx(p)
      else if(inside .and. this%point_distance(p, p_center) > radius_nodes) then
        inside = .false.
        last = this%dim_to_linear_idx(p) - 1
        sections(n_sections) = array_index(first, last)
      end if

      if(inside .and. (p(1) == this%axes(1)%n_nodes .or. p(1) == p_high(1))) then
        inside = .false.
        last = this%dim_to_linear_idx(p)
        sections(n_sections) = array_index(first, last)
      end if

      if(all(p == p_high)) exit
      call this%increment_idx(p_low, p_high, p)
    end do
    if(n_sections /= n_sections_upper_lim) then
      ! Resize to actual dimensions
      allocate(sections_temp(n_sections))
      sections_temp = sections(1:n_sections)
      deallocate(sections)
      call move_alloc(sections_temp, sections)
    end if
  end function get_circle_sections

  pure function get_circle_nodes(this, p_center, radius) result(nodes)
  class(regular_grid), intent(in) :: this
    integer, intent(in), dimension(:) :: p_center
    real(wp), intent(in) :: radius
    real(wp) :: radius_nodes
    integer, dimension(:), allocatable :: nodes, nodes_temp
    integer, dimension(this%n_dims) :: p_low, p_high, p
    integer :: dim, n_nodes_upper_lim, n_nodes, i
    character(*), parameter :: proc_name = 'get_circle_nodes'

    if(.not. this%is_equispaced()) then
      call error_msg(proc_name, 'Not defined for non-equispaced grids (spacings: ' &
        // convert_to_char(this%axes%spacing), module_name_opt = module_name)
    end if
    radius_nodes = radius/this%axes(1)%spacing
    if(radius_nodes < 0.0_wp) then
      allocate(nodes(0))
      return
    else if(radius_nodes < 1.0_wp) then
      nodes = this%dim_to_linear_idx(p_center)
      return
    else if(all(radius_nodes > this%axes%n_nodes &
      .or. (this%axes%bound_cond%is_periodic .and. radius_nodes >= 0.5_wp*this%axes%n_nodes))) then
      nodes = [(i, i=1, this%n_nodes_total)]
      return
    end if

    select case(this%n_dims)
    case(1)
      n_nodes_upper_lim = ceiling(2*radius_nodes)
    case(2)
      n_nodes_upper_lim = ceiling(radius_nodes**2*pi)
    case(3)
      n_nodes_upper_lim = ceiling(4*radius_nodes**3*pi/3)
    case default
      n_nodes_upper_lim = 0
      call error_msg(proc_name, 'Not defined for the given amount of dimensions (' // convert_to_char(this%n_dims) // ')')
    end select
    allocate(nodes(n_nodes_upper_lim))

    where(radius_nodes < 0.5_wp*this%axes%n_nodes)
      p_low = mod_idx(p_center - ceiling(radius_nodes), this%axes%n_nodes)
      p_high = mod_idx(p_center + ceiling(radius_nodes), this%axes%n_nodes)
    elsewhere
      p_low = 1
      p_high = this%axes%n_nodes
    end where
    p = p_low
    n_nodes = 0
    do
      if(this%point_distance(p, p_center) <= radius_nodes) then
        n_nodes = n_nodes + 1
        if(n_nodes > n_nodes_upper_lim) then
          call error_msg(proc_name,&
            'Number of nodes (n_nodes) is greater than the assumed upper limit (n_nodes_upper_lim)')
        end if
        nodes(n_nodes) = this%dim_to_linear_idx(p)
      end if
      if(all(p == p_high)) return
      do dim = 1, this%n_dims
        if(p(dim) /= p_high(dim)) then
          p(dim) = mod_idx(p(dim) + 1, this%axes(dim)%n_nodes)
          exit
        else
          p(dim) = p_low(dim)
        end if
      end do
    end do
    if(n_nodes /= n_nodes_upper_lim) then
      ! Resize to actual dimensions
      allocate(nodes_temp(n_nodes))
      nodes_temp = nodes(1:n_nodes)
      deallocate(nodes)
      call move_alloc(nodes_temp, nodes)
    end if
  end function get_circle_nodes

  pure function get_circle_mask(this, p_center, radius) result(mask)
  class(regular_grid), intent(in) :: this
    integer, intent(in), dimension(:) :: p_center
    real(wp), intent(in) :: radius
    real(wp) :: radius_nodes
    logical, dimension(this%n_nodes_total) :: mask
    integer, dimension(this%n_dims) :: p_low, p_high, p
    integer :: dim
    character(*), parameter :: proc_name = 'get_circle_mask'

    if(.not. this%is_equispaced()) then
      call error_msg(proc_name, 'Not defined for non-equispaced grids (spacings: ' &
        // convert_to_char(this%axes%spacing), module_name_opt = module_name)
    end if
    radius_nodes = radius/this%axes(1)%spacing
    mask = .false.
    if(radius_nodes < 0.0_wp) return
    mask(this%dim_to_linear_idx(p_center)) = .true.
    if(radius_nodes < 1.0_wp) return
    if(all(radius_nodes > this%axes%n_nodes &
      .or. (this%axes%bound_cond%is_periodic .and. radius_nodes >= 0.5_wp*this%axes%n_nodes))) then
      mask = .true.
      return
    end if

    where(radius_nodes < 0.5_wp*this%axes%n_nodes)
      p_low = mod_idx(p_center - ceiling(radius_nodes), this%axes%n_nodes)
      p_high = mod_idx(p_center + ceiling(radius_nodes), this%axes%n_nodes)
    elsewhere
      p_low = 1
      p_high = this%axes%n_nodes
    end where
    p = p_low
    do
      if(this%point_distance(p, p_center) <= radius_nodes) then
        mask(this%dim_to_linear_idx(p)) = .true.
      end if
      if(all(p == p_high)) return
      do dim = 1, this%n_dims
        if(p(dim) /= p_high(dim)) then
          p(dim) = mod_idx(p(dim) + 1, this%axes(dim)%n_nodes)
          exit
        else
          p(dim) = p_low(dim)
        end if
      end do
    end do
  end function get_circle_mask

  pure subroutine increment_idx_full_range(this, p)
  class(regular_grid), intent(in) :: this
    integer, intent(inout) :: p(:)
    integer :: dim

    do dim = 1, this%n_dims
      if(p(dim) /= this%axes(dim)%n_nodes) then
        p(dim) = p(dim) + 1
        exit
      else
        p(dim) = 1
      end if
    end do
  end subroutine increment_idx_full_range

  pure subroutine increment_idx_limits(this, p_low, p_high, p)
  class(regular_grid), intent(in) :: this
    integer, intent(in), dimension(:) :: p_low, p_high
    integer, intent(inout) :: p(:)
    integer :: dim

    do dim = 1, this%n_dims
      if(p(dim) /= p_high(dim)) then
        p(dim) = mod_idx(p(dim) + 1, this%axes(dim)%n_nodes)
        exit
      else
        p(dim) = p_low(dim)
      end if
    end do
  end subroutine increment_idx_limits

  pure real(wp) function real_point_distance(this, p_1, p_2) result(dist)
  class(regular_grid), intent(in) :: this
    real(wp), intent(in), dimension(:) :: p_1, p_2
    real(wp), dimension(this%n_dims) :: dist_dim
    character(*), parameter :: proc_name = 'real_point_distance'

    if(any(p_1 < 0.0_wp .or. p_1 > this%axes%length)) then
      call error_msg(proc_name, 'Point "p_1" outside grid dimensions', module_name_opt = module_name)
    end if
    if(any(p_2 < 0.0_wp .or. p_2 > this%axes%length)) then
      call error_msg(proc_name, 'Point "p_2" outside grid dimensions', module_name_opt = module_name)
    end if
    dist_dim = abs(p_2 - p_1)
    where(this%axes%bound_cond%is_periodic)
      where(dist_dim > 0.5_wp*this%axes%length)
        dist_dim = this%axes%length - dist_dim
      end where
    end where
    dist = norm2(dist_dim)
  end function real_point_distance

  pure elemental integer recursive function offset_single_dim_idx(&
      dim_idx, offset, bc_type_min, bc_type_max, max_idx) result(offset_idx)
    integer, intent(in) :: dim_idx, offset, max_idx
    integer(int8), intent(in) :: bc_type_min, bc_type_max
    character(*), parameter :: proc_name = 'offset_dim_idx'

    offset_idx = dim_idx + offset
    if(offset_idx > max_idx) then
      select case(bc_type_max)
      case(bc_periodic)
        offset_idx = mod_idx(offset_idx, max_idx)
      case(bc_neumann)
        offset_idx = 2*max_idx - offset_idx
        offset_idx = offset_single_dim_idx(offset_idx, 0, bc_type_min, bc_type_max, max_idx)
      case(bc_dirichlet)
        offset_idx = max_idx
      case default
        call error_msg(proc_name, 'Not defined for the given upper boundary condition type = '&
          // convert_to_char(int(bc_type_max)), module_name_opt = module_name)
      end select
    else if(offset_idx < 1) then
      select case(bc_type_max)
      case(bc_periodic)
        offset_idx = mod_idx(offset_idx, max_idx)
      case(bc_neumann)
        offset_idx = -(offset_idx + 1)
        offset_idx = offset_single_dim_idx(offset_idx, 0, bc_type_min, bc_type_max, max_idx)
      case(bc_dirichlet)
        offset_idx = 1
      case default
        call error_msg(proc_name, 'Not defined for the given lower boundary condition type = '&
          // convert_to_char(int(bc_type_min)), module_name_opt = module_name)
      end select
    end if
  end function offset_single_dim_idx

  pure subroutine offset_dim_idx(this, dim_idx_start, offset, dim_idx_offset)
  class(regular_grid), intent(in) :: this
    integer, intent(in), dimension(:) :: dim_idx_start, offset
    integer, intent(out), dimension(:) :: dim_idx_offset

    dim_idx_offset = offset_single_dim_idx(&
      dim_idx_start, offset, this%axes%bound_cond%type_min, this%axes%bound_cond%type_max, this%axes%n_nodes)
  end subroutine offset_dim_idx

  pure function get_surface_points(this, mask) result(surf_p)
  class(regular_grid), intent(in) :: this
    logical, intent(in) :: mask(:)
    type(int_vector), allocatable, dimension(:) :: surf_p, surf_p_temp
    integer :: node, n_surf, dim, offset_val, offset_node
    integer, allocatable, dimension(:) :: p, offset, p_offset
    character(*), parameter :: proc_name = 'get_surface_points'

    if(size(mask) /= this%n_nodes_total) call error_msg(proc_name, 'Size of "mask" inconsistent with number of grid nodes')
    if(all(mask) .or. .not. any(mask)) then
      allocate(surf_p(0))
      return
    end if

    n_surf = 0
    allocate(surf_p(count(mask)))
    allocate(p(this%n_dims))
    allocate(p_offset(this%n_dims))
    allocate(offset(this%n_dims))
    do concurrent(node = 1:this%n_nodes_total, mask(node))
      p = this%linear_to_dim_idx(node)
      outer: do dim = 1, this%n_dims
        offset = 0
        do offset_val = -1, 1, 2
          offset(dim) = offset_val
          call this%offset_dim_idx(p, offset, p_offset)
          offset_node = this%dim_to_linear_idx(p_offset)
          if(.not. mask(offset_node)) then
            n_surf = n_surf + 1
            surf_p(n_surf)%v = p
            exit outer
          end if
        end do
      end do outer
    end do
    ! Resize to actual dimensions
    allocate(surf_p_temp(n_surf))
    surf_p_temp = surf_p(1:n_surf)
    deallocate(surf_p)
    call move_alloc(surf_p_temp, surf_p)
  end function get_surface_points
end module mod_regular_grid
