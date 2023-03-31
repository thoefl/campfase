! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_solid_shape
  use mod_base_functions
  use mod_regular_grid
  implicit none
  private

  character(*), parameter :: module_name = 'mod_solid_shape'

  integer(int8), parameter, public :: shape_infinite = 0
  integer(int8), parameter, public :: shape_plane = 1
  integer(int8), parameter, public :: shape_round = 2

  type :: solid_shape_axis
    integer(int8) :: shape
    real(wp) :: length
  end type solid_shape_axis

  type, public :: solid_shape
    type(solid_shape_axis), allocatable :: axes(:)
  contains
    procedure, public :: get_solid_mask
  end type solid_shape

contains

  pure function get_solid_mask(this, grid) result(solid)
  class(solid_shape), intent(in) :: this
    type(regular_grid), intent(in) :: grid
    integer :: node, n_round_dims
    integer, dimension(size(this%axes)) :: p, p_low, p_high, p_box_low, p_box_high, p_dim, p_dim_low, p_dim_high
    real(wp), dimension(size(this%axes)) :: p_center, n_empty_nodes
    real(wp) :: radius
    logical, dimension(grid%n_nodes_total) :: solid
    logical, dimension(size(this%axes)) :: round_dims
    character(*), parameter :: proc_name = 'get_solid_mask'

    if(any(this%axes%length <= 0)) then
      call error_msg(proc_name, 'Shape axis length has to be greater than 0', module_name_opt = module_name)
    end if
    solid = .false.

    round_dims = (this%axes%shape == shape_round)
    n_round_dims = count(round_dims)
    p_center = (grid%axes%length + 1)/2.0_wp

    where(this%axes%shape == shape_infinite)
      n_empty_nodes = 0
      p_box_low = 1
      p_box_high = grid%axes%n_nodes
    elsewhere
      n_empty_nodes = max(grid%axes%n_nodes - this%axes%length, 0.0_wp)
      p_box_low = 1 + ceiling(n_empty_nodes/2)
      p_box_high = grid%axes%n_nodes - ceiling(n_empty_nodes/2)
    end where

    if(n_round_dims == 1) then
      call error_msg(proc_name, 'Number of round dimensions cannot be exactly 1', module_name_opt = module_name)
    else if(n_round_dims < 0) then
      call error_msg(proc_name, 'Number of round dimensions cannot be negative', module_name_opt = module_name)
    else if(n_round_dims == 0) then
      p = p_box_low
      do
        node = grid%dim_to_linear_idx(p)
        solid(node) = .true.
        if(all(p == p_box_high)) exit
        call grid%increment_idx(p_box_low, p_box_high, p)
      end do
    else if(n_round_dims > 1) then
      if(.not. all_equal(pack(this%axes%length, round_dims))) then
        call error_msg(proc_name,&
          'Currently only supports circular round shapes, i.e. length has to be equal for all round dimensions', &
          module_name_opt = module_name)
      end if
      radius = minval(this%axes%length, round_dims)/2.0_wp
      if(radius <= 0.0_wp) call error_msg(proc_name, 'Radius <= 0')

      p_low = max(floor(p_center - radius), 1)
      p_high = min(ceiling(p_center + radius), grid%axes%n_nodes)

      where(.not. round_dims)
        p_low = p_box_low
        p_high = p_box_low
        p_center = p_box_low
        p_dim_high = p_box_high
      end where

      p = p_low

      if(.not. grid%is_equispaced()) then
        call error_msg(proc_name, 'Currently only supports equispaced grids', module_name_opt = module_name)
      end if
      do
        if(grid%point_distance(p*grid%axes%spacing, p_center) <= radius) then
          ! If all dimensions are round, set just one point
          if(n_round_dims == grid%n_dims) then
            node = grid%dim_to_linear_idx(p)
            solid(node) = .true.
            ! Else, replicate in non-round dimensions
          else
            p_dim_low = p
            where(round_dims)
              p_dim_high = p
            end where
            p_dim = p_dim_low
            do
              node = grid%dim_to_linear_idx(p_dim)
              solid(node) = .true.
              if(all(p_dim == p_dim_high)) exit
              call grid%increment_idx(p_dim_low, p_dim_high, p_dim)
            end do
          end if
        end if
        if(all(p == p_high)) exit
        call grid%increment_idx(p_low, p_high, p)
      end do
    end if
  end function get_solid_mask
end module mod_solid_shape
