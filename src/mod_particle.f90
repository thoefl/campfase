! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_particle
  use mod_base_functions
  use mod_regular_grid
  implicit none
  private

  type, public :: spheric_particle
    type(regular_grid), pointer :: grid
    real(wp) :: radius
    integer, dimension(:), allocatable :: coordinates
    integer :: linear_idx, id
  contains
    procedure, public :: init
    procedure, public :: get_surface_area
    procedure, public :: get_volume
  end type spheric_particle

contains

  subroutine init(this, grid, id, radius, coordinates)
    type(regular_grid), intent(in), target :: grid
  class(spheric_particle), intent(inout) :: this
    integer, intent(in) :: id, coordinates(:)
    real(wp), intent(in) :: radius

    if(radius < 0) error stop 'ERROR: particle_init: Radius cannot be negative!'
    if(any(coordinates < 0)) error stop 'ERROR: particle_init: Coordinates cannot be negative!'

    this%grid => grid
    this%id = id
    this%radius = radius
    this%coordinates = coordinates
    this%linear_idx = grid%dim_to_linear_idx(coordinates)
  end subroutine init

  pure elemental real(wp) function get_surface_area(this) result(surface_area)
  class(spheric_particle), intent(in) :: this
    character(*), parameter :: proc_name = 'get_surface_area'

    select case(this%grid%n_dims)
    case(1)
      surface_area = 2
    case(2)
      surface_area = 2*this%radius*pi
    case(3)
      surface_area = 4*this%radius**2*pi
    case default
      surface_area = 0
      call error_msg(proc_name, 'Not defined for given amount of spatial dimensions: ' // convert_to_char(this%grid%n_dims))
    end select
  end function get_surface_area

  pure elemental real(wp) function get_volume(this) result(volume)
  class(spheric_particle), intent(in) :: this
    character(*), parameter :: proc_name = 'get_volume'

    select case(this%grid%n_dims)
    case(1)
      volume = 2*this%radius
    case(2)
      volume = this%radius**2*pi
    case(3)
      volume = 4/3.0_wp*this%radius**3*pi
    case default
      volume = 0
      call error_msg(proc_name, 'Not defined for given amount of spatial dimensions: ' // convert_to_char(this%grid%n_dims))
    end select
  end function get_volume

end module mod_particle
