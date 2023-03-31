! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_piecewise_linear
  use mod_base_functions
  implicit none
  private

  character(*), parameter :: module_name = 'mod_piecewise_linear'

  type, public :: piecewise_linear
    private
    integer :: n_val
    real(wp), allocatable, dimension(:) :: x_known, y_known, dy_dx
    logical :: extrapolate = .false.
  contains
    procedure, public :: init
    procedure, public :: get_value_at
    procedure, public :: get_slope_at
    procedure, public :: get_low_idx
    procedure, public :: get_x_known
    procedure, public :: get_n_val
  end type piecewise_linear

contains

  pure subroutine init(this, x_known, y_known, extrapolate_opt)
  class(piecewise_linear), intent(inout) :: this
    real(wp), intent(in), dimension(:) :: x_known, y_known
    logical, intent(in), optional :: extrapolate_opt
    character(*), parameter :: proc_name = 'init'

    if(size(x_known) /= size(y_known)) then
      call error_msg(proc_name, 'Size of x_known and y_known is different', module_name_opt = module_name)
    else if(size(x_known) <= 1) then
      call error_msg(proc_name, 'Size of x_known/y_known has to be greater than 1', module_name_opt = module_name)
    else if(any(x_known(2:) <= x_known(:size(x_known) - 1))) then
      call error_msg(proc_name, 'x_known given not monotonically increasing', module_name_opt = module_name)
    end if
    if(present(extrapolate_opt)) this%extrapolate = extrapolate_opt
    this%x_known = x_known
    this%y_known = y_known
    this%n_val = size(x_known)
    this%dy_dx = (y_known(2:) - y_known(:size(y_known) - 1))/(x_known(2:) - x_known(:size(x_known) - 1))
  end subroutine init

  pure elemental real(wp) function get_value_at(this, x) result(y)
  class(piecewise_linear), intent(in) :: this
    real(wp), intent(in) :: x
    integer :: low_idx
    character(*), parameter :: proc_name = 'get_value_at'

    low_idx = this%get_low_idx(x)
    if(low_idx == 0) then
      y = this%y_known(1) + (x - this%x_known(1))*this%dy_dx(1)
    else if(low_idx == this%n_val) then
      y = this%y_known(this%n_val) + (x - this%x_known(this%n_val))*this%dy_dx(this%n_val - 1)
    else
      y = this%y_known(low_idx) + (x - this%x_known(low_idx))*this%dy_dx(low_idx)
    end if
  end function get_value_at

  pure elemental real(wp) function get_slope_at(this, x) result(y)
  class(piecewise_linear), intent(in) :: this
    real(wp), intent(in) :: x
    integer :: low_idx
    character(*), parameter :: proc_name = 'get_slope_at'

    low_idx = this%get_low_idx(x)
    if(low_idx == 0) then
      y = this%dy_dx(1)
    else if(low_idx == this%n_val) then
      y = this%dy_dx(this%n_val - 1)
    else
      y = this%dy_dx(low_idx)
    end if
  end function get_slope_at

  pure elemental integer function get_low_idx(this, x) result(low_idx)
  class(piecewise_linear), intent(in) :: this
    real(wp), intent(in) :: x
    character(*), parameter :: proc_name = 'get_low_idx'

    low_idx = last_true(x >= this%x_known)
    if(low_idx == 0 .or. low_idx >= this%n_val) then
      if(.not. this%extrapolate) then
        call error_msg(proc_name,&
          'x given (' // convert_to_char(x) // ') outside interpolation range (' &
          // convert_to_char(this%x_known(1)) // ' - ' // convert_to_char(this%x_known(this%n_val)) &
          // ') (extrapolation not enabled)', module_name_opt = module_name)
      end if
    end if
  end function get_low_idx

  pure function get_x_known(this) result(x_known)
  class(piecewise_linear), intent(in) :: this
    real(wp), dimension(size(this%x_known)) :: x_known

    x_known = this%x_known
  end function get_x_known

  pure integer function get_n_val(this) result(n_val)
  class(piecewise_linear), intent(in) :: this

    n_val = this%n_val
  end function get_n_val
end module mod_piecewise_linear
