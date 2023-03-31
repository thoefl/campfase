! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_rational_number
  use iso_fortran_env, only: int64
  use mod_base_functions, only: operator(//)
  implicit none
  private

  type, public :: rational_number
    integer(int64) :: dividend, divisor
  contains
    procedure, public :: simplify
  end type rational_number

  public operator (+)
  public operator (*)
  public operator (-)
  public operator (/)

  interface operator (+)
    module procedure add_int
  end interface

  interface operator (*)
    module procedure multiply_int
  end interface

  interface operator (-)
    module procedure subtract_int
  end interface

  interface operator (/)
    module procedure divide_int
  end interface

contains

  pure elemental subroutine simplify(this)
  class(rational_number), intent(inout) :: this
    integer(int64) :: factor

    factor = greatest_common_divisor(this%dividend, this%divisor)
    this%dividend = this%dividend//factor
    this%divisor = this%divisor//factor
  end subroutine simplify

  pure elemental type(rational_number) function add_int(frac, i) result(res)
    type(rational_number), intent(in) :: frac
    integer, intent(in) :: i

    res%dividend = frac%dividend + frac%divisor*i
    res%divisor = frac%divisor
    call res%simplify()
  end function add_int

  pure elemental type(rational_number) function subtract_int(frac, i) result(res)
    type(rational_number), intent(in) :: frac
    integer, intent(in) :: i

    res = frac + (-i)
  end function subtract_int

  pure elemental type(rational_number) function multiply_int(frac, i) result(res)
    type(rational_number), intent(in) :: frac
    integer, intent(in) :: i

    res%dividend = frac%dividend*i
    res%divisor = frac%divisor
    call res%simplify()
  end function multiply_int

  pure elemental type(rational_number) function divide_int(frac, i) result(res)
    type(rational_number), intent(in) :: frac
    integer, intent(in) :: i

    res%dividend = frac%dividend
    res%divisor = frac%divisor*i
    call res%simplify()
  end function divide_int

  pure elemental integer(int64) recursive function greatest_common_divisor(a, b) result(div)
    integer(int64), intent(in) :: a, b
    integer(int64) :: remainder

    if(a == b) then
      div = a
      return
    else if(b > a) then
      div = greatest_common_divisor(b, a)
    else
      remainder = mod(abs(a), abs(b))
      if(remainder == 0) then
        div = b
      else
        div = greatest_common_divisor(b, remainder)
      end if
    end if
  end function greatest_common_divisor
end module mod_rational_number
