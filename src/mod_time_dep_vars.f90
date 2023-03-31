! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_time_dep_vars
  use mod_base_functions
  use mod_regular_grid
  implicit none
  private

  character(*), parameter :: module_name = 'mod_time_dep_vars'
  real(wp), parameter :: time_uninitialized = -huge(0.0_wp)
  real(wp), parameter :: time_invalid = -1.0_wp

  type, private :: time_dependent_variable
    real(wp) :: time = time_uninitialized
    character(:), allocatable :: var_name
    type(regular_grid) :: grid
  contains
    procedure, public :: time_is_equal
    procedure, public :: time_is_greater_equal
    procedure, public :: time_diff_within_tolerance
    procedure, public :: check_is_valid
    procedure, public :: check_time_is_equal
    procedure, public :: check_time_is_greater_equal
    procedure, public :: check_time_diff_within_tolerance
    procedure, public :: init
    procedure, public :: output_image => var_output_image
    procedure, public :: invalidate
  end type time_dependent_variable

  type, public, extends(time_dependent_variable) :: time_dependent_integer
    !dir$ attributes align : 64 :: vals
    integer, allocatable, dimension(:) :: vals
  end type time_dependent_integer

  type, public, extends(time_dependent_variable) :: time_dependent_real
    !dir$ attributes align : 64 :: vals
    real(wp), allocatable, dimension(:) :: vals
  contains
    procedure, public :: coordinate_shift
  end type time_dependent_real

  type, public, extends(time_dependent_variable) :: time_dependent_complex
    !dir$ attributes align : 64 :: vals
    complex(wp), allocatable, dimension(:) :: vals
  end type time_dependent_complex

contains

  pure subroutine init(this, grid, var_name)
  class(time_dependent_variable), intent(inout) :: this
    type(regular_grid), intent(in), target :: grid
    character(*), intent(in) :: var_name
    character(*), parameter :: proc_name = 'init'

    this%var_name = trim(var_name)
    this%grid = grid
    this%time = time_invalid
    select type(this)
    class is(time_dependent_real)
      allocate(this%vals(grid%n_nodes_total))
      this%vals = 0
    class is(time_dependent_integer)
      allocate(this%vals(grid%n_nodes_total))
      this%vals = 0
    class is(time_dependent_complex)
      allocate(this%vals(grid%n_nodes_total))
      this%vals = 0
    class default
      call error_msg(proc_name, 'Not defined for the given type', module_name_opt=module_name)
    end select
  end subroutine init

  pure elemental logical function time_is_equal(a, b)
  class(time_dependent_variable), intent(in) :: a, b

    time_is_equal = a%time == b%time .and. a%time >= 0.0_wp
  end function time_is_equal

  pure elemental logical function time_is_greater_equal(a, b)
  class(time_dependent_variable), intent(in) :: a, b

    time_is_greater_equal = a%time >= b%time
  end function time_is_greater_equal

  pure elemental logical function time_diff_within_tolerance(a, b, time_diff)
  class(time_dependent_variable), intent(in) :: a, b
    real(wp), intent(in) :: time_diff
    real(wp) :: time_tol

    time_tol = 1E6*max(spacing(a%time), spacing(b%time), spacing(time_diff))        
    time_diff_within_tolerance = abs(a%time - b%time - time_diff) <= time_tol
  end function time_diff_within_tolerance

  pure elemental subroutine check_is_valid(this, caller_name)
  class(time_dependent_variable), intent(in) :: this
    character(*), intent(in) :: caller_name
    character(:), allocatable :: err_str

    if(this%time < 0) then
      err_str = trim(img_head()) // ' ERROR: ' // caller_name // ': '
      if(this%time == time_uninitialized) then
        err_str = err_str // 'variable ' // this%var_name // ' has not been initialized!'
      else if(this%time == time_invalid) then
        err_str = err_str // 'variable ' // this%var_name // ' has been set as invalid!'
      else
        err_str = err_str // 'unknown error (time for variable ' // this%var_name // ' is negative)'
      end if
      error stop err_str
    end if
  end subroutine check_is_valid

  pure elemental subroutine check_time_is_equal(a, b, caller_name)
  class(time_dependent_variable), intent(in) :: a, b
    character(*), intent(in) :: caller_name
    character(:), allocatable :: err_str

    if(.not. a%time_is_equal(b)) then
      err_str = trim(img_head()) // ' ERROR: ' // caller_name // ': variable times differ: ' //&
        a%var_name // ' = ' // convert_to_char(a%time) // ' /= ' //&
        b%var_name // ' = ' // convert_to_char(b%time)
      error stop err_str
    end if
  end subroutine check_time_is_equal

  pure elemental subroutine check_time_is_greater_equal(a, b, caller_name)
  class(time_dependent_variable), intent(in) :: a, b
    character(*), intent(in) :: caller_name
    character(:), allocatable :: err_str

    if(a%time < b%time) then
      err_str = trim(img_head()) // ' ERROR: ' // caller_name // ': variable time smaller: ' //&
        a%var_name // ' = ' // convert_to_char(a%time) // ' < ' //&
        b%var_name // ' = ' // convert_to_char(b%time)
      error stop err_str
    end if
  end subroutine check_time_is_greater_equal

  pure elemental subroutine check_time_diff_within_tolerance(a, b, time_diff, caller_name)
  class(time_dependent_variable), intent(in) :: a, b
    real(wp), intent(in) :: time_diff
    character(*), intent(in) :: caller_name
    character(:), allocatable :: err_str

    if(.not. a%time_diff_within_tolerance(b, time_diff)) then
      err_str = trim(img_head()) // ' ERROR: ' // caller_name // ': variable time difference greater than tolerance: ' //&
        a%var_name // ' - ' // b%var_name // ' - dt = ' // convert_to_char(a%time - b%time - time_diff) // ', tol: ' //&
        convert_to_char(1E6*max(spacing(a%time), spacing(b%time), spacing(time_diff)))
      error stop err_str
    end if
  end subroutine check_time_diff_within_tolerance

  subroutine var_output_image(this, name_suffix_opt, a_min_opt, a_max_opt, co_img_opt)
  class(time_dependent_variable), intent(in) :: this
    character(*), intent(in), optional :: name_suffix_opt
    real(wp), intent(in), optional :: a_min_opt, a_max_opt
    logical, intent(in), optional :: co_img_opt
    character(*), parameter :: real_suffix = '_real'
    character(*), parameter :: im_suffix = '_imaginary'
    character(*), parameter :: proc_name = 'output_image'
    character(:), allocatable :: name

    call this%check_is_valid(proc_name)
    name = this%var_name
    if(present(name_suffix_opt)) then
      if(len(name_suffix_opt) > 0) then
        name = name // '_' // name_suffix_opt
      end if
    end if
    select type(this)
    class is(time_dependent_real)
      call output_image(this%vals, name, this%grid%axes(1)%n_nodes, a_min_opt, a_max_opt, co_img_opt)
    class is(time_dependent_integer)
      ! FIXME: allow choice of a_min/a_max
      call output_image(this%vals, name, this%grid%axes(1)%n_nodes, co_img_opt=co_img_opt)
    class is(time_dependent_complex)
      call output_image(this%vals%re, name // real_suffix, this%grid%axes(1)%n_nodes, a_min_opt, a_max_opt, co_img_opt)
      call output_image(this%vals%im, name // im_suffix, this%grid%axes(1)%n_nodes, a_min_opt, a_max_opt, co_img_opt)
    end select
  end subroutine var_output_image

  pure subroutine invalidate(this)
  class(time_dependent_variable), intent(inout) :: this

    this%time = time_invalid
  end subroutine invalidate

  pure subroutine coordinate_shift(this, shift)
  class(time_dependent_real), intent(inout) :: this
    real(wp), intent(in) :: shift(:)
    real(wp), allocatable :: original_vals(:)
    integer :: dim, linear_idx, linear_idx_offset
    real(wp) :: shift_factor
    integer, dimension(this%grid%n_dims) :: offset, dim_idx_offset, dim_idx
    character(*), parameter :: proc_name = 'coord_shift'

    if(any(shift < 0 .or. shift > 1)) then
      call error_msg(proc_name, 'Shift value has to be between 0 and 1', module_name_opt = module_name)
    end if
    if(size(shift) /= this%grid%n_dims) then
      call error_msg(proc_name, &
        'Shift dimensions different from variable dimensions: ' // this%var_name, module_name_opt = module_name)
    end if
    allocate(original_vals, source = this%vals)
    dim_idx = 1
    offset = 0
    this%vals = 0
    do
      shift_factor = product(merge(shift, 1 - shift, offset == 0))
      do concurrent(linear_idx = 1:this%grid%n_nodes_total)
        dim_idx = this%grid%linear_to_dim_idx(linear_idx)
        call this%grid%offset_dim_idx(dim_idx, offset, dim_idx_offset)
        linear_idx_offset = this%grid%dim_to_linear_idx(dim_idx_offset)
        this%vals(linear_idx) = this%vals(linear_idx) + shift_factor*original_vals(linear_idx_offset)
      end do
      if(all(offset == 1)) exit
      do dim = 1, this%grid%n_dims
        if(offset(dim) == 0) then
          offset(dim) = 1
          exit
        else if(offset(dim) == 1) then
          offset(dim) = 0
        else
          call error_msg(proc_name, 'Offset should be either 0 or 1')
        end if
      end do
    end do
    deallocate(original_vals)
  end subroutine coordinate_shift
end module mod_time_dep_vars
