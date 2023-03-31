! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_base_functions
  use, intrinsic :: iso_fortran_env, only:  int8, int16, int32, int64, real32, real64, input_unit, output_unit, error_unit, &
    event_type, team_type, lock_type
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_positive_inf, ieee_negative_inf, ieee_value
  use mod_realtype
  use mod_array_index

  implicit none

  character(*), parameter :: default_character_format = 'a'
  character(*), parameter :: default_integer_format = 'i0'
  character(*), parameter :: default_logical_format = 'l1'
  character(*), parameter :: default_real_format = 'es12.5'
  character(*), parameter :: list_value_format = '(*(a," = ",g0.4,:,", "))'
  character(*), parameter :: list_format = '(*(g0.4,:," "))'
  character(*), parameter :: csv_format = '(*(g0.20,:,", "))'

  real(wp), parameter :: pi = 4.0_wp*ATAN(1.0_wp)

  real(wp), parameter :: delta_num = 1E6*epsilon(1.0_wp)

  integer, parameter, private :: mpi_cache_max = 128

  interface co_any
    module procedure :: co_any_scalar
    module procedure :: co_any_array
  end interface co_any

  interface co_all
    module procedure :: co_all_scalar
    module procedure :: co_all_array
  end interface co_all

  interface co_maxloc
    module procedure :: co_maxloc_scalar
    module procedure :: co_maxloc_array
  end interface co_maxloc

  interface all_equal
    module procedure :: all_equal_real
    module procedure :: all_equal_int
  end interface all_equal

  interface co_all_equal
    module procedure :: co_all_equal_real
    module procedure :: co_all_equal_int
    module procedure :: co_all_equal_int64
  end interface co_all_equal

  interface output_image
    module procedure :: output_real_image
    module procedure :: output_logical_image
    module procedure :: output_integer_image
  end interface output_image

  interface convert_to_char
    module procedure :: int_to_char
    module procedure :: int64_to_char
    module procedure :: real_to_char
    module procedure :: realarray_to_char
    module procedure :: realmatrix_to_char
    module procedure :: logical_to_char
    module procedure :: intarray_to_char
  end interface convert_to_char

  interface operator(//)
    module procedure pure_int_int_divide, pure_int64_int_divide, pure_int64_int64_divide
  end interface

  interface option_or_default
    module procedure :: option_or_default_real
    module procedure :: option_or_defalt_integer
    module procedure :: option_or_defalt_logical
    module procedure :: option_or_defalt_character
  end interface option_or_default

  interface cartesian_product
    module procedure :: cartesian_product_int
    module procedure :: cartesian_product_string
  end interface cartesian_product

  interface sign_is_equal
    module procedure :: sign_is_equal_real
  end interface sign_is_equal

  interface resize_array
    module procedure :: resize_real_array, resize_int_array
  end interface resize_array

  interface concatenate_image_arrays
    module procedure :: concatenate_image_int_arrays
  end interface concatenate_image_arrays

  type :: int_vector
    !dir$ attributes align : 64 :: v
    integer, dimension(:), allocatable :: v
  end type int_vector

  type :: int_ragged_matrix
    type(int_vector), allocatable :: row(:)
  end type int_ragged_matrix

  type :: real_vector
    !dir$ attributes align : 64 :: vals
    real(wp), dimension(:), allocatable :: vals
  end type real_vector

  type :: complex_vector
    !dir$ attributes align : 64 :: vals
    complex(wp), dimension(:), allocatable :: vals
  end type complex_vector

  type :: string
    character(:), allocatable :: c
  end type string

  type :: string_vector
    type(string), allocatable :: s(:)
  end type string_vector

contains

  pure elemental integer function mod_idx(val, nidx)
    integer, intent(in) :: val, nidx

    mod_idx = modulo(val - 1, nidx) + 1
  end function mod_idx

  pure integer function first_true(a)
    logical, intent(in) :: a(:)
    integer :: i

    first_true = 0
    do i = 1,size(a)
      if(a(i)) then
        first_true = i
        return
      end if
    end do
  end function first_true

  pure integer function last_true(a)
    logical, intent(in) :: a(:)
    integer :: i

    last_true = 0
    do i = size(a), 1, -1
      if(a(i)) then
        last_true = i
        return
      end if
    end do
  end function last_true

  pure integer function nth_true(a, n)
    logical, intent(in) :: a(:)
    integer, intent(in) :: n
    integer :: i, j

    nth_true = 0
    j = 0
    do i = 1, size(a)
      if(a(i)) then
        j = j + 1
        if(j == n) then
          nth_true = i
          return
        end if
      end if
    end do
  end function nth_true


  ! sleep for ms milliseconds using wall time
  subroutine millisleep(ms)
    integer, intent(in) :: ms
    integer(int64) :: t_wallcount, t_wallrate, t_wallcount_finish

    call system_clock(t_wallcount, t_wallrate)
    t_wallcount_finish = t_wallcount + ms*t_wallrate/1000
    do
      call system_clock(t_wallcount)
      if(t_wallcount >= t_wallcount_finish) exit
    end do
  end subroutine millisleep

  pure function img_head(step_opt)
    integer(int64), intent(in), optional :: step_opt
    character(128) :: img_head

    if(present(step_opt)) then
      write(img_head, '(a,i2,a,i2,a,i0,a)') 'Image ', this_image(), '/', num_images(), ', step ', step_opt, ': '
    else
      write(img_head, '(a,i2,a,i2,a)') 'Image ', this_image(), '/', num_images(), ': '
    end if
  end function img_head

  subroutine img_printhead(step_opt)
    integer(int64), intent(in), optional :: step_opt
    character(32) :: max_img, format_str
    integer :: max_img_chars

    write(max_img, '(i0)') num_images()
    max_img_chars = len_trim(max_img)
    write(format_str, '(a,i0,a)') '(2(a,i', max_img_chars, '))'
    write(output_unit, format_str, advance = 'no') 'Image ', this_image(), '/', num_images()
    if(present(step_opt)) then
      write(output_unit, '(a,i0)', advance = 'no') ', step ', step_opt
    end if
    write(output_unit, '(a)', advance = 'no') ': '
  end subroutine img_printhead

  pure subroutine string_insert(s, end_pos, s_in)
    character(*), intent(inout) :: s
    integer, intent(inout) :: end_pos
    character(*), intent(in) :: s_in
    integer :: s_len

    s_len = len(s_in)
    if(end_pos + s_len - 1 > len(s)) then
      error stop 'ERROR: string_insert: Maximum string length exceeded!'
    end if
    s(end_pos:end_pos + s_len - 1) = s_in
    call inc(end_pos, s_len)
  end subroutine string_insert

  subroutine newline()
    write (output_unit, '(a)', advance = 'yes') ''
  end subroutine newline

  pure subroutine inc(a, inc_opt)
    integer, intent(inout) :: a
    integer, intent(in), optional :: inc_opt

    if(present(inc_opt)) then
      a = a + inc_opt
    else
      a = a + 1
    end if
  end subroutine inc

  pure subroutine dec(a)
    integer, intent(inout) :: a

    a = a - 1
  end subroutine dec

  pure integer function count_unique(a)
    integer, intent(in) :: a(:)
    integer :: i, j

    count_unique = size(a)
    do i = 1, size(a)
      do j = i + 1, size(a)
        if(a(i) == a(j)) count_unique = count_unique - 1
      end do
    end do
  end function count_unique

  pure real(wp) function findroot_bisect(func, lbound_init, ubound_init, x_tol, y_target_opt, y_tol_opt) result(x)
    real(wp), intent(in) :: x_tol, lbound_init, ubound_init
    real(wp), intent(in), optional :: y_target_opt, y_tol_opt
    real(wp) :: y_tol
    integer, parameter :: max_iter = 100
    interface
      pure real(wp) function func(x)
        import :: wp
        real(wp), intent(in) :: x
      end function func
    end interface
    integer :: i
    real(wp) :: y, y_target, lbound, ubound, sign_y_ubound, sign_y_lbound, sign_y

    if(present(y_target_opt)) then
      y_target = y_target_opt
    else
      y_target = 0.0_wp
    end if
    if(present(y_tol_opt)) then
      y_tol = y_tol_opt
    else
      y_tol = 10*spacing(y_target)
    end if

    lbound = lbound_init
    ubound = ubound_init
    sign_y_ubound = sign(1.0_wp, func(ubound) - y_target)
    sign_y_lbound = sign(1.0_wp, func(lbound) - y_target)
    if(sign_y_ubound == sign_y_lbound) error stop 'ERROR: findroot_bisect: Upper and lower bound given have same result sign!'
    do i = 1, max_iter
      x = (ubound + lbound)/2
      y = func(x) - y_target
      if(abs(y) <= y_tol) return
      sign_y = sign(1.0_wp, y)
      if(sign_y == sign_y_ubound) then
        ubound = x
      else if(sign_y == sign_y_lbound) then
        lbound = x
      else
        error stop 'ERROR: findroot_bisect: Sign of y is not equal to both upper and lower bound sign!'
      end if
      if(abs(ubound - lbound) <= x_tol) return
    end do
    error stop 'ERROR: findroot_bisect failed to converge!'
  end function findroot_bisect

  pure real(wp) function findroot_newton(func, dfunc_dx, x_init, x_tol, lbound_init, ubound_init, y_target_opt) result(x)
    real(wp), intent(in) :: x_init, x_tol, lbound_init, ubound_init
    real(wp), intent(in), optional :: y_target_opt
    real(wp), parameter :: y_tol = epsilon(0.0_wp)
    integer, parameter :: max_iter = 100
    interface
      pure real(wp) function func(x)
        import :: wp
        real(wp), intent(in) :: x
      end function func
      pure real(wp) function dfunc_dx(x)
        import :: wp
        real(wp), intent(in) :: x
      end function dfunc_dx
    end interface
    integer :: i
    real(wp) :: y, dy_dx, x_next, y_target, lbound, ubound

    if(present(y_target_opt)) then
      y_target = y_target_opt
    else
      y_target = 0.0_wp
    end if

    x = x_init
    x_next = huge(x)
    lbound = lbound_init
    ubound = ubound_init
    do i = 1, max_iter
      dy_dx = dfunc_dx(x)
      if(abs(dy_dx) <= y_tol) exit
      y = func(x) - y_target
      x_next = x - y/dy_dx
      if(x_next > ubound) x_next = 0.5_wp*(x + ubound)
      if(x_next < lbound) x_next = 0.5_wp*(x + lbound)
      if(sign(1.0_wp, y) /= sign(1.0_wp, dy_dx)) then
        lbound = x
      else
        ubound = x
      end if
      if(abs(x - x_next) <= x_tol) exit
      x = x_next
    end do
    if(abs(x - x_next) > x_tol) then
      error stop 'ERROR: findroot_newton failed to converge!'
    end if
  end function findroot_newton

  pure function linspace(startval, stopval, nval)
    real(wp), intent(in) :: startval, stopval
    integer, intent(in) :: nval
    real(wp) :: linspace(nval)
    integer :: i

    linspace = [(startval + real(i,wp)/(nval - 1)*(stopval - startval), i=0, nval-1)]
  end function linspace

  pure logical function all_equal_real(val) result(all_equal)
    real(wp), intent(in) :: val(:)

    all_equal = all(val - val(1) == 0)
  end function all_equal_real

  pure logical function all_equal_int(val) result(all_equal)
    integer, intent(in) :: val(:)

    all_equal = all(val - val(1) == 0)
  end function all_equal_int

  pure integer function right_image()
    if(this_image() == num_images()) then
      right_image = 1
    else
      right_image = this_image() + 1
    end if
  end function right_image

  pure integer function left_image()
    if(this_image() == 1) then
      left_image = num_images()
    else
      left_image = this_image() - 1
    end if
  end function left_image

  pure integer function next_image(img)
    integer, intent(in) :: img
    if(img == num_images()) then
      next_image = 1
    else
      next_image = img + 1
    end if
  end function next_image

  pure integer function prev_image(img)
    integer, intent(in) :: img
    if(img == 1) then
      prev_image = num_images()
    else
      prev_image = img - 1
    end if
  end function prev_image

  ! Could be implemented using co_reduce and pairwise all/any instead, but MPI_REDUCE seems to have a memory leak bug (?) when not
  ! directly followed by a call to MPI_BARRIER, see:
  ! https://stackoverflow.com/questions/33754220/mpi-reduce-causing-memory-leak
  subroutine co_any_scalar(a)
    logical, intent(inout) :: a
    integer(int8) :: a_int
    
    a_int = merge(1_int8, 0_int8, a)
    call co_max(a_int)
    a = (a_int == 1_int8)
  end subroutine co_any_scalar

  subroutine co_any_array(a)
    logical, intent(inout), contiguous :: a(:)
    integer(int8), dimension(:), allocatable :: a_int
    integer :: i, remaining, cache_size
   
    i = 1
    do
      remaining = size(a) - i + 1
      if(remaining == 0) return
      cache_size = min(remaining, mpi_cache_max)
      if(.not. allocated(a_int)) then
        allocate(a_int(cache_size))
      else if(cache_size < mpi_cache_max) then
        deallocate(a_int)
        allocate(a_int(cache_size))
      end if
      a_int = merge(1_int8, 0_int8, a(i:i + cache_size - 1))
      call co_max(a_int)
      a(i:i + cache_size - 1) = (a_int == 1_int8)
      i = i + cache_size
    end do
  end subroutine co_any_array

  ! See comment for co_any
  subroutine co_all_scalar(a)
    logical, intent(inout) :: a
    integer(int8) :: a_int
    
    a_int = merge(1_int8, 0_int8, a)
    call co_min(a_int)
    a = (a_int == 1_int8)
  end subroutine co_all_scalar

  subroutine co_all_array(a)
    logical, intent(inout), contiguous :: a(:)
    integer(int8), dimension(:), allocatable :: a_int
    integer :: i, remaining, cache_size
   
    i = 1
    do
      remaining = size(a) - i + 1
      if(remaining == 0) return
      cache_size = min(remaining, mpi_cache_max)
      if(.not. allocated(a_int)) then
        allocate(a_int(cache_size))
      else if(cache_size < mpi_cache_max) then
        deallocate(a_int)
        allocate(a_int(cache_size))
      end if
      a_int = merge(1_int8, 0_int8, a(i:i + cache_size - 1))
      call co_min(a_int)
      a(i:i + cache_size - 1) = (a_int == 1_int8)
      i = i + cache_size
    end do
  end subroutine co_all_array

  integer function co_maxloc_scalar(a) result(co_maxloc)
    real(wp), intent(in) :: a
    real(wp) :: a_max

    a_max = a
    call co_max(a_max)
    if(a == a_max) then
      co_maxloc = this_image()
    else
      co_maxloc = num_images() + 1
    end if
    call co_min(co_maxloc)
  end function co_maxloc_scalar

  function co_maxloc_array(a) result(co_maxloc)
    real(wp), intent(in), contiguous :: a(:)
    real(wp), dimension(:), allocatable :: a_max
    integer, dimension(size(a)) :: co_maxloc
    integer :: i, remaining, cache_size
   
    i = 1
    do
      remaining = size(a) - i + 1
      if(remaining == 0) exit
      cache_size = min(remaining, mpi_cache_max)
      if(.not. allocated(a_max)) then
        allocate(a_max(cache_size))
      else if(cache_size < mpi_cache_max) then
        deallocate(a_max)
        allocate(a_max(cache_size))
      end if
      a_max = a(i:i + cache_size - 1)
      call co_max(a_max)
      co_maxloc(i:i + cache_size - 1) = merge(this_image(), num_images() + 1, a(i:i + cache_size - 1) == a_max)
      i = i + cache_size
    end do
    call co_min(co_maxloc)
  end function co_maxloc_array

  pure subroutine shift_matrix(a, nx, shift)
    real(wp), intent(inout), dimension(:) :: a
    integer, intent(in) :: nx, shift
    integer :: i, ny, x_start, x_end

    ny = size(a)/nx
    do concurrent(i = 1:ny)
      x_start = (i - 1)*nx + 1
      x_end = i*nx
      a(x_start:x_end) = cshift(a(x_start:x_end), shift)
    end do
  end subroutine shift_matrix

  pure real(wp) function mean(a, mask_opt)
    real(wp), intent(in) :: a(:)
    logical, intent(in), optional :: mask_opt(:)
    logical :: mask(size(a))

    if(present(mask_opt)) then
      mask = mask_opt
    else
      mask = .true.
    end if

    if(any(mask)) then
      mean = sum(a, mask)/count(mask)
    else
      mean = 0.0_wp
    end if
  end function mean

  pure elemental real(wp) function harmonic(n)
    integer, intent(in) :: n
    integer :: i

    harmonic = 1
    ! Start with the smallest values (largest n) to preserve accuracy
    do i = n, 2, -1
      harmonic = harmonic + 1.0_wp/i
    end do
  end function harmonic

  pure elemental integer function pure_int_int_divide(i, j) result(div)
    integer, intent(in) :: i, j
    character(*), parameter :: proc_name = 'pure_int_int_divide'

    if(mod(i,j) /= 0) then
      call error_msg(proc_name, 'Pure integer division impossible: ' // convert_to_char(i) // '/' // convert_to_char(j))
    end if
    div = i/j
  end function pure_int_int_divide

  pure elemental integer(int64) function pure_int64_int_divide(i, j) result(div)
    integer(int64), intent(in) :: i
    integer, intent(in) :: j
    character(*), parameter :: proc_name = 'pure_int64_int_divide'

    if(mod(i, int(j, int64)) /= 0) then
      call error_msg(proc_name, 'Pure integer division impossible: ' // convert_to_char(i) // '/' // convert_to_char(j))
    end if
    div = i/j
  end function pure_int64_int_divide

  pure elemental integer(int64) function pure_int64_int64_divide(i, j) result(div)
    integer(int64), intent(in) :: i, j
    character(*), parameter :: proc_name = 'pure_int64_int64_divide'

    if(mod(i,j) /= 0) then
      call error_msg(proc_name, 'Pure integer division impossible: ' // convert_to_char(i) // '/' // convert_to_char(j))
    end if
    div = i/j
  end function pure_int64_int64_divide

  pure real(wp) function norm2_error(approx, exact)
    real(wp), intent(in) :: approx(:), exact(:)

    if(any(abs(exact) > delta_num)) then
      norm2_error = norm2(exact - approx)/norm2(exact)
    else
      norm2_error = 0.0_wp
    end if
  end function norm2_error

  pure elemental real(wp) function rel_error(a, b)
    real(wp), intent(in) :: a, b

    if(a /= 0 .or. b /= 0) then
      rel_error = abs(a - b)/max(abs(a), abs(b))
    else
      rel_error = 0.0_wp
    end if
  end function rel_error

  pure elemental real(wp) function abs_error(a, b)
    real(wp), intent(in) :: a, b

    abs_error = abs(a - b)
  end function abs_error

  pure logical function is_square_matrix(a)
    real(wp), intent(in) :: a(:,:)

    is_square_matrix = size(a,1) == size(a,2)
  end function is_square_matrix

  pure function matrix_diag(a)
    real(wp), intent(in) :: a(:,:)
    real(wp) :: matrix_diag(size(a,1))
    integer :: i

    if(.not. is_square_matrix(a)) then
      error stop 'ERROR: matrix_diag: input is not a square 2d matrix!'
    end if
    do concurrent(i = 1:size(a,1))
      matrix_diag(i) = a(i,i)
    end do
  end function matrix_diag

  pure integer function character_count(string, char)
    character(*), intent(in) :: string
    character, intent(in) :: char
    integer :: i

    character_count = 0
    do concurrent(i = 1:len(string))
      if(string(i:i) == char) call inc(character_count)
    end do
  end function character_count

  pure logical function string_contains_character(string, char)
    character(*), intent(in) :: string
    character, intent(in) :: char
    integer :: i

    string_contains_character = .false.
    do i = 1, len(string)
      if(string(i:i) == char) then
        string_contains_character = .true.
        return
      end if
    end do
  end function string_contains_character

  pure subroutine split_string(str, delimiter, split_str)
    character(*), intent(in) :: str
    character, intent(in) :: delimiter
    type(string), intent(out), allocatable :: split_str(:)
    integer :: n_delim, n_chars, i, delim_count, split_start_idx
    integer, dimension(:), allocatable :: delim_idx

    if(delimiter /= ' ') then
      n_chars = len_trim(str)
    else
      n_chars = len(str)
    end if
    n_delim = character_count(str, delimiter)
    if(n_delim == 0) then
      allocate(split_str(1))
      split_str(1)%c = str
      return
    end if
    allocate(delim_idx(n_delim))
    delim_count = 0
    do i = 1, n_chars
      if(str(i:i) == delimiter) then
        call inc(delim_count)
        delim_idx(delim_count) = i
      end if
    end do
    if(delim_count /= n_delim) then
      error stop 'ERROR: split_str: Number of delimiters found does not match delimiter character count!'
    end if
    allocate(split_str(n_delim + 1))
    split_start_idx = 1
    do i = 1, n_delim
      if(split_start_idx >= delim_idx(i)) then
        split_str(i)%c = ''
      else
        split_str(i)%c = str(split_start_idx:delim_idx(i) - 1)
      end if
      split_start_idx = delim_idx(i) + 1
    end do
    if(split_start_idx >= n_chars) then
      split_str(n_delim + 1)%c = ''
    else
      split_str(n_delim + 1)%c = str(split_start_idx:n_chars)
    end if
  end subroutine split_string

  pure function filename_without_extension(filename)
    character(*), intent(in) :: filename
    character(:), allocatable :: filename_without_extension
    type(string), allocatable :: split_filename(:)
    integer :: n_parts

    call split_string(filename, '.', split_filename)
    n_parts = size(split_filename)
    if(n_parts <= 2) then
      filename_without_extension = join_string(split_filename(:n_parts - 1), '.')
    else
      filename_without_extension = split_filename(1)%c
    end if
  end function filename_without_extension

  pure function delete_char_from_string(string, char) result(new_string)
    character(*), intent(in) :: string
    character, intent(in) :: char
    character(:), allocatable :: new_string
    integer :: i, i_new
    character(*), parameter :: proc_name = 'delete_char_from_string'

    allocate(character(len(string) - character_count(string, char)) :: new_string)
    i_new = 1
    do i = 1, len(string)
      if(string(i:i) /= char) then
        if(i_new > len(new_string)) call error_msg(proc_name, 'New string length exceeded')
        new_string(i_new:i_new) = string(i:i)
        i_new = i_new + 1
      end if
    end do
  end function delete_char_from_string

  pure function get_unique_chars(string) result(unique_chars)
    character(*), intent(in) :: string
    character(:), allocatable :: unique_chars
    character :: char
    integer :: i

    unique_chars = string
    i = 1
    do
      char = unique_chars(i:i)
      if(character_count(unique_chars, char) > 1) then
        unique_chars = unique_chars(:i) // delete_char_from_string(unique_chars(i+1:), char)
      end if
      if(i >= len(unique_chars)) exit
      i = i + 1
    end do
  end function get_unique_chars

  pure function get_chars_not_in_list(string, list) result(not_in_list)
    character(*), intent(in) :: string, list
    character(:), allocatable :: list_unique, not_in_list
    integer :: i, not_len

    list_unique = get_unique_chars(list)
    allocate(character(len(string)) :: not_in_list)
    not_len = 0
    do i = 1, len(string)
      if(.not. string_contains_character(list, string(i:i))) then
        not_len = not_len + 1
        not_in_list(not_len:not_len) = string(i:i)
      end if
    end do
    not_in_list = not_in_list(:not_len)
  end function get_chars_not_in_list

  pure function join_string(split_str, delimiter)
    type(string), intent(in) :: split_str(:)
    character, intent(in) :: delimiter
    character(:), allocatable :: join_string
    integer :: n_parts, n_chars_join, i, char_start_idx, delim_idx
    integer, dimension(:), allocatable :: n_chars_part

    n_parts = size(split_str)
    if(n_parts == 1) then
      allocate(character(len(split_str(1)%c)) :: join_string)
      join_string = split_str(1)%c
      return
    end if
    allocate(n_chars_part(n_parts))
    do concurrent(i = 1:n_parts)
      n_chars_part(i) = len(split_str(i)%c)
    end do
    n_chars_join = sum(n_chars_part) + n_parts - 1
    allocate(character(n_chars_join) :: join_string)
    char_start_idx = 1
    do i = 1, n_parts
      delim_idx = char_start_idx + n_chars_part(i)
      join_string(char_start_idx:delim_idx - 1) = split_str(i)%c
      if(i < n_parts) then
        join_string(delim_idx:delim_idx) = delimiter
        char_start_idx = delim_idx + 1
      end if
    end do
  end function join_string

  subroutine output_real_to_csv(a, filename)
    real(wp), intent(in) :: a(:)
    character(*), intent(in) :: filename
    character(*), parameter :: csv_format = '(*(es24.16e3,:,", "))'

    open(unit=10, access='sequential', action='write', status='replace', file=filename, form='formatted') 
    write(10, csv_format) a
    close(10)
  end subroutine output_real_to_csv

  real(wp) function get_walltime(t_wallcount_start, t_wallrate) result(t_wall)
    integer(int64), intent(in) :: t_wallcount_start, t_wallrate
    integer(int64) :: t_wallcount

    call system_clock(t_wallcount)
    t_wall = real((t_wallcount - t_wallcount_start), wp)/t_wallrate
  end function

  impure elemental logical function co_all_equal_real(a)
    real(wp), intent(in) :: a
    real(wp) :: a_co_min, a_co_max

    a_co_min = a
    a_co_max = a
    call co_max(a_co_max)
    call co_min(a_co_min)
    co_all_equal_real = (a_co_min == a_co_max)
  end function co_all_equal_real

  impure elemental logical function co_all_equal_int(a)
    integer, intent(in) :: a
    integer :: a_co_min, a_co_max

    a_co_min = a
    a_co_max = a
    call co_max(a_co_max)
    call co_min(a_co_min)
    co_all_equal_int = (a_co_min == a_co_max)
  end function co_all_equal_int

  impure elemental logical function co_all_equal_int64(a)
    integer(int64), intent(in) :: a
    integer(int64) :: a_co_min, a_co_max

    a_co_min = a
    a_co_max = a
    call co_max(a_co_max)
    call co_min(a_co_min)
    co_all_equal_int64 = (a_co_min == a_co_max)
  end function co_all_equal_int64

  subroutine print_message(string, step_opt, co_img_opt)
    character(len = *), intent(in) :: string
    integer(int64), intent(in), optional :: step_opt
    integer, intent(in), optional :: co_img_opt
    integer :: co_img

    if(present(co_img_opt)) then
      co_img = co_img_opt
    else
      co_img = this_image()
    end if
    if(this_image() == co_img) then
      call img_printhead(step_opt)
      write(output_unit, '(a)') string
    end if
  end subroutine print_message

  pure real(wp) function real_derf(a)
    real(wp), intent(in) :: a

    real_derf = 2/sqrt(pi)*exp(-a**2)
  end function real_derf

  ! Wrapper function needed to pass intrinsic erf as a dummy procedure
  pure real(wp) function real_erf(a)
    real(wp), intent(in) :: a

    real_erf = erf(a)
  end function real_erf

  pure elemental real(wp) function erfinv(a)
    real(wp), intent(in) :: a

    if(a <= -1.0_wp) then
      erfinv = ieee_value(1.0_wp, ieee_negative_inf)
    else if(a >= 1.0_wp) then
      erfinv = ieee_value(1.0_wp, ieee_positive_inf)
    else
      erfinv = findroot_newton(real_erf, real_derf, a, 100*epsilon(1.0_wp), -12.0_wp, 12.0_wp, y_target_opt = a)
    end if
  end function erfinv

  pure real(wp) function true_ratio(a)
    logical, intent(in) :: a(:)

    true_ratio = real(count(a), wp)/size(a)
  end function true_ratio

  pure subroutine error_msg(caller_name_opt, message_opt, step_opt, module_name_opt)
    character(*), intent(in), optional :: caller_name_opt, message_opt, module_name_opt
    integer(int64), intent(in), optional :: step_opt
    character(:), allocatable :: error_message

    error_message = trim(img_head(step_opt)) // ' ERROR: '
    if(present(module_name_opt)) error_message = error_message // module_name_opt // ': '
    if(present(caller_name_opt)) error_message = error_message // caller_name_opt // ': '
    if(present(message_opt)) then
      error_message = error_message // message_opt
    else
      error_message = error_message // 'Unspecified error!'
    end if
    error stop error_message
  end subroutine error_msg

  subroutine warning_msg(caller_name_opt, message_opt, step_opt)
    character(*), intent(in), optional :: caller_name_opt, message_opt
    integer(int64), intent(in), optional :: step_opt
    character(:), allocatable :: warning_message

    warning_message = trim(img_head(step_opt)) // ' WARNING: '
    if(present(caller_name_opt)) warning_message = warning_message // caller_name_opt // ': '
    if(present(message_opt)) then
      warning_message = warning_message // message_opt
    else
      warning_message = warning_message // 'Unspecified warning!'
    end if
    write(output_unit, '(a)') warning_message
  end subroutine warning_msg

  subroutine start_critical_in_order()
    if(this_image() > 1) sync images(left_image())
  end subroutine start_critical_in_order

  subroutine end_critical_in_order()
    if(this_image() < num_images()) sync images(right_image())
  end subroutine end_critical_in_order

  pure function get_lock_filename(filename) result(lock_filename)
    character(*), intent(in) :: filename
    character(:), allocatable :: lock_filename

    lock_filename = filename // '.lock'
  end function get_lock_filename

  subroutine wait_for_file_access(filename)
    character(*), intent(in) :: filename
    integer, parameter :: wait_ms = 10
    integer, parameter :: max_wait_ms = 1000
    integer, parameter :: max_it = max_wait_ms/wait_ms
    character(*), parameter :: proc_name = 'wait_for_file_access'
    character(:), allocatable :: lock_filename
    logical :: is_locked
    integer :: i

    lock_filename = get_lock_filename(filename)
    do i = 1, max_it
      inquire(file = lock_filename, exist = is_locked)
      if(.not. is_locked) return
      call millisleep(10)
    end do
    call error_msg(proc_name, 'Could not access file ' // filename)
  end subroutine wait_for_file_access

  subroutine lock_file_access(filename, lock_file_unit)
    character(*), intent(in) :: filename
    integer, intent(out) :: lock_file_unit
    character(:), allocatable :: lock_filename
    integer :: open_stat
    character(*), parameter :: proc_name = 'lock_file_access'

    lock_filename = get_lock_filename(filename)
    open(status = 'NEW', newunit = lock_file_unit, file = lock_filename, iostat = open_stat)
    if(open_stat /= 0) call error_msg(proc_name, 'Could not create lock file for file ' // filename)
  end subroutine lock_file_access

  subroutine unlock_file_access(filename, lock_file_unit)
    character(*), intent(in) :: filename
    integer, intent(in) :: lock_file_unit
    character(:), allocatable :: lock_filename
    integer :: close_stat
    character(*), parameter :: proc_name = 'unlock_file_access'

    lock_filename = get_lock_filename(filename)
    close(status = 'DELETE', unit = lock_file_unit, iostat = close_stat)
    if(close_stat /= 0) call error_msg(proc_name, 'Could not delete lock file for file ' // filename)
  end subroutine unlock_file_access

  subroutine exclusive_open(filename, file_unit, lock_file_unit, readonly_opt)
    character(*), intent(in) :: filename
    integer, intent(out) :: file_unit, lock_file_unit
    logical, intent(in), optional :: readonly_opt
    logical :: readonly
    character(*), parameter :: proc_name = 'exclusive_open'
    character(9) :: open_mode
    integer :: open_stat

    if(present(readonly_opt)) then
      readonly = readonly_opt
    else
      readonly = .false.
    end if

    if(readonly) then
      open_mode = 'READ'
    else
      open_mode = 'READWRITE'
    end if

    call wait_for_file_access(filename)
    call lock_file_access(filename, lock_file_unit)
    open(status = 'UNKNOWN', action = trim(open_mode), newunit = file_unit, file = filename, iostat = open_stat)
    if(open_stat /= 0) call error_msg(proc_name, 'Could not open file ' // filename)
  end subroutine exclusive_open

  subroutine exclusive_close(filename, file_unit, lock_file_unit)
    character(*), intent(in) :: filename
    integer, intent(in) :: file_unit, lock_file_unit
    character(*), parameter :: proc_name = 'exclusive_close'
    integer :: close_stat

    close(unit = file_unit, status = 'KEEP', iostat = close_stat)
    if(close_stat /= 0) call error_msg(proc_name, 'Could not close file ' // filename)
    call unlock_file_access(filename, lock_file_unit)
  end subroutine exclusive_close

  pure subroutine unique_ints(ints, ints_unique)
    integer, dimension(:), intent(in) :: ints
    integer, dimension(:), intent(out), allocatable :: ints_unique
    logical, dimension(size(ints)) :: select
    integer :: i, j

    select = .true.
    do j = 1, size(ints) - 1
      do i = j + 1, size(ints)
        if(ints(i) == ints(j)) then
          select(j) = .false.
          exit
        end if
      end do
    end do
    ints_unique = pack(ints, select)
  end subroutine unique_ints

  pure elemental real(wp) function option_or_default_real(option, default) result(arg)
    real(wp), intent(in), optional :: option
    real(wp), intent(in) :: default

    if(present(option)) then
      arg = option
    else
      arg = default
    end if
  end function option_or_default_real

  pure integer elemental function option_or_defalt_integer(option, default) result(arg)
    integer, intent(in), optional :: option
    integer, intent(in) :: default

    if(present(option)) then
      arg = option
    else
      arg = default
    end if
  end function option_or_defalt_integer

  pure logical elemental function option_or_defalt_logical(option, default) result(arg)
    logical, intent(in), optional :: option
    logical, intent(in) :: default

    if(present(option)) then
      arg = option
    else
      arg = default
    end if
  end function option_or_defalt_logical

  pure function option_or_defalt_character(option, default) result(arg)
    character(*), intent(in), optional :: option
    character(*), intent(in) :: default
    character(:), allocatable :: arg

    if(present(option)) then
      arg = option
    else
      arg = default
    end if
  end function option_or_defalt_character

  subroutine output_real_image(a, name, res_x, a_min_opt, a_max_opt, co_img_opt)
    real(wp), intent(in) :: a(:)
    character(*), intent(in) :: name
    integer, intent(in) :: res_x
    real(wp), optional, intent(in) :: a_min_opt, a_max_opt
    logical, optional, intent(in) :: co_img_opt
    integer :: res_y, y_line, pgm_unit
    character(1) :: a_byte(res_x)
    character(128) :: filename
    real(wp) :: a_min, a_max

    res_y = size(a)//res_x
    filename = get_pgm_filename(name, co_img_opt)
    a_max = option_or_default(a_max_opt, maxval(a))
    a_min = option_or_default(a_min_opt, minval(a))

    call create_pgm_file(filename, res_x, res_y, pgm_unit, &
      comment_opt = 'min = ' // convert_to_char(a_min) // ', max = ' // convert_to_char(a_max))

    do y_line = 1, res_y
      if(a_max - a_min < 128*spacing(max(abs(a_max), abs(a_min)))) then
        a_byte = char(0)
      else
        a_byte = char(max(min(nint(255*(a((y_line - 1)*res_x + 1:y_line*res_x) - a_min)/(a_max - a_min)), 255), 0))
      end if
      write(pgm_unit) a_byte
    end do
    close(pgm_unit)
  end subroutine output_real_image

  subroutine output_logical_image(a, name, res_x, co_img_opt)
    logical, intent(in) :: a(:)
    character(*), intent(in) :: name
    integer, intent(in) :: res_x
    logical, optional, intent(in) :: co_img_opt
    character(128) :: filename
    integer :: res_y, y_line, pgm_unit
    character(1) :: a_byte(res_x)

    res_y = size(a)//res_x
    filename = get_pgm_filename(name, co_img_opt)
    call create_pgm_file(filename, res_x, res_y, pgm_unit, 1)
    do y_line = 1, res_y
      where(a((y_line - 1)*res_x + 1:y_line*res_x))
        a_byte = char(1)
      elsewhere
        a_byte = char(0)
      end where
      write(pgm_unit) a_byte
    end do
    close(pgm_unit)
  end subroutine output_logical_image

  subroutine output_integer_image(a, name, res_x, a_min_opt, a_max_opt, co_img_opt)
    integer, intent(in) :: a(:)
    character(*), intent(in) :: name
    integer, intent(in) :: res_x
    integer, optional, intent(in) :: a_min_opt, a_max_opt
    logical, optional, intent(in) :: co_img_opt
    integer :: res_y, y_line, pgm_unit, a_max_out
    character(1) :: a_byte(res_x)
    character(128) :: filename
    logical :: rescale_needed
    integer :: a_min, a_max

    a_max = option_or_default(a_max_opt, maxval(a))
    a_min = option_or_default(a_min_opt, minval(a))
    res_y = size(a)//res_x
    filename = get_pgm_filename(name, co_img_opt)

    rescale_needed = a_min /= 0 .or. a_max < 0 .or. a_max > 255
    if(rescale_needed) then
      a_max_out = max(min(255, a_max - a_min), 1)
    else
      a_max_out = a_max
    end if

    call create_pgm_file(filename, res_x, res_y, pgm_unit, a_max_out)

    do y_line = 1, res_y
      if(a_max == a_min) then
        a_byte = char(0)
      else if(rescale_needed) then
        a_byte = char(max(min(a_max_out*(a((y_line - 1)*res_x + 1:y_line*res_x) - a_min)/(a_max - a_min), a_max_out), 0))
      else
        a_byte = char(max(min(a((y_line - 1)*res_x + 1:y_line*res_x), 255), 0))
      end if
      write(pgm_unit) a_byte
    end do
    close(pgm_unit)
  end subroutine output_integer_image

  pure character(128) function get_pgm_filename(name, co_img_opt) result(filename)
    character(*), intent(in) :: name
    logical, optional, intent(in) :: co_img_opt
    character(*), parameter :: pgm_extension = '.pgm'

    filename = name//'.pgm'
    if(present(co_img_opt)) then
      if(co_img_opt) write (filename, '(a,a4,i2.2,a4)') trim(name), '_img', this_image(), '.pgm'
    end if
  end function get_pgm_filename

  subroutine create_pgm_file(filename, res_x, res_y, pgm_unit, max_val_opt, comment_opt)
    character(*), intent(in) :: filename
    integer, intent(in) :: res_x, res_y
    integer, intent(out) :: pgm_unit
    integer, intent(in), optional :: max_val_opt
    character(*), intent(in), optional :: comment_opt
    integer :: max_val

    if(present(max_val_opt)) then
      max_val = max_val_opt
    else
      max_val = 255
    end if

    open(file=trim(filename), newunit=pgm_unit, form='unformatted', access='stream', status='replace')
    write(pgm_unit) 'P5' // NEW_LINE('A')
    if(present(comment_opt)) then
      write(pgm_unit) '#' // comment_opt // NEW_LINE('A')
    end if
    write(pgm_unit) convert_to_char(res_x) // ' ' // convert_to_char(res_y)
    write(pgm_unit) NEW_LINE('A') // convert_to_char(max_val) // NEW_LINE('A')

  end subroutine create_pgm_file

  pure function int_to_char(i) result(c)
    integer, intent(in) :: i
    character(:), allocatable :: c

    allocate(character(128) :: c)
    write(c, '(i0)') i
    c = trim(c)
  end function int_to_char

  pure function int64_to_char(i) result(c)
    integer(int64), intent(in) :: i
    character(:), allocatable :: c

    allocate(character(128) :: c)
    write(c, '(i0)') i
    c = trim(c)
  end function int64_to_char

  pure function real_to_char(a) result(c)
    real(wp), intent(in) :: a
    character(:), allocatable :: c

    allocate(character(128) :: c)
    write(c, '(g0.4)') a
    c = trim(c)
  end function real_to_char

  pure function realarray_to_char(a) result(c)
    real(wp), intent(in) :: a(:)
    character(:), allocatable :: c
    integer :: i

    c = real_to_char(a(1))
    do i = 2,size(a)
      c = c // ' ' // real_to_char(a(i))
    end do
  end function realarray_to_char

  pure function realmatrix_to_char(a) result(c)
    real(wp), intent(in) :: a(:,:)
    character(:), allocatable :: c
    integer :: i, j

    c = realarray_to_char(a(:,1))
    do j = 1,size(a,1)
      do i = 1,size(a,2)
        c = c // ' ' // real_to_char(a(j,i))
      end do
    end do
  end function realmatrix_to_char

  pure function logical_to_char(a) result(c)
    logical, intent(in) :: a
    character(1) :: c

    if(a) then
      c = 'T'
    else
      c = 'F'
    end if
  end function logical_to_char

  pure function intarray_to_char(i) result(c)
    integer, intent(in) :: i(:)
    character(:), allocatable :: c
    integer :: idx

    c = int_to_char(i(1))
    do idx = 2,size(i)
      c = c // ' ' // int_to_char(i(idx))
    end do
  end function intarray_to_char

  pure function concatenate_character_array(c, separator_opt) result(cat)
    character(*), dimension(:), intent(in) :: c
    character, intent(in), optional :: separator_opt
    character(:), allocatable :: cat
    character(12) :: format_string
    character :: separator
    integer :: cat_length

    separator = option_or_default(separator_opt, ' ')
    cat_length = (len(c(1)) + 1)*size(c)
    allocate(character(cat_length) :: cat)
    write(format_string, '(*(a))') '(*(a,:,"', separator, '"))'
    write(cat, format_string) c
  end function concatenate_character_array

  pure function cartesian_product_int(a) result(prod)
    type(int_vector), intent(in) :: a(:)
    type(int_vector), allocatable :: prod(:)
    integer, allocatable :: a_idx(:)
    integer :: prod_len, n_dims, dim, prod_idx

    prod_len = 1
    n_dims = size(a)
    allocate(a_idx(n_dims))
    a_idx = 1
    do dim = 1, n_dims
      prod_len = prod_len*size(a(dim)%v)
    end do
    allocate(prod(prod_len))
    do prod_idx = 1, prod_len
      allocate(prod(prod_idx)%v(n_dims))
      do dim = 1, n_dims
        prod(prod_idx)%v(dim) = a(dim)%v(a_idx(dim))
      end do
      do dim = 1, n_dims
        if(a_idx(dim) /= size(a(dim)%v)) then
          a_idx(dim) = a_idx(dim) + 1
          exit
        else
          a_idx(dim) = 1
        end if
      end do
    end do
  end function cartesian_product_int

  pure function cartesian_product_string(a) result(prod)
    type(string_vector), intent(in) :: a(:)
    type(string_vector), allocatable :: prod(:)
    integer, allocatable :: a_idx(:)
    integer :: prod_len, n_dims, dim, prod_idx

    prod_len = 1
    n_dims = size(a)
    allocate(a_idx(n_dims))
    a_idx = 1
    do dim = 1, n_dims
      prod_len = prod_len*size(a(dim)%s)
    end do
    allocate(prod(prod_len))
    do prod_idx = 1, prod_len
      allocate(prod(prod_idx)%s(n_dims))
      do dim = 1, n_dims
        prod(prod_idx)%s(dim)%c = a(dim)%s(a_idx(dim))%c
      end do
      do dim = 1, n_dims
        if(a_idx(dim) /= size(a(dim)%s)) then
          a_idx(dim) = a_idx(dim) + 1
          exit
        else
          a_idx(dim) = 1
        end if
      end do
    end do
  end function cartesian_product_string

  pure elemental type(string) function concatenate_string_vector(a, separator, discard_empty_opt) result(cat)
    type(string_vector), intent(in) :: a
    character(*), intent(in) :: separator
    logical, intent(in), optional :: discard_empty_opt
    logical :: discard_empty
    integer :: cat_len, a_idx
    integer, allocatable, dimension(:) :: str_len
    type(array_index) :: char_idx

    discard_empty = option_or_default(discard_empty_opt, .false.)
    allocate(str_len(size(a%s)))
    do a_idx = 1, size(a%s)
      str_len(a_idx) = len(a%s(a_idx)%c)
    end do
    if(.not. discard_empty) then
      cat_len = sum(str_len) + (size(a%s) - 1)*len(separator)
    else
      cat_len = max(sum(str_len) + (count(str_len > 0) - 1)*len(separator), 0)
    end if
    allocate(character(cat_len) :: cat%c)
    char_idx = array_index(1, 1)
    do a_idx = 1, size(a%s)
      if(str_len(a_idx) == 0 .and. discard_empty) cycle
      if(char_idx%first > 1 .or. (a_idx > 1 .and. .not. discard_empty) .and. len(separator) > 0) then
        char_idx%last = char_idx%first + len(separator) - 1
        cat%c(char_idx%first:char_idx%last) = separator
        char_idx%first = char_idx%last + 1
      end if
      if(str_len(a_idx) == 0) cycle
      char_idx%last = char_idx%first + str_len(a_idx) - 1
      cat%c(char_idx%first:char_idx%last) = a%s(a_idx)%c
      char_idx%first = char_idx%last + 1
    end do
  end function concatenate_string_vector

  pure logical function wildcard_string_match(str, substr) result(match)
    character(*), intent(in) :: str, substr
    integer :: match_idx, prev_match_idx
    type(string), allocatable :: substr_split(:)
    integer :: i
    character, parameter :: wildcard = '*'
    character(*), parameter :: proc_name = 'wildcard_string_match'

    match = .true.
    call split_string(substr, wildcard, substr_split)
    prev_match_idx = 1
    if(index(str, wildcard) > 0) call error_msg(proc_name, 'Given string contains wildcard characters')
    do i = 1, size(substr_split)
      if(len(substr_split(i)%c) == 0) cycle
      match_idx = index(str(prev_match_idx:), substr_split(i)%c)
      if(match_idx < 1 &
        .or. (i == 1 .and. match_idx /= 1) &
        .or. (i == size(substr_split) .and. match_idx + len(substr_split(i)%c) - 1 /= len(str))) then
        match = .false.
        return
      end if
      prev_match_idx = prev_match_idx + match_idx + len(substr_split(i)%c) - 1
    end do
  end function wildcard_string_match

  pure elemental logical function sign_is_equal_real(a, b) result(sign_is_equal)
    real(wp), intent(in) :: a, b

    sign_is_equal = ((a < 0 .and. b < 0) .or. (a >= 0 .and. b >= 0))
  end function sign_is_equal_real

  pure elemental real(wp) function nearest_divisor_less_or_equal(a, target_val)
    real(wp), intent(in) :: a, target_val
    character(*), parameter :: proc_name = 'nearest_divisor_less_or_equal'

    if(a == 0) call error_msg(proc_name, '"a" must not be zero')
    if(target_val == 0) call error_msg(proc_name, '"target_val" must not be zero')
    if(.not. sign_is_equal(a, target_val)) call error_msg(proc_name, 'Signs of "a" and "target_val" must be equal')
    nearest_divisor_less_or_equal = target_val/(ceiling(target_val/a))
  end function nearest_divisor_less_or_equal

  pure subroutine img_partition(nval, nval_per_img, val_lo, val_hi, img_opt, nimg_opt)
    integer, intent(in) :: nval
    integer, optional, intent(in) :: img_opt, nimg_opt
    real(wp), intent(out) :: nval_per_img
    integer, intent(out) :: val_lo, val_hi
    integer :: img, nimg

    nimg = option_or_default(nimg_opt, num_images())
    img = option_or_default(img_opt, this_image())
    nval_per_img = real(nval, wp)/nimg
    val_lo = nint((img - 1)*nval_per_img) + 1
    val_hi = nint(img*nval_per_img)
  end subroutine img_partition

  subroutine resize_real_array(a, target_size)
    real(wp), intent(inout), allocatable :: a(:)
    integer, intent(in) :: target_size
    real(wp), allocatable :: a_temp(:)

    if(size(a) /= target_size) then
      allocate(a_temp(target_size))
      a_temp = a(1:target_size)
      deallocate(a)
      call move_alloc(a_temp, a)
    end if
  end subroutine resize_real_array

  subroutine resize_int_array(a, target_size)
    integer, intent(inout), allocatable :: a(:)
    integer, intent(in) :: target_size
    integer, allocatable :: a_temp(:)

    if(size(a) /= target_size) then
      allocate(a_temp(target_size))
      a_temp = a(1:target_size)
      deallocate(a)
      call move_alloc(a_temp, a)
    end if
  end subroutine resize_int_array

  subroutine concatenate_image_int_arrays(a)
    integer, intent(inout), allocatable :: a(:)
    integer, allocatable :: a_local(:)
    integer :: total_size, start_idx, end_idx
    integer, allocatable :: local_size(:)

    allocate(local_size(num_images()))
    local_size = 0
    local_size(this_image()) = size(a)
    call co_sum(local_size)
    total_size = sum(local_size)
    call move_alloc(a, a_local)
    allocate(a(total_size))
    a = 0
    if(local_size(this_image()) > 0) then
      start_idx = sum(local_size(1:this_image() - 1)) + 1
      end_idx = start_idx + local_size(this_image()) - 1
      a(start_idx:end_idx) = a_local
    end if
    call co_sum(a)
  end subroutine concatenate_image_int_arrays

  subroutine get_allocatable_command_argument(arg_no, arg)
    integer, intent(in) :: arg_no
    character(:), allocatable, intent(inout) :: arg
    integer :: arg_len

    if(allocated(arg)) deallocate(arg)
    call get_command_argument(arg_no, length = arg_len)
    allocate(character(arg_len) :: arg)
    call get_command_argument(arg_no, arg)
  end subroutine get_allocatable_command_argument

end module mod_base_functions
