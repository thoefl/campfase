! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_spectral

  !TODO: could in-place-transforms increase FFT performance?
  !      Theoretically possible, but benchmark results show no clear benefit for larger FFT sizes
  !TODO: Could other FFT libraries provide benefit, e.g. FFTE (written in Fortran)?

  use mod_globals, only: verbosity_level
  use mod_base_functions
  use mod_axis
  use mod_regular_grid
  use mod_time_dep_vars
  use mod_fcdata
  use mod_parameters
  use mod_fftw

  implicit none
  private

  public :: fft_load_wisdom

  type, public :: fft_field
    !dir$ attributes align : 64 :: k2_fft, k4_fft, fft_vals_backup
    type(time_dependent_real) :: real, real_laplacian
    type(time_dependent_complex) :: fft, fft_prev, fft_laplacian
    real(wp), allocatable, dimension(:) :: k2_fft, k4_fft
    type(real_vector), allocatable, dimension(:) :: k_fft, one_dim_real
    type(time_dependent_real), allocatable, dimension(:) :: real_grad
    type(complex_vector), allocatable, dimension(:) :: one_dim_fft
    complex(wp), allocatable, dimension(:) :: fft_vals_backup
    type(regular_grid), private :: fft_grid
    logical :: keep_prev_fft, backup_backward_transform
    type(c_ptr) :: fprefs_all_dim, bprefs_all_dim
    type(c_ptr), allocatable, dimension(:) :: fprefs_one_dim, bprefs_one_dim
    character(:), allocatable, private :: var_name
    integer(int64) :: n_forward_alldim = 0
    integer(int64) :: n_backward_alldim = 0
  contains
    procedure, public :: init => fft_field_init
    procedure, public :: reset => fft_field_reset
    procedure, public :: forward_alldim => fft_forward_alldim
    procedure, public :: backward_alldim => fft_backward_alldim
    procedure, public :: calc_grad => fft_calc_grad
    procedure, public :: calc_real_laplacian => fft_calc_real_laplacian
    procedure, public :: calc_fft_laplacian => fft_calc_fft_laplacian
    procedure, public :: nabla_square => fft_nabla_square
    procedure, public :: convolve => fft_convolve
    procedure, public :: swap_fft => fft_swap_fft
    procedure, public :: fft_allocate_check => fft_fft_allocate_check
  end type fft_field

  public :: fft_grid_from_real, fft_field_broadcast, fft_nabla_product, spectral_solve

  interface fft_field_broadcast
    module procedure :: fft_field_array_broadcast, fft_field_single_broadcast
  end interface fft_field_broadcast

contains

  pure elemental type(regular_grid) function fft_grid_from_real(real_grid) result(fft_grid)
    type(regular_grid), intent(in) :: real_grid
    type(axis), allocatable :: fft_axes(:)
    integer :: dim, n_nodes_fft
    character(*), parameter :: fft_axes_suffix = '_fft'

    allocate(fft_axes(real_grid%n_dims))
    do dim = 1, real_grid%n_dims
      associate(fft_axis => fft_axes(dim), &
          real_axis => real_grid%axes(dim))
        if(dim == 1) then
          n_nodes_fft = (real_axis%n_nodes//2) + 1
        else
          n_nodes_fft = real_axis%n_nodes
        end if
        fft_axis = axis(n_nodes_fft, 2*pi/real_axis%length, real_axis%bound_cond, real_axis%label // fft_axes_suffix)
      end associate
    end do

    call fft_grid%init(fft_axes)
  end function fft_grid_from_real

  subroutine fft_field_init(this, real_grid, var_name, keep_prev_fft_opt)
  class(fft_field), intent(inout) :: this
    type(regular_grid), intent(in) :: real_grid
    logical, intent(in), optional :: keep_prev_fft_opt
    type(real_vector), dimension(real_grid%n_dims) :: k2_fft_dim
    integer :: i, j, dim, node
    integer, allocatable, dimension(:) :: node_coords
    integer, parameter :: two_dimensional(2) = [1, 2], x_dimension(1) = [1], y_dimension(1) = [2]
    character(*), intent(in) :: var_name

    this%backup_backward_transform = .not. backwards_fft_preserves_input()
    call this%real%init(real_grid, trim(var_name))
    this%var_name = trim(var_name)
    this%keep_prev_fft = option_or_default(keep_prev_fft_opt, .false.)

    allocate(this%one_dim_real(real_grid%n_dims))
    allocate(this%one_dim_fft(real_grid%n_dims))
    allocate(this%k_fft(real_grid%n_dims))
    allocate(this%real_grad(real_grid%n_dims))
    allocate(this%fprefs_one_dim(real_grid%n_dims))
    allocate(this%bprefs_one_dim(real_grid%n_dims))

    this%fft_grid = fft_grid_from_real(real_grid)

    do dim = 1, real_grid%n_dims
      allocate(this%k_fft(dim)%vals((real_grid%axes(dim)%n_nodes//2) + 1))
      this%k_fft(dim)%vals = [(i*this%fft_grid%axes(dim)%spacing, i=0, (real_grid%axes(dim)%n_nodes//2) - 1), 0.0_wp]

      allocate(k2_fft_dim(dim)%vals(this%fft_grid%axes(dim)%n_nodes))
      if(dim == 1) then
        k2_fft_dim(dim)%vals = [(i*this%fft_grid%axes(dim)%spacing, i=0, this%fft_grid%axes(dim)%n_nodes - 1)]
      else
        k2_fft_dim(dim)%vals = cshift([(i*this%fft_grid%axes(dim)%spacing, i=-max(this%fft_grid%axes(dim)%n_nodes/2, 1) + 1,&
          this%fft_grid%axes(dim)%n_nodes/2)], this%fft_grid%axes(dim)%n_nodes/2 - 1)
      end if
    end do
    allocate(this%k2_fft(this%fft_grid%n_nodes_total))
    allocate(this%k4_fft(this%fft_grid%n_nodes_total))
    this%k2_fft = 0.0_wp

    do concurrent(node = 1:this%fft_grid%n_nodes_total)
      node_coords = this%fft_grid%linear_to_dim_idx(node)
      do concurrent(dim = 1:this%fft_grid%n_dims)
        this%k2_fft(node) = this%k2_fft(node) + k2_fft_dim(dim)%vals(node_coords(dim))**2
      end do
    end do
    do concurrent(i = 1:this%fft_grid%axes(1)%n_nodes, j = 1:real_grid%axes(2)%n_nodes)
      this%k2_fft(i + (j - 1)*this%fft_grid%axes(1)%n_nodes) = k2_fft_dim(1)%vals(i)**2 + k2_fft_dim(2)%vals(j)**2
    end do
    this%k4_fft = this%k2_fft**2

  end subroutine fft_field_init

  pure elemental subroutine fft_field_reset(this)
  class(fft_field), intent(inout) :: this
    integer :: dim

    call this%real%invalidate()
    call this%fft%invalidate()
    call this%fft_prev%invalidate()
    call this%fft_laplacian%invalidate()
    call this%real_laplacian%invalidate()
    if(allocated(this%real_grad)) then
      do dim = 1, this%real%grid%n_dims
        call this%real_grad(dim)%invalidate()
      end do
    end if
  end subroutine fft_field_reset

  pure subroutine fft_swap_fft(this)
  class(fft_field), intent(inout) :: this
    complex(wp), allocatable, dimension(:) :: temp
    real(wp) :: fft_time, fft_prev_time

    ! Store FFT times
    fft_time = this%fft%time
    fft_prev_time = this%fft_prev%time
    ! Invalidate all variables other than the current and previous FFT
    call this%reset()

    ! Perform the swap
    call move_alloc(this%fft%vals, temp)
    call move_alloc(this%fft_prev%vals, this%fft%vals)
    call move_alloc(temp, this%fft_prev%vals)

    ! Restore FFT times
    this%fft%time = fft_prev_time
    this%fft_prev%time = fft_time
  end subroutine fft_swap_fft

  subroutine fft_forward_alldim(this)
  class(fft_field), intent(inout) :: this
    character(*), parameter :: proc_name = 'fft_forward_alldim'
    !dir$ vector aligned

    call this%fft_allocate_check()
    call this%real%check_is_valid(proc_name)
    call this%real%check_time_is_greater_equal(this%fft, proc_name)
    if(this%fft%time_is_equal(this%real)) return
    call fftw_execute_dft_r2c(this%fprefs_all_dim, this%real%vals, this%fft%vals)
    this%fft%time = this%real%time
    this%n_forward_alldim = this%n_forward_alldim + 1
    !call fft_errcheck(fft_err)
  end subroutine fft_forward_alldim

  subroutine fft_fft_allocate_check(this)
  class(fft_field), intent(inout) :: this
    real(wp), dimension(:), allocatable :: vals_backup
    !dir$ vector aligned

    if(.not. allocated(this%fft%vals)) then
      call this%fft%init(this%fft_grid, this%var_name // '_fft')
      if(this%keep_prev_fft) call this%fft_prev%init(this%fft_grid, this%var_name // '_fft_prev')
      if(this%backup_backward_transform) allocate(this%fft_vals_backup(this%fft_grid%n_nodes_total))
      vals_backup = this%real%vals
      call fft_prefsetup_oop(this%real%vals, this%fft%vals, this%real%grid%axes,&
        this%fprefs_all_dim, this%bprefs_all_dim, .true., this%fft%var_name)
      this%real%vals = vals_backup
    end if  
  end subroutine fft_fft_allocate_check

  subroutine fft_backward_alldim(this)
  class(fft_field), intent(inout) :: this
    character(*), parameter :: proc_name = 'fft_backward_alldim'
    !dir$ vector aligned

    call this%fft%check_is_valid(proc_name)
    call this%fft%check_time_is_greater_equal(this%real, proc_name)
    if(this%real%time_is_equal(this%fft)) return
    if(this%backup_backward_transform) this%fft_vals_backup = this%fft%vals
    call fftw_execute_dft_c2r(this%bprefs_all_dim, this%fft%vals, this%real%vals)
    if(this%backup_backward_transform) this%fft%vals = this%fft_vals_backup
    this%real%vals = this%real%vals/this%real%grid%n_nodes_total
    this%real%time = this%fft%time
    this%n_backward_alldim = this%n_backward_alldim + 1
    !call fft_errcheck(fft_err)
  end subroutine fft_backward_alldim

  subroutine fft_calc_grad(this, grad_dim_opt)
  class(fft_field), intent(inout) :: this
    integer, intent(in), optional :: grad_dim_opt
    integer :: i, grad_dim, dim, node_start, node_end, node_stride
    integer, dimension(:), allocatable :: grad_dim_list
    integer, dimension(this%real%grid%n_dims) :: node_start_coords, node_end_coords
    character(*), parameter :: proc_name = 'fft_calc_grad'
    !dir$ vector aligned

    call this%real%check_is_valid(proc_name)
    call this%real%check_time_is_greater_equal(this%fft, proc_name)

    if(present(grad_dim_opt)) then
      grad_dim_list = [grad_dim_opt]
    else
      grad_dim_list = [(i, i=1, this%real%grid%n_dims)]
    end if

    !FIXME: use associates to clean this up
    do i = 1, size(grad_dim_list)
      grad_dim = grad_dim_list(i)
      call this%real%check_time_is_greater_equal(this%real_grad(grad_dim), proc_name)
      if(this%real_grad(grad_dim)%time_is_equal(this%real)) return
      if(.not. allocated(this%one_dim_fft(grad_dim)%vals)) then
        allocate(this%one_dim_real(grad_dim)%vals(this%real%grid%axes(grad_dim)%n_nodes))
        allocate(this%one_dim_fft(grad_dim)%vals(this%real%grid%axes(grad_dim)%n_nodes/2 + 1))
        call fft_prefsetup_oop(this%one_dim_real(grad_dim)%vals, this%one_dim_fft(grad_dim)%vals,&
          [this%real%grid%axes(grad_dim)], this%fprefs_one_dim(grad_dim), this%bprefs_one_dim(grad_dim), .false.,&
          this%var_name // '_fft_1d_dim_' // convert_to_char(grad_dim))
      end if
      if(.not. allocated(this%real_grad(grad_dim)%vals)) then
        call this%real_grad(grad_dim)%init(this%real%grid, this%var_name // '_real_grad_' // convert_to_char(grad_dim))
      end if

      if(grad_dim <= 0 .or. grad_dim > this%real%grid%n_dims) then
        error stop 'ERROR: ' // proc_name // ': grad_dim has to be inside the closed interval [1, n_dims]!'
      else if(grad_dim == 1) then
        node_stride = 1
      else
        node_stride = product(this%real%grid%axes(:grad_dim - 1)%n_nodes)
      end if

      node_start_coords = 1
      outer: do
        node_end_coords = node_start_coords
        node_end_coords(grad_dim) = this%real%grid%axes(grad_dim)%n_nodes
        node_start = this%real%grid%dim_to_linear_idx(node_start_coords)
        node_end = this%real%grid%dim_to_linear_idx(node_end_coords)
        this%one_dim_real(grad_dim)%vals = this%real%vals(node_start:node_end:node_stride)
        call fftw_execute_dft_r2c(this%fprefs_one_dim(grad_dim), this%one_dim_real(grad_dim)%vals,&
          this%one_dim_fft(grad_dim)%vals)
        !call fft_errcheck(fft_err)
        this%one_dim_fft(grad_dim)%vals = (0.0_wp,1.0_wp)*this%k_fft(grad_dim)%vals*this%one_dim_fft(grad_dim)%vals
        call fftw_execute_dft_c2r(this%bprefs_one_dim(grad_dim), this%one_dim_fft(grad_dim)%vals,&
          this%one_dim_real(grad_dim)%vals)
        !call fft_errcheck(fft_err)
        this%real_grad(grad_dim)%vals(node_start:node_end:node_stride) =&
          this%one_dim_real(grad_dim)%vals/this%real%grid%axes(grad_dim)%n_nodes

        do dim = 1, this%real%grid%n_dims
          if(dim == grad_dim) cycle
          if(node_start_coords(dim) /= this%real%grid%axes(dim)%n_nodes) then
            node_start_coords(dim) = node_start_coords(dim) + 1
            cycle outer
          else
            node_start_coords(dim) = 1
          end if
        end do
        exit
      end do outer
      this%real_grad(grad_dim)%time = this%real%time
    end do
  end subroutine fft_calc_grad

  subroutine fft_calc_fft_laplacian(this)
  class(fft_field), intent(inout) :: this
    character(*), parameter :: proc_name = 'fft_calc_fft_laplacian'
    integer :: i
    !dir$ vector aligned

    call this%forward_alldim()
    call this%fft%check_is_valid(proc_name)
    if(.not. allocated(this%fft_laplacian%vals)) then
      call this%fft_laplacian%init(this%fft%grid, this%var_name // '_fft_laplacian')
    end if
    if(this%fft_laplacian%time_is_equal(this%real)) return
    do concurrent(i = 1:this%fft%grid%n_nodes_total)
      this%fft_laplacian%vals(i) = -this%k2_fft(i)*this%fft%vals(i)
    end do
    this%fft_laplacian%time = this%fft%time
  end subroutine fft_calc_fft_laplacian

  subroutine fft_calc_real_laplacian(this)
  class(fft_field), intent(inout) :: this
    character(*), parameter :: proc_name = 'fft_calc_real_laplacian'
    !dir$ vector aligned

    call this%calc_fft_laplacian()
    call this%fft_laplacian%check_is_valid(proc_name)
    if(.not. allocated(this%real_laplacian%vals)) then
      call this%real_laplacian%init(this%real%grid, this%var_name // '_real_laplacian')
    end if
    if(this%real_laplacian%time_is_equal(this%real)) return
    call fftw_execute_dft_c2r(this%bprefs_all_dim, this%fft_laplacian%vals, this%real_laplacian%vals)
    this%real_laplacian%vals = this%real_laplacian%vals/this%real%grid%n_nodes_total
    this%real_laplacian%time = this%fft_laplacian%time
    !call fft_errcheck(fft_err)
  end subroutine fft_calc_real_laplacian

  pure elemental real(wp) function fft_nabla_square(this, node)
  class(fft_field), intent(in) :: this
    integer, intent(in) :: node
    integer :: dim
    !dir$ vector aligned

    fft_nabla_square = 0.0_wp
    do concurrent(dim = 1:this%real%grid%n_dims)
      fft_nabla_square = fft_nabla_square + this%real_grad(dim)%vals(node)**2
    end do
  end function fft_nabla_square

  pure elemental real(wp) function fft_nabla_product(a, b, i)
  class(fft_field), intent(in) :: a, b
    integer, intent(in) :: i
    integer :: dim
    !dir$ vector aligned

    fft_nabla_product = 0.0_wp
    do concurrent(dim = 1:a%real%grid%n_dims)
      fft_nabla_product = fft_nabla_product + a%real_grad(dim)%vals(i)*b%real_grad(dim)%vals(i)
    end do
  end function fft_nabla_product

  subroutine fft_convolve(this, a)
  class(fft_field), intent(inout) :: this
  class(fft_field), intent(in) :: a
    !dir$ vector aligned

    this%fft%vals = this%fft%vals*a%fft%vals
  end subroutine fft_convolve

  subroutine fft_field_array_broadcast(field, field_name, step, img_offset_opt)
  class(fft_field), intent(inout) :: field(:)
    character(*), intent(in), optional :: field_name
    integer(int64), intent(in), optional :: step
    integer, intent(in), optional :: img_offset_opt
    integer :: i, img_offset, broadcast_img
    character(*), parameter :: proc_name = 'fft_field_array_broadcast'

    if(present(img_offset_opt)) then
      img_offset = img_offset_opt
    else
      img_offset = 0
    end if
    if(present(field_name) .and. verbosity_level >= 5) call print_message('Waiting for ' // field_name // ' sync...', step)
    do i = 1, size(field)
      broadcast_img = img_offset + i
      if(this_image() == broadcast_img) then
        call field(i)%real%check_time_is_greater_equal(field(i)%fft, proc_name)
      end if
      call co_broadcast(field(i)%real%vals, broadcast_img)
      call co_broadcast(field(i)%real%time, broadcast_img)
    end do
    if(present(field_name) .and. verbosity_level >= 5) call print_message(field_name // ' synced!', step, 1)
  end subroutine fft_field_array_broadcast

  subroutine fft_field_single_broadcast(field, field_name, step, broadcast_img_opt)
  class(fft_field), intent(inout) :: field
    character(*), intent(in), optional :: field_name
    integer(int64), intent(in), optional :: step
    integer, intent(in), optional :: broadcast_img_opt
    integer :: broadcast_img
    character(*), parameter :: proc_name = 'fft_field_single_broadcast'

    if(present(broadcast_img_opt)) then
      broadcast_img = broadcast_img_opt
    else
      broadcast_img = 1
    end if
    if(present(field_name) .and. verbosity_level >= 5) call print_message('Waiting for ' // field_name // ' sync...', step)
    if(this_image() == broadcast_img) then
      call field%real%check_time_is_greater_equal(field%fft, proc_name)
    end if
    call co_broadcast(field%real%vals, broadcast_img)
    call co_broadcast(field%real%time, broadcast_img)
    if(present(field_name) .and. verbosity_level >= 5) call print_message(field_name // ' synced!', step, 1)
  end subroutine fft_field_single_broadcast    

  ! Solves equation of the form del_a/del_t = D*(a_coeff*Delta(a) + mu) with ffts as input
  pure subroutine spectral_solve(a, mu, D, dt, use_bdfab, bdfab_stage, a_coeff_opt)
  class(fft_field), intent(inout) :: a, mu
    real(wp), intent(in) :: D, dt
    real(wp), intent(in), optional :: a_coeff_opt
    logical, intent(in) :: use_bdfab
    integer, intent(in) :: bdfab_stage
    real(wp) :: a_coeff
    character(*), parameter :: proc_name = 'spectral_solve'
    !dir$ vector aligned

    call a%fft%check_time_is_greater_equal(a%real, proc_name)
    call mu%fft%check_time_is_greater_equal(mu%real, proc_name)
    call a%fft%check_time_is_equal(mu%fft, proc_name)
    if(use_bdfab .and. bdfab_stage > 0) then
      call a%fft%check_time_diff_within_tolerance(a%fft_prev, dt, proc_name)
      call mu%fft%check_time_diff_within_tolerance(mu%fft_prev, dt, proc_name)
    end if

    a_coeff = option_or_default(a_coeff_opt, 1.0_wp)

    if(use_bdfab) then
      if(bdfab_stage == 0) then
        a%fft_prev%vals = a%fft%vals
        mu%fft_prev%vals = mu%fft%vals
        a%fft_prev%time = a%fft%time
        mu%fft_prev%time = mu%fft%time
        a%fft%vals = (dt*D*(mu%fft%vals + a%k2_fft*stab_coeff_bdfab*a%fft%vals) + a%fft%vals)&
          /(1 + dt*a%k2_fft*D*(a_coeff + stab_coeff_bdfab))
      else if(bdfab_stage == 1) then
        a%fft%vals = (4*a%fft%vals - a%fft_prev%vals + 2*dt*D*(2*mu%fft%vals - mu%fft_prev%vals&
          + a%k2_fft*stab_coeff_bdfab*(2*a%fft%vals - a%fft_prev%vals)))/(3 + 2*dt*a%k2_fft*D*(a_coeff + stab_coeff_bdfab))
      else
        a%fft_prev%vals = (4*a%fft%vals - a%fft_prev%vals + 2*dt*D*(2*mu%fft%vals - mu%fft_prev%vals&
          + a%k2_fft*stab_coeff_bdfab*(2*a%fft%vals - a%fft_prev%vals)))/(3 + 2*dt*a%k2_fft*D*(a_coeff + stab_coeff_bdfab))
        a%fft_prev%time = a%fft%time
        call a%swap_fft()
      end if
    else
      a%fft%vals =&
        (dt*D*(mu%fft%vals + a%k2_fft*stab_coeff_be*a%fft%vals) + a%fft%vals)/(1 + dt*a%k2_fft*D*(a_coeff + stab_coeff_be))
    end if
    a%fft%time = a%fft%time + dt
  end subroutine

  subroutine fft_load_wisdom(real_grid)
    type(regular_grid), intent(in) :: real_grid
    integer :: n_fft_threads
    character(:), allocatable :: wisdom_filename_all_dim
    type(string), allocatable :: wisdom_filename_one_dim(:)
    integer, dimension(:), allocatable :: unique_dims
    character(*), parameter :: proc_name = 'fft_load_wisdom'
    logical :: wisdom_exists
    integer :: fftw_status, i

    n_fft_threads = fftw_planner_nthreads()
    call unique_ints(real_grid%axes%n_nodes, unique_dims)
    allocate(wisdom_filename_one_dim(size(unique_dims)))
    wisdom_filename_all_dim = 'fftw_wisdom_'
    do i = 1, real_grid%n_dims
      wisdom_filename_all_dim = wisdom_filename_all_dim // convert_to_char(real_grid%axes(i)%n_nodes)
      if(i < real_grid%n_dims) wisdom_filename_all_dim = wisdom_filename_all_dim // 'x'
    end do
    wisdom_filename_all_dim = wisdom_filename_all_dim // '_' // convert_to_char(n_fft_threads) // '_threads.txt'
    do i = 1, size(unique_dims)
      wisdom_filename_one_dim(i)%c = 'fftw_wisdom_' // convert_to_char(unique_dims(i)) // '_' &
        // convert_to_char(n_fft_threads) // '_threads.txt'
    end do
    critical
    inquire(file = wisdom_filename_all_dim, exist = wisdom_exists)
    if(wisdom_exists) then
      fftw_status = fftw_import_wisdom_from_filename(wisdom_filename_all_dim // c_null_char)
      if(fftw_status == 0) call warning_msg(proc_name, 'Could not load FFTW wisdom from file ' // wisdom_filename_all_dim)
    end if
    do i = 1, size(unique_dims)
      inquire(file = wisdom_filename_one_dim(i)%c, exist = wisdom_exists)
      if(wisdom_exists) then
        fftw_status = fftw_import_wisdom_from_filename(wisdom_filename_one_dim(i)%c // c_null_char)
        if(fftw_status == 0) then
          call warning_msg(proc_name, 'Could not load FFTW wisdom from file ' // wisdom_filename_one_dim(i)%c)
        end if
      end if
    end do
    end critical
  end subroutine fft_load_wisdom

  subroutine fft_prefsetup_oop(real_data, fft_data, real_axes, fft_fprefs, fft_bprefs, preserve_input, var_name)
    real(wp), intent(inout), dimension(:) :: real_data
    complex(wp), intent(inout), dimension(:) :: fft_data
    type(axis), intent(in) :: real_axes(:)
    logical, intent(in) :: preserve_input
    type(c_ptr), intent(out) :: fft_fprefs, fft_bprefs
    character(*), intent(in) :: var_name
    type(axis), dimension(size(real_axes)) :: fft_axes
    integer, dimension(size(real_axes)) :: real_size_inv, fft_size_inv
    integer :: i, fft_fflags, fft_bflags, fftw_status, n_dims, n_fft_threads
    character(:), allocatable :: wisdom_filename
    logical :: wisdom_exists
    character(*), parameter :: proc_name = 'fft_prefsetup_oop'
    integer, parameter :: max_read_attempts = 100

    !FFTW plan creation is not thread safe!
    n_fft_threads = fftw_planner_nthreads()
    wisdom_filename = 'fftw_wisdom_'
    n_dims = size(real_axes)
    do i = 1, n_dims - 1
      wisdom_filename = wisdom_filename // convert_to_char(real_axes(i)%n_nodes) // 'x'
    end do
    wisdom_filename = wisdom_filename // convert_to_char(real_axes(n_dims)%n_nodes)
    wisdom_filename = wisdom_filename // '_' // convert_to_char(n_fft_threads) // '_threads.txt'

    fft_axes = real_axes
    !FIXME: should check for odd node count
    fft_axes(1)%n_nodes = fft_axes(1)%n_nodes/2 + 1
    real_size_inv = real_axes(n_dims:1:-1)%n_nodes
    fft_size_inv = fft_axes(n_dims:1:-1)%n_nodes

    if(preserve_input) then
      fft_fflags = ior(FFTW_EXHAUSTIVE, FFTW_PRESERVE_INPUT)
    else
      fft_fflags = ior(FFTW_EXHAUSTIVE, FFTW_DESTROY_INPUT)
    end if
    fft_bflags = FFTW_EXHAUSTIVE

    if(verbosity_level >= 1) call print_message('Creating forward FFT plan for ' // var_name // '...')
    fft_fprefs = fftw_plan_many_dft_r2c(&
      rank = n_dims,&
      n = real_size_inv,&
      howmany = 1,&
      in = real_data,&
      inembed = real_size_inv,&
      istride = 1,&
      idist = 0,&
      out = fft_data,&
      onembed = fft_size_inv,&
      ostride = 1,&
      odist = 0,&
      flags = fft_fflags)

    if(verbosity_level >= 1) call print_message('Creating backward FFT plan for ' // var_name // '...')
    fft_bprefs = fftw_plan_many_dft_c2r(&
      rank = n_dims,&
      n = real_size_inv,&
      howmany = 1,&
      in = fft_data,&
      inembed = fft_size_inv,&
      istride = 1,&
      idist = 0,&
      out = real_data,&
      onembed = real_size_inv,&
      ostride = 1,&
      odist = 0,&
      flags = fft_bflags)

    if(verbosity_level >= 1) call print_message('FFT plan creation finished for ' // var_name // '!')
    inquire(file = wisdom_filename, exist = wisdom_exists)
    if(.not. wisdom_exists) then
      fftw_status = fftw_export_wisdom_to_filename(wisdom_filename // c_null_char)
      if(fftw_status == 0) call error_msg(proc_name, 'Could not save FFTW wisdom to file ' // wisdom_filename)
    end if
  end subroutine fft_prefsetup_oop

  ! FFTW does not offer any out of place complex-to-real transforms that preserve input (i.e. the complex part or FFT)
  ! MKL however has this option, so we need to check whether the complex part will be destroyed to back it up, if necessary.
  ! In a way, this function checks if we are running FFTW or Intel MKL
  logical function backwards_fft_preserves_input()
    type(c_ptr) :: plan
    complex(wp) :: complex_dummy(4)
    real(wp) :: real_dummy(4)
    integer, parameter, dimension(2) :: dummy_size = [2,2]

    !$OMP CRITICAL
    plan = fftw_plan_many_dft_c2r(&
      rank = 2,&
      n = dummy_size,&
      howmany = 1,&
      in = complex_dummy,&
      inembed = dummy_size,&
      istride = 1,&
      idist = 0,&
      out = real_dummy,&
      onembed = dummy_size,&
      ostride = 1,&
      odist = 0,&
      flags = ior(FFTW_EXHAUSTIVE, FFTW_PRESERVE_INPUT))
    !$OMP END CRITICAL

    backwards_fft_preserves_input = C_ASSOCIATED(plan)
    if(backwards_fft_preserves_input) call fftw_destroy_plan(plan)
  end function backwards_fft_preserves_input
  !
  !pure subroutine fourier_cont(&
  !    a, cont_mat_l, cont_mat_r, cont_bc)
  !    real(wp), contiguous, intent(inout) :: a(:)
  !    real(wp), intent(in) :: cont_mat_l(:,:), cont_mat_r(:,:), cont_bc(:)
  !    integer :: i, fc_xstart, fc_xend, fc_dl, fc_dr
  !    !dir$ assume_aligned a:64
  !    
  !    fc_dl = size(cont_mat_l, 2) + 1
  !    fc_dr = size(cont_mat_r, 2) + 1
  !    do i = 1,ny
  !        fc_xstart = (i - 1)*nx + rb_node + 1
  !        fc_xend = fc_xstart + fc_nodes - 1
  !        !fc_xend = (i - 1)*nx + rb_node + fc_nodes
  !        a(fc_xstart:fc_xend) = &
  !            matmul(cont_mat_l, a((i-1)*nx + 2:(i-1)*nx + fc_dl))&
  !            + matmul(cont_mat_r, a(fc_xstart - fc_dr:fc_xstart - 2))&
  !            + cont_bc
  !    end do
  !end subroutine fourier_cont

  ! FIXME: reimplement fourier continuation
  !pure subroutine fourier_cont_var_lb(a, fc_data, surf_x_node)
  !    type(fft_field), intent(inout) :: a
  !    type(fourier_cont_data), intent(in) :: fc_data
  !    real(wp), contiguous, intent(in) :: surf_x_node(:)
  !    real(wp) :: bc_offset
  !    integer :: node_y, fc_xstart, fc_xend, surf_x_nearest_node
  !    
  !    do concurrent(node_y = 1:ny)
  !        surf_x_nearest_node = nint(surf_x_node(node_y))
  !        bc_offset = surf_x_node(node_y) - surf_x_nearest_node
  !        fc_xend = (node_y - 1)*nx + surf_x_nearest_node - 1
  !        fc_xstart = fc_xend - fc_nodes + 1
  !        a%real%vals(fc_xstart:fc_xend) = &
  !            matmul(cont_mat_l_factor(fc_data%d_l, bc_offset)*fc_data%cont_mat_l, a%real%vals(fc_xend + 2:fc_xend + fc_data%d_l))&
  !            + matmul(fc_data%cont_mat_r, a%real%vals(fc_xstart - fc_data%d_r:fc_xstart - 2))&
  !            + cont_bc_l_factor(fc_data%d_l, bc_offset)*fc_data%cont_bc_l*a%bc%val_min + fc_data%cont_bc_r*a%bc%val_max
  !    end do
  !end subroutine fourier_cont_var_lb

  !pure subroutine fc_clear(a, nx_int)
  !    real(wp), intent(inout) :: a(:)
  !    integer, intent(in) :: nx_int
  !    integer :: node_y
  !    
  !    if(nx_int < nx) then
  !        do concurrent(node_y = 1:ny)
  !            a((node_y-1)*nx + nx_int + 1:node_y*nx) = epsilon(0.0_wp)
  !        end do
  !    end if
  !end subroutine fc_clear

  !pure real(wp) elemental function cont_bc_l_factor(fc_dl, bc_offset)
  !    integer, intent(in) :: fc_dl
  !    real(wp), intent(in) :: bc_offset
  !    integer :: i
  !    
  !    cont_bc_l_factor = 1.0_wp
  !    do concurrent(i = 1:fc_dl - 1)
  !        cont_bc_l_factor = cont_bc_l_factor*(1 - bc_offset/i)
  !    end do
  !    cont_bc_l_factor = 1.0_wp/cont_bc_l_factor
  !end function cont_bc_l_factor
  !
  !pure function cont_mat_l_factor(fc_dl, bc_offset)
  !    integer, intent(in) :: fc_dl
  !    real(wp), intent(in) :: bc_offset
  !    real(wp) :: cont_mat_l_factor(fc_nodes, fc_dl - 1)
  !    integer :: i, j
  !    
  !    do concurrent(i = 1:fc_nodes, j = 1:fc_dl - 1)
  !        if(abs(bc_offset) > epsilon(1.0_wp)) then
  !            cont_mat_l_factor(i, j) = 1/(1 - bc_offset/j) + j/((fc_nodes - i + 1)*(j/bc_offset - 1))
  !        else
  !            cont_mat_l_factor(i, j) = 1
  !        end if
  !    end do
  !end function cont_mat_l_factor
end module mod_spectral
