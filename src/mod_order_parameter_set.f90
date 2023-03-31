! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_order_parameter_set
  use mod_base_functions
  use mod_parameters
  use mod_regular_grid
  use mod_time_dep_vars
  use mod_spectral
  use mod_phase_properties
  use mod_phase_interface
  use mod_phase_set
  use mod_initial_field_settings
  use mod_random
  implicit none
  private

  public :: equi_eta

  real(wp), parameter :: eta_zero = epsilon(1.0_wp)

  character(*), parameter :: module_name = 'mod_order_parameter_set'

  type, extends(fft_field), public :: order_parameter_t
    type(phase_interface) :: interface_data
    integer :: idx, order_idx, phase_idx
  contains
    procedure, public :: sharp_to_diffuse_interfaces
    procedure, public :: get_highest_order_idx_active
    procedure, public :: get_n_eta_active
    procedure, public :: is_active
  end type order_parameter_t

  type, public :: order_parameter_set
    integer :: n_eta
    type(regular_grid) :: grid
    type(order_parameter_t), allocatable :: eta
    type(time_dependent_real) :: sqsum, quadsum
    type(time_dependent_integer) :: eta_labels
  contains
    procedure, public :: init
    procedure, public :: set_if_param
    procedure, public :: calc_sqsum
    procedure, public :: calc_quadsum
    procedure, public :: calc_phi
    procedure, public :: check_eta_is_valid
    procedure, public :: apply_initial_conditions
    procedure, public :: eta_to_eta_labels
    procedure, public :: eta_labels_to_eta
    procedure, public :: reset
  end type order_parameter_set

contains

  subroutine init(this, grid, phase_properties, time_now)
  class(order_parameter_set), intent(inout) :: this
    type(regular_grid), intent(in) :: grid
    type(phase_properties_t), intent(in), target :: phase_properties(:)
    logical, allocatable :: corresponding_phase(:)
    real(wp), intent(in) :: time_now
    character(*), parameter :: proc_name = 'init'

    allocate(this%eta)
    this%n_eta = num_images()
    this%eta%idx = this_image()
    this%grid = grid
    corresponding_phase = this%eta%idx >= phase_properties%idx%eta%first .and. this%eta%idx <= phase_properties%idx%eta%last
    if(count(corresponding_phase) == 0) then
      call error_msg(proc_name, 'Phase corresponding to eta_idx ' // convert_to_char(this%eta%idx) // ' not found' &
        , module_name_opt = module_name)
    else if(count(corresponding_phase) > 1) then
      call error_msg(proc_name, 'Multiple phases defined for eta_idx ' // convert_to_char(this%eta%idx) // &
        ' (check phase eta index range)', module_name_opt = module_name)
    else if(count(corresponding_phase) == 1) then
      this%eta%phase_idx = first_true(corresponding_phase)
      this%eta%order_idx = this%eta%idx - phase_properties(this%eta%phase_idx)%idx%eta%first + 1
    else
      call error_msg(proc_name, module_name_opt = module_name)
    end if
    call this%eta%init(&
      grid, eta_name // '_' // trim(phase_properties(this%eta%phase_idx)%name) // '_' &
      // convert_to_char(this%eta%order_idx), keep_prev_fft_opt = .true.)
    call this%eta_labels%init(grid, eta_idx_name)
    this%eta%real%time = time_now
    call this%sqsum%init(grid, sqsum_name)
    call this%quadsum%init(grid, sqsum_name)
    this%sqsum%vals = 0
    this%quadsum%vals = 0
  end subroutine init

  pure elemental subroutine set_if_param(this, interface_data)
  class(order_parameter_set), intent(inout) :: this
    type(phase_interface), intent(in) :: interface_data

    this%eta%interface_data = interface_data
  end subroutine set_if_param

  subroutine calc_sqsum(this)
  class(order_parameter_set), intent(inout) :: this
    character(*), parameter :: proc_name = 'calc_sqsum'
    !dir$ vector aligned

    if(this%sqsum%time_is_equal(this%eta%real)) return
    this%sqsum%vals = this%eta%real%vals**2
    call co_sum(this%sqsum%vals)
    this%sqsum%time = this%eta%real%time
  end subroutine calc_sqsum

  subroutine calc_quadsum(this)
  class(order_parameter_set), intent(inout) :: this
    character(*), parameter :: proc_name = 'calc_quadsum'
    !dir$ vector aligned

    if(this%quadsum%time_is_equal(this%eta%real)) return
    this%quadsum%vals = this%eta%real%vals**4
    call co_sum(this%quadsum%vals)
    this%quadsum%time = this%eta%real%time
  end subroutine calc_quadsum

  subroutine calc_phi(this, phases)
  class(order_parameter_set), intent(inout) :: this
    type(phase_set), intent(inout) :: phases
    character(*), parameter :: proc_name = 'calc_phi'
    !dir$ vector aligned

    if(.not. allocated(this%eta)) call error_msg(proc_name, 'Eta is not allocated', module_name_opt = module_name)
    call this%calc_sqsum()
    if(this%eta%phase_idx < phases%n_phases) then
      associate(phi => phases%phi(this%eta%phase_idx))
        if(phi%time_is_equal(this%eta%real)) return
        phi%vals = this%eta%real%vals**2/this%sqsum%vals
        phi%time = this%eta%real%time
      end associate
    end if
  end subroutine calc_phi

  subroutine check_eta_is_valid(this, caller_name)
  class(order_parameter_set), intent(in) :: this
    character(*), intent(in) :: caller_name

    if(any(this%eta%real%vals < eta_min) .or. any(this%eta%real%vals > eta_max)) then
      call error_msg(caller_name, this%eta%real%var_name &
        // ' outside valid value range (' // convert_to_char(eta_min) // ' - ' // convert_to_char(eta_max) // ')!')
    end if
  end subroutine check_eta_is_valid  

  subroutine sharp_to_diffuse_interfaces(this)
  class(order_parameter_t), intent(inout) :: this
    type(fft_field) :: equi_conv
    real(wp) :: val_sum
    integer :: node

    call equi_conv%init(this%real%grid, 'equi_conv')
    equi_conv%real%time = this%real%time
    do concurrent(node = 1:this%real%grid%n_nodes_total)
      equi_conv%real%vals(node) = equi_deta(this%real%grid%point_distance(1, node), this%interface_data%get_width())
    end do
    val_sum = sum(equi_conv%real%vals)
    do concurrent(node = 1:this%real%grid%n_nodes_total)
      equi_conv%real%vals(node) = equi_conv%real%vals(node)/val_sum
    end do
    call this%forward_alldim()
    call equi_conv%forward_alldim()
    call this%convolve(equi_conv)
    ! Need to invalidate eta, otherwise backward FFT calculation will be skipped
    call this%real%invalidate()
    call this%backward_alldim()
  end subroutine sharp_to_diffuse_interfaces

  pure elemental real(wp) function equi_eta(x, if_width)
    real(wp), intent(in) :: x, if_width

    equi_eta = 0.5_wp*(1 - tanh(2*x/if_width))
  end function equi_eta

  pure elemental real(wp) function equi_phi(x, if_width)
    real(wp), intent(in) :: x, if_width
    real(wp) :: eta_eq

    eta_eq = equi_eta(x, if_width)
    equi_phi = eta_eq**2/(2*eta_eq*(eta_eq - 1) + 1)
  end function equi_phi

  pure elemental real(wp) function equi_deta(x, if_width)
    real(wp), intent(in) :: x, if_width

    equi_deta = (1 - tanh(2*x/if_width)**2)/if_width
  end function equi_deta

  pure elemental real(wp) function equi_dphi(x, if_width)
    real(wp), intent(in) :: x, if_width
    real(wp) :: eta_eq, dphi_deta

    eta_eq = equi_eta(x, if_width)
    dphi_deta = 2*eta_eq*(eta_eq - 1)/(2*eta_eq*(eta_eq - 1) + 1)**2
    equi_dphi = dphi_deta*equi_deta(x,if_width)
  end function equi_dphi

  subroutine apply_initial_conditions(this, init_field)
  class(order_parameter_set), intent(inout) :: this
    type(initial_field_settings), intent(in) :: init_field
    real(wp) :: particle_r, solid_fraction
    logical, allocatable, dimension(:) :: solid
    integer :: i, node, n_solid
    integer, allocatable, dimension(:) :: p, invalid_eta
    real(wp), allocatable, dimension(:) :: center_coordinates, precipitate_coordinates, precipitates
    !type(team_type) :: phase_team
    character(*), parameter :: proc_name = 'apply_initial_conditions'
    real(wp), parameter :: free_coordinate = -1.0_wp

    !form team(this%eta%phase_idx, phase_team)
    call this%eta%real%invalidate()
    this%eta_labels%time = 0.0_wp
    solid = init_field%solid_part%get_solid_mask(this%grid)
    n_solid = count(solid)
    solid_fraction = real(n_solid, wp)/this%grid%n_nodes_total
    allocate(p(this%grid%n_dims))
    allocate(precipitate_coordinates(this%grid%n_dims))
    p = 1
    precipitate_coordinates = free_coordinate
    center_coordinates = (this%grid%axes%n_nodes + 1)/2.0_wp

    particle_r = 0
    if(init_field%n_precipitates > 0) then
      select case(this%grid%n_dims)
      case(1)
        particle_r = n_solid*init_field%precipitate_phase%equi_fraction/2
      case(2)
        particle_r = sqrt(n_solid*init_field%precipitate_phase%equi_fraction/(init_field%n_precipitates*pi))
      case(3)
        particle_r = (0.75_wp*n_solid*init_field%precipitate_phase%equi_fraction &
          /(init_field%n_precipitates*pi))**(1/3.0_wp)
      case default
        call error_msg(proc_name, 'Undefined for given number of dimensions: ' // convert_to_char(this%grid%n_dims))
      end select
    end if

    select case (init_field%eta_type)
    case(init_monocrystal)
      this%eta_labels%vals = init_field%matrix_phase%idx%eta%first
    case(init_polycrystal)
      this%eta_labels%vals = -1
      !change team(phase_team)
      !if(team_number() == init_field%matrix_phase%idx%phase) then
      call init_field%polycrystal%construct(this%grid, init_field%matrix_phase%idx, this%eta_labels%vals)
      where(this%eta_labels%vals == 0)
        this%eta_labels%vals = init_field%pore_phase%idx%eta%first
      end where
      !end if
      !end team
      call co_max(this%eta_labels%vals)
    case(init_single_boundary_pin)
      if(init_field%matrix_phase%n_order < 2) then
        call error_msg(proc_name, 'Initial eta condition "single_boundary_pin" requires at least 2 possible crystal ' &
          // 'orientations in the matrix phase', module_name_opt = module_name)
      end if
      do
        node = this%grid%dim_to_linear_idx(p)
        if(p(1) > center_coordinates(1)) then
          this%eta_labels%vals(node) = init_field%matrix_phase%idx%eta%first
        else
          this%eta_labels%vals(node) = init_field%matrix_phase%idx%eta%first + 1
        end if
        if(all(p == this%grid%axes%n_nodes)) exit
        call this%grid%increment_idx(p)
      end do
      precipitate_coordinates(1) = center_coordinates(1)
    case(init_quad_boundary_pin)
      if(init_field%matrix_phase%n_order < 4) then
        call error_msg(proc_name, 'Initial eta condition "quad_boundary_pin" requires at least 4 possible crystal ' &
          // 'orientations in the matrix phase', module_name_opt = module_name)
      end if
      if(this%grid%n_dims <= 1) then
        call error_msg(proc_name, 'Quadruple boundary is only possible for 2 or more this%grid dimensions')
      end if
      do
        node = this%grid%dim_to_linear_idx(p)
        if(p(1) > center_coordinates(1)) then
          if(p(2) > center_coordinates(2)) then
            this%eta_labels%vals(node) = init_field%matrix_phase%idx%eta%first
          else
            this%eta_labels%vals(node) = init_field%matrix_phase%idx%eta%first + 1
          end if
        else
          if(p(2) > center_coordinates(2)) then
            this%eta_labels%vals(node) = init_field%matrix_phase%idx%eta%first + 2
          else
            this%eta_labels%vals(node) = init_field%matrix_phase%idx%eta%first + 3
          end if
        end if
        if(all(p == this%grid%axes%n_nodes)) exit
        call this%grid%increment_idx(p)
      end do
      precipitate_coordinates(1:2) = center_coordinates(1:2)
    case(init_spinodal)
      if(this_image() == 1) then
        do i = 1, this%grid%n_nodes_total
          this%eta_labels%vals(i) = random_int(1, this%n_eta)
        end do
      end if
      call co_broadcast(this%eta_labels%vals, 1)
    case(init_spinodal_dense)
      if(this_image() == 1) then
        do i = 1, this%grid%n_nodes_total
          this%eta_labels%vals(i) = random_int(2, this%n_eta)
        end do
      end if
      call co_broadcast(this%eta_labels%vals, 1)
    case default
      call error_msg(proc_name, 'Unknown initial eta condition', module_name_opt = module_name)
    end select

    where(.not. solid)
      this%eta_labels%vals = init_field%pore_phase%idx%eta%first
    end where

    invalid_eta = pack(this%eta_labels%vals, (this%eta_labels%vals > this%n_eta) .or. (this%eta_labels%vals <= 0))
    if(size(invalid_eta) > 0) then
      call error_msg(proc_name,&
        'Invalid eta label value encountered after voronoi construction: ' // convert_to_char(invalid_eta))
    end if

    if(all(precipitate_coordinates /= free_coordinate) .and. init_field%n_precipitates > 1) then
      call error_msg(proc_name,&
        'All precipitate coordinates fixed (i.e. single possible location, but more than one precipitate requested')
    end if

    if(this_image() == 1) call this%eta_labels%output_image(name_suffix_opt = 'initial')
    call this%eta_labels_to_eta()
    call this%eta%real%output_image(name_suffix_opt = 'initial_diffuse')

    if(init_field%n_precipitates > 0) then
      if(init_field%n_precipitates > 1) call error_msg(proc_name, 'n_precipitates > 1 currently under construction')
      call this%eta_labels%invalidate()
      allocate(precipitates(this%eta%real%grid%n_nodes_total))
      precipitates = 0.0_wp
      do concurrent(node = 1:this%eta%real%grid%n_nodes_total)
        p = this%eta%real%grid%linear_to_dim_idx(node)
        precipitates(node) = this%eta%real%grid%point_distance(center_coordinates, real(p, wp)) - particle_r
      end do
      precipitates = equi_eta(precipitates, this%eta%interface_data%get_width())
      if(this%eta%idx == init_field%precipitate_phase%idx%eta%first) then
        this%eta%real%vals = min(this%eta%real%vals + precipitates, 1.0_wp)
      else
        this%eta%real%vals = max(this%eta%real%vals - precipitates, 0.0_wp)
      end if
      call this%eta%fft%invalidate()
      call this%eta%real%output_image(name_suffix_opt = 'initial_precipitate')
      call this%check_eta_is_valid(proc_name)
    end if

    if(init_field%spinodal_noise > 0.0_wp) then
      do node = 1, this%eta%real%grid%n_nodes_total
        this%eta%real%vals(node) = this%eta%real%vals(node) &
          + random_real(-init_field%spinodal_noise, init_field%spinodal_noise)
      end do
    end if
  end subroutine apply_initial_conditions

  subroutine eta_to_eta_labels(this)
  class(order_parameter_set), intent(inout) :: this
    !real(wp), parameter :: eta_threshold = 1/3.0_wp

    if(this%eta_labels%time_is_equal(this%eta%real)) return
    this%eta_labels%vals = co_maxloc(this%eta%real%vals)
    !where(this%eta%real%vals > eta_threshold)
    !    this%eta_labels%vals = this%eta%idx
    !elsewhere
    !    this%eta_labels%vals = 0
    !end where
    !call co_max(this%eta_labels%vals)
    this%eta_labels%time = this%eta%real%time
  end subroutine eta_to_eta_labels

  subroutine eta_labels_to_eta(this, eta_shift_opt)
  class(order_parameter_set), intent(inout) :: this
    real(wp), intent(in), optional :: eta_shift_opt(:)
    character(*), parameter :: proc_name = 'eta_labels_to_eta'

    if(this%eta_labels%time_is_equal(this%eta%real)) return
    where(this%eta_labels%vals == this%eta%idx)
      this%eta%real%vals = 1.0_wp
    elsewhere
      this%eta%real%vals = 0.0_wp
    end where
    this%eta%real%time = this%eta_labels%time
    if(present(eta_shift_opt)) call this%eta%real%coordinate_shift(eta_shift_opt)
    call this%eta%sharp_to_diffuse_interfaces()
    call this%check_eta_is_valid(proc_name)
  end subroutine eta_labels_to_eta

  pure logical function is_active(this)
  class(order_parameter_t), intent(in) :: this

    is_active = this%order_idx == 1 .or. any(abs(this%real%vals) > eta_zero)
  end function is_active

  integer function get_n_eta_active(this) result(n_eta_active)
  class(order_parameter_t), intent(in) :: this

    if(this%is_active()) then
      n_eta_active = 1
    else
      n_eta_active = 0
    end if
    call co_sum(n_eta_active)
  end function get_n_eta_active

  integer function get_highest_order_idx_active(this, phase_idx) result(highest_order_idx_active)
  class(order_parameter_t), intent(in) :: this
    integer, intent(in) :: phase_idx

    if(this%is_active() .and. this%phase_idx == phase_idx) then
      highest_order_idx_active = this%order_idx
    else
      highest_order_idx_active = 0
    end if
    call co_max(highest_order_idx_active)
  end function get_highest_order_idx_active

  pure elemental subroutine reset(this)
  class(order_parameter_set), intent(inout) :: this

    call this%sqsum%invalidate()
    call this%quadsum%invalidate()
    call this%eta_labels%invalidate()
    if(allocated(this%eta)) call this%eta%reset()
  end subroutine reset
end module mod_order_parameter_set
