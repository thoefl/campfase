! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_json
  use json_module
  use hdf5, only: H5T_NATIVE_DOUBLE, H5T_NATIVE_REAL
  use mod_globals, only: verbosity_level
  use mod_base_functions
  use mod_axis
  use mod_boundary_condition
  use mod_regular_grid
  use mod_stopwatch
  use mod_solid_shape
  use mod_polycrystal
  use mod_initial_field_settings
  use mod_parameters
  use mod_species_properties
  use mod_phase_properties
  use mod_phase_interface
  use mod_time_step
  use mod_fft_osc
  use mod_spectral, only: fft_grid_from_real
  use mod_nucleation
  use mod_constants
  use mod_unit_factors
  use mod_arrhenius_model
  use mod_hdf5, only: hdf_file
  implicit none
  private

  integer, parameter :: op_invalid = -1
  integer, parameter :: op_unset = 0
  integer, parameter :: op_multiply = 1
  integer, parameter :: op_divide = 2
  integer, parameter :: op_power = 3
  character(*), parameter :: module_name = 'mod_json'

  public :: json_read_grid, json_read_species, json_read_phase_properties, json_read_interfacial_data, json_read_time_step, &
    json_read_fft_oscillation_data, json_read_timer_set_data, json_read_hdf_output_data, json_read_grain_reindex_data, &
    json_read_nucleation_data, json_read_pressure_temperature_data, json_read_chemical_energy_data, &
    json_read_restore_file_settings, json_read_benchmark_settings, json_read_fourier_continuation_settings, &
    json_read_check_settings, json_read_initial_conditions, json_read_opencalphad_settings, json_calculate_units, json_file

contains

  pure elemental integer(int8) function boundary_condition_string_to_int(bc_type_string) result(bc_type_int)
    character(*), intent(in) :: bc_type_string
    character(*), parameter :: proc_name = 'boundary_condition_string_to_int'

    select case(bc_type_string)
    case('periodic')
      bc_type_int = bc_periodic
    case('dirichlet')
      bc_type_int = bc_dirichlet
    case('neumann')
      bc_type_int = bc_neumann
    case default
      bc_type_int = 0
      call error_msg(proc_name, 'Undefined boundary condition type')
    end select
  end function

  subroutine json_calculate_units(json, units)
    type(json_file), intent(inout) :: json
    type(unit_factor_t), intent(out) :: units
    character(:), allocatable :: unit_name
    real(wp) :: spacing
    character(*), parameter :: proc_name = 'json_calculate_units'

    call json%get('grid.axes(1).spacing.unit', unit_name)
    if(unit_name /= 'meter') call error_msg(proc_name, 'Units other than meter for axes spacing currently not implemented')
    call json%get('grid.axes(1).spacing.value', spacing)
    units = unit_factor_t(V_m = 1E-5_wp, G_m = 1E3_wp, D = 1E-14_wp, meter = spacing, kelvin = 1.0_wp)
  end subroutine json_calculate_units

  subroutine json_read_grid(json, constants, units, grid)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(regular_grid), intent(out) :: grid
    integer :: n_dims, dim, n_nodes
    type(axis), allocatable :: axes(:)
    character(:), allocatable :: bc_type_string, label
    real(wp) :: bc_val_min, bc_val_max, spacing
    integer(int8) :: bc_type_min, bc_type_max
    type(boundary_condition) :: bc
    character(*), parameter :: proc_name = 'json_read_grid'

    call json%info('grid.axes', n_children = n_dims)

    allocate(axes(n_dims))

    do dim = 1, n_dims
      call json%get('grid.axes(' // convert_to_char(dim) // ').nodes', n_nodes)
      call json_get_dimensionless_value(json, 'grid.axes(' // convert_to_char(dim) // ').spacing', constants, units, spacing)
      call json%get('grid.axes(' // convert_to_char(dim) // ').label', label)
      call json%get('grid.axes(' // convert_to_char(dim) // ').boundary_condition.minimum.value', bc_val_min)
      call json%get('grid.axes(' // convert_to_char(dim) // ').boundary_condition.maximum.value', bc_val_max)
      call json%get('grid.axes(' // convert_to_char(dim) // ').boundary_condition.minimum.type', bc_type_string)
      bc_type_min = boundary_condition_string_to_int(bc_type_string)
      call json%get('grid.axes(' // convert_to_char(dim) // ').boundary_condition.maximum.type', bc_type_string)
      bc_type_max = boundary_condition_string_to_int(bc_type_string)
      bc = boundary_condition(bc_type_min, bc_type_max, bc_val_min, bc_val_max)
      axes(dim) = axis(n_nodes, spacing, bc, label)
    end do
    call grid%init(axes)

  end subroutine json_read_grid

  subroutine json_read_species(json, constants, units, species, n_species)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(species_properties_t), allocatable, intent(out) :: species(:)
    character(:), allocatable :: species_root, species_name_fetch
    integer, intent(out) :: n_species
    real(wp) :: A, E_a
    integer :: species_idx
    character(*), parameter :: proc_name = 'json_read_species'

    call json%info('species', n_children = n_species)
    allocate(species(n_species))
    do species_idx = 1, n_species
      species_root = 'species(' // convert_to_char(species_idx) // ')'
      call json%get(species_root // '.name', species_name_fetch)
      species(species_idx)%name = species_name_fetch
      call json%get(species_root // '.instant_diffusion', species(species_idx)%instant_diffusion)
      call json_get_dimensionless_value(json, species_root // '.diffusion_coefficient.A', constants, units, A)
      call json_get_dimensionless_value(json, species_root // '.diffusion_coefficient.E_a', constants, units, E_a)
      species(species_idx)%diffusion_coefficient = arrhenius_model(A, E_a, constants%R_gas)
    end do
  end subroutine json_read_species

  subroutine json_read_phase_properties(json, constants, units, phases, n_phases)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(phase_properties_t), allocatable, intent(out) :: phases(:)
    character(:), allocatable :: phase_root, string_fetch
    integer, intent(out) :: n_phases
    integer :: phase_idx
    character(*), parameter :: proc_name = 'json_read_phase_properties'

    call json%info('phases', n_children = n_phases)
    allocate(phases(n_phases))
    do phase_idx = 1, n_phases
      phase_root = 'phases(' // convert_to_char(phase_idx) // ')'
      call json%get(phase_root // '.name', string_fetch)
      phases(phase_idx)%name = string_fetch
      call json_get_dimensionless_value(json,&
        phase_root // '.V_m', constants, units, phases(phase_idx)%V_m)
      call json_get_dimensionless_value(&
        json, phase_root // '.mobility_factor', constants, units, phases(phase_idx)%mobility_factor)
      call json%get(phase_root // '.is_compset', phases(phase_idx)%is_compset)
      call json%get(phase_root // '.calphad_use_phase_conc', phases(phase_idx)%calphad_use_phase_conc)
      call json%get(phase_root // '.default_comp', string_fetch)
      phases(phase_idx)%default_comp = string_fetch
      call json%get(phase_root // '.brightness', phases(phase_idx)%brightness)
    end do
  end subroutine json_read_phase_properties

  subroutine json_read_interfacial_data(json, constants, units, interfacial_data)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(phase_interface), intent(out) :: interfacial_data
    real(wp) :: width, energy, gamma
    character(*), parameter :: proc_name = 'json_read_interfacial_data'


    call json_get_dimensionless_value(json, 'interface.width', constants, units, width)
    call json_get_dimensionless_value(json, 'interface.energy', constants, units, energy)
    call json%get('interface.gamma', gamma)
    call interfacial_data%init(width, energy, gamma)
  end subroutine json_read_interfacial_data

  subroutine json_read_time_step(json, constants, units, sim_clock, wall_clock, dt)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(system_stopwatch), intent(in), target :: wall_clock
    type(manual_stopwatch), intent(in), target :: sim_clock
    type(time_step), intent(out) :: dt
    real(wp) :: bdfab_switch_min_factor
    character(*), parameter :: proc_name = 'json_read_time_step'


    call json_get_dimensionless_value(json, 'time_step.initial', constants, units, dt%init_val)
    dt%val = dt%init_val
    call json%get('time_step.prestep_c', dt%prestep_c)
    call json%get('time_step.bdfab.enabled', dt%bdfab_enabled)
    call json%get('time_step.bdfab.force', dt%force_bdfab)
    dt%use_bdfab = dt%bdfab_enabled
    call json%get('time_step.adaptive.enabled', dt%adaptive%enabled)
    if(dt%adaptive%enabled) then
      call json_read_timer_set_data(json, constants, units, sim_clock, wall_clock, 'time_step.adaptive', dt%adaptive%timer)
      call json%get('time_step.adaptive.rel_error_tolerance', dt%adaptive%rel_err_tol)
      call json%get('time_step.adaptive.max_change_factor', dt%adaptive%max_factor)
      call json%get('time_step.adaptive.min_change_factor', dt%adaptive%min_factor)
      call json%get('time_step.adaptive.bdfab_min_factor_for_switch_to_euler', bdfab_switch_min_factor)
      dt%adaptive%bdfab_switch_rel_error_max = dt%adaptive%rel_err_tol/bdfab_switch_min_factor
      call json%get('time_step.adaptive.richardson_steps', dt%adaptive%n_steps_richardson)
      call json%get('time_step.adaptive.use_bdfab', dt%adaptive%use_bdfab)
      if(dt%adaptive%min_factor > 1.0_wp) call error_msg(proc_name, 'time_step.adaptive.min_change_factor must be <= 1.0')
      if(dt%adaptive%max_factor < 1.0_wp) call error_msg(proc_name, 'time_step.adaptive.max_change_factor must be >= 1.0')
      dt%adaptive%rel_err_max = dt%adaptive%rel_err_tol/dt%adaptive%min_factor
    end if
  end subroutine json_read_time_step

  !FIXME: Data should ideally not be allocated for variables that are not checked
  subroutine json_read_fft_oscillation_data(json, real_grid, var_name, osc)
    type(json_file), intent(inout) :: json
    type(regular_grid), intent(in) :: real_grid
    type(regular_grid) :: fft_grid
    character(*), intent(in) :: var_name
    type(fft_grid_max_freq_oscillations), allocatable, intent(inout) :: osc
    real(wp) :: errtol
    integer :: error_nsave, n_vars, var, n_osc_max
    character(:), allocatable :: json_var_name
    character(*), parameter :: proc_name = 'json_read_fft_oscillation_data'

    call json%info('fft_oscillation_check.variables', n_children = n_vars)
    do var = 1, n_vars
      call json%get('fft_oscillation_check.variables(' // convert_to_char(var) // ')', json_var_name)
      if(wildcard_string_match(var_name, json_var_name)) then
        if(.not. allocated(osc)) allocate(osc)
        fft_grid = fft_grid_from_real(real_grid)
        call json%get('fft_oscillation_check.tolerance', errtol)
        call json%get('fft_oscillation_check.max_count', n_osc_max)
        call json%get('fft_oscillation_check.n_save', error_nsave)
        call osc%init(var_name, real_grid, fft_grid, errtol, n_osc_max, error_nsave)
        exit
      end if
    end do
  end subroutine json_read_fft_oscillation_data

  !FIXME: dstep_trigger should be int64, but json-fortran built with default kind
  !       (should not be a problem in most cases)
  subroutine json_read_timer_set_data(json, constants, units, sim_clock, wall_clock, set_name, set, late_allowed_opt)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(system_stopwatch), intent(in), target :: wall_clock
    type(manual_stopwatch), intent(in), target :: sim_clock
    character(*), intent(in) :: set_name
    type(timer_set), intent(inout) :: set
    logical, intent(in), optional :: late_allowed_opt
    real(wp) :: dt_trigger
    integer :: dstep_trigger, n_timers, i
    character(:), allocatable :: timer_name
    character(*), parameter :: proc_name = 'json_read_timer_set_data'

    call json%info(set_name // '.by', n_children=n_timers)
    do i = 1, n_timers
      call json%info(set_name // '.by(' // convert_to_char(i) // ')', name=timer_name)
      select case(timer_name)
      case('wall_clock')
        allocate(set%wall_timer)
        call json_get_dimensionless_value(json, set_name // '.by.wall_clock', constants, units, dt_trigger)
        ! Restore dimension for wall time
        dt_trigger = dt_trigger*units%second
        call set%wall_timer%init(wall_clock, dt_trigger)
      case('sim_clock')
        allocate(set%sim_timer)
        call json_get_dimensionless_value(json, set_name // '.by.sim_clock', constants, units, dt_trigger)
        call set%sim_timer%init(sim_clock, dt_trigger)
      case('step')
        allocate(set%step_timer)
        call json%get(set_name // '.by.step.count', dstep_trigger)
        call set%step_timer%init(int(dstep_trigger, int64), late_allowed_opt)
      case default
        call error_msg(proc_name, 'Unknown timer type in ' // set_name // ': ' // timer_name)
      end select
    end do
  end subroutine json_read_timer_set_data

  subroutine json_read_hdf_output_data(json, constants, units, grid, start_date, start_time, sim_clock, wall_clock, hdf)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(regular_grid), intent(in) :: grid
    character(8), intent(in) :: start_date
    character(10), intent(in) :: start_time
    type(system_stopwatch), intent(in), target :: wall_clock
    type(manual_stopwatch), intent(in), target :: sim_clock
    type(hdf_file), intent(out) :: hdf
    character(:), allocatable :: sim_name, invalid_chars
    character(*), parameter :: hdf_output_root = 'output.hdf'
    character(*), parameter :: proc_name = 'json_read_hdf_output_data'

    ! Spatial plus one temporal dimension
    hdf%n_dims = grid%n_dims + 1
    hdf%maxdims = [maxval(grid%axes%n_nodes), grid%axes(2:)%n_nodes, hdf_max_saves]
    hdf%startdims = [grid%axes%n_nodes, 1]
    hdf%arraychunk = hdf%startdims
    allocate(hdf%scalarchunk(grid%n_dims))
    hdf%scalarchunk = 1

    call json_read_timer_set_data(json, constants, units, sim_clock, wall_clock, hdf_output_root, hdf%output)
    call json_read_timer_set_data(&
      json, constants, units, sim_clock, wall_clock, hdf_output_root // '.eta_backup', hdf%eta_output)

    if(kind(1.0_wp) == 8) then
      hdf%realtype = H5T_NATIVE_DOUBLE
    else if (kind(1.0_wp) == 4) then
      hdf%realtype = H5T_NATIVE_REAL
    else
      call error_msg(proc_name, 'Unrecognized working precision (wp) for type real, '&
        // 'should be either sp (real32) or dp (real64)')
    end if

    call json%get('simulation_name', sim_name)
    if(verify(trim(sim_name), valid_filename_chars) /= 0) then
      invalid_chars = get_chars_not_in_list(sim_name, valid_filename_chars)
      call error_msg(proc_name, 'simulation_name contains invalid character(s): ' // invalid_chars)
    end if
    if(len_trim(sim_name) > 0) then
      hdf%filename = start_date // '_' // start_time(1:6) // '_' // trim(sim_name) // '.h5'
    else
      hdf%filename = start_date // '_' // start_time(1:6) // '.h5'
    end if
    call json%get(hdf_output_root // '.extended', hdf%extended)
    call json%get(hdf_output_root // '.compression_level', hdf%compression_level)
    call json%get(hdf_output_root // '.save_on_dt_min_factor', hdf%save_on_dt_min_factor)
  end subroutine json_read_hdf_output_data

  subroutine json_read_grain_reindex_data(json, constants, units, sim_clock, wall_clock, reindex_data)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(system_stopwatch), intent(in), target :: wall_clock
    type(manual_stopwatch), intent(in), target :: sim_clock
    type(grain_reindex_settings) :: reindex_data
    character(*), parameter :: json_root = 'grain_reindex'

    call json_read_timer_set_data(json, constants, units, sim_clock, wall_clock, json_root, reindex_data%timer)
    call json%get(json_root // '.for_performance.enabled', reindex_data%for_performance%enabled)
    call json%get(json_root // '.for_performance.n_eta', reindex_data%for_performance%n_eta)
    call json%get(json_root // '.for_performance.n_eta_min_active', reindex_data%for_performance%n_eta_min_active)
    call json%get(json_root // '.for_performance.fill_ratio', reindex_data%for_performance%fill_ratio)
    call json%get(json_root // '.eta_separators.intra', reindex_data%separators%intra)
    call json%get(json_root // '.eta_separators.inter', reindex_data%separators%inter)
    call json%get(json_root // '.eta_separators.fill', reindex_data%separators%fill)
    call json%get(json_root // '.eta_separators.relative_tolerance', reindex_data%separators%rel_tol)
  end subroutine json_read_grain_reindex_data

  subroutine json_read_nucleation_data(json, constants, units, sim_clock, wall_clock, nucleation)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(system_stopwatch), intent(in), target :: wall_clock
    type(manual_stopwatch), intent(in), target :: sim_clock
    type(nucleation_data), intent(out) :: nucleation
    integer :: n_children, i
    character(:), allocatable :: rate_name
    character(*), parameter :: json_root = 'nucleation'
    character(*), parameter :: proc_name = 'json_read_nucleation_data'

    call json_read_timer_set_data(json, constants, units, sim_clock, wall_clock, json_root, nucleation%timer)
    call json_read_timer_set_data(json, constants, units, sim_clock, wall_clock, json_root // '.seed_hold', &
      nucleation%expiration_template, late_allowed_opt = .true.)
    call json_get_dimensionless_value(json, json_root // '.radius', constants, units, nucleation%r_nodes)
    call json%get(json_root // '.n_nuclei_max', nucleation%n_nuclei_max)
    call json%info(json_root // '.manual_rate', n_children = n_children)
    do i = 1, n_children
      call json%info(json_root // '.manual_rate(' // convert_to_char(i) // ')', name=rate_name)
      select case(rate_name)
      case('A')
        allocate(nucleation%manual_prefactor)
        call json_get_dimensionless_value(json, json_root // '.manual_rate.' // rate_name, constants, units, &
          nucleation%manual_prefactor)
      case('E_a_factor')
        allocate(nucleation%manual_activation_energy)
        call json%get(json_root // '.manual_rate.' // rate_name, nucleation%manual_activation_energy)
      case default
        call error_msg(proc_name, &
          'Unrecognized factor in ' // json_root // '.manual_rate: ' // rate_name)
      end select
    end do
  end subroutine json_read_nucleation_data

  subroutine json_read_pressure_temperature_data(json, constants, units, temperature, pressure)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(piecewise_linear), intent(out) :: temperature
    real(wp), intent(out) :: pressure
    integer :: n_profile_steps, i
    real(wp), allocatable, dimension(:) :: temp_val, time_val
    real(wp) :: temp_rate, time_hold, temp_init
    character(:), allocatable :: step_type, step_root
    character(*), parameter :: proc_name = 'json_read_pressure_temperature_data'

    call json_get_dimensionless_value(json, 'initial_conditions.pressure', constants, units, pressure)
    call json_get_dimensionless_value(json, 'initial_conditions.temperature', constants, units, temp_init)
    call json%info('temperature_profile', n_children = n_profile_steps)
    if(n_profile_steps == 0) then
      temp_val = [temp_init, temp_init]
      time_val = [0.0_wp, huge(0.0_wp)]
    else
      allocate(temp_val(n_profile_steps + 2))
      allocate(time_val(n_profile_steps + 2))
      temp_val(1) = temp_init
      time_val(1) = 0.0_wp
      do i = 1, n_profile_steps
        step_root = 'temperature_profile(' // convert_to_char(i) // ')'
        call json%get(step_root // '.type', step_type)
        if(step_type == 'ramp') then
          call json_get_dimensionless_value(json, step_root // '.target', constants, units, temp_val(i + 1))
          if(temp_val(i + 1) == temp_val(i)) then
            call error_msg(proc_name, 'Ramp target temperature is equal to temperature given in previous step')
          end if
          call json_get_dimensionless_value(json, step_root // '.rate', constants, units, temp_rate)
          if(.not. sign_is_equal(temp_val(i + 1) - temp_val(i), temp_rate)) then
            call error_msg(proc_name, 'Ramp target temperature is inconsistent with sign of rate given')
          end if
          time_val(i + 1) = time_val(i) + (temp_val(i + 1) - temp_val(i))/temp_rate
        else if(step_type == 'hold') then
          call json_get_dimensionless_value(json, step_root // '.duration', constants, units, time_hold)
          if(time_hold <= 0.0_wp) call error_msg(proc_name, 'Hold time has to be greater than 0')
          temp_val(i + 1) = temp_val(i)
          time_val(i + 1) = time_val(i) + time_hold
        else
          call error_msg(proc_name, 'Unknown temperature profile step type ' // step_type)
        end if
      end do
      temp_val(size(temp_val)) = temp_val(size(temp_val) - 1)
      time_val(size(time_val)) = huge(0.0_wp)
    end if
    call temperature%init(time_val, temp_val)
  end subroutine json_read_pressure_temperature_data

  subroutine json_read_chemical_energy_data(json, constants, units, chem)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(chemical_energy_model), intent(out) :: chem
    character(*), parameter :: json_root = 'chemical_energy'

    call json_get_dimensionless_value(json, json_root // '.parabolic_coefficient.minimum', constants, units, chem%p_coeff_min)
    call json_get_dimensionless_value(json, json_root // '.parabolic_coefficient.maximum', constants, units, chem%p_coeff_max)
  end subroutine json_read_chemical_energy_data

  subroutine json_read_restore_file_settings(json, restore)
    type(json_file), intent(inout) :: json
    type(restore_file_settings), intent(out) :: restore
    character(*), parameter :: json_root = 'restore_file'

    call json%get(json_root // '.enabled', restore%enabled)
    if(restore%enabled) then
      call json%get(json_root // '.c_only', restore%c_only)
      call json%get(json_root // '.restore_idx', restore%idx)
      call json%get(json_root // '.hdf_filename', restore%hdf_filename)
    end if
  end subroutine json_read_restore_file_settings

  subroutine json_read_benchmark_settings(json, benchmark_name, bench)
    type(json_file), intent(inout) :: json
    character(*), intent(in) :: benchmark_name
    type(benchmark_settings), intent(out) :: bench
    character(*), parameter :: json_root = 'benchmark'

    call json%get(json_root // '.' // benchmark_name // '.enabled', bench%enabled)
    call json%get(json_root // '.' // benchmark_name // '.count', bench%count)
  end subroutine json_read_benchmark_settings

  subroutine json_read_fourier_continuation_settings(json, f_cont)
    type(json_file), intent(inout) :: json
    type(fourier_continuation_settings), intent(out) :: f_cont
    character(*), parameter :: json_root = 'fourier_continuation'

    call json%get(json_root // '.n_nodes', f_cont%n_nodes)
    call json%get(json_root // '.w', f_cont%w)
    call json%get(json_root // '.lazy_boundary', f_cont%lazy_boundary)
  end subroutine json_read_fourier_continuation_settings

  subroutine json_read_check_settings(json, check)
    type(json_file), intent(inout) :: json
    type(check_settings), intent(out) :: check
    character(*), parameter :: json_root = 'check'

    call json%get(json_root // '.debug_level', check%debug_level)
    call json%get(json_root // '.equilibrium', check%equilibrium)
    call json%get(json_root // '.energy_model', check%energy_model)
  end subroutine json_read_check_settings

  subroutine json_read_initial_conditions(json, constants, units, grid, species, phases, init_field)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(regular_grid), intent(in) :: grid
    type(species_properties_t), intent(inout) :: species(:)
    type(phase_properties_t), intent(in), target :: phases(:)
  class(initial_field_settings), intent(out) :: init_field
    integer :: n_eta_shift_dims, n_finite_dims, dim, finite_dim, balance_idx, i, species_idx, n_x_init, n_unset, phase_idx
    character(:), allocatable :: type_string, shape_string, axis_label, x_init_name, phase_name
    character(*), parameter :: proc_name = 'json_read_initial_conditions'
    real(wp), parameter :: unset = -huge(0.0_wp)

    balance_idx = 0
    species%x_init = unset
    call json%info('initial_conditions.atomic_fractions', n_children = n_x_init)
    if(n_x_init < size(species) - 1) then
      call error_msg(proc_name, &
        'Initial atomic fractions have to be defined for at least n_species - 1 species (underdetermined)')
    else if(n_x_init > size(species)) then
      call error_msg(proc_name, &
        'Number of initial atomic fractions in "initial_conditions" greater than amount of species defined in "species"')
    end if
    do i = 1, n_x_init
      call json%info('initial_conditions.atomic_fractions(' // convert_to_char(i) // ')', name=x_init_name)
      do species_idx = 1, size(species)
        if(trim(species(species_idx)%name) == trim(x_init_name)) then
          call json%get('initial_conditions.atomic_fractions.' // x_init_name, species(species_idx)%x_init)
          exit
        else if(species_idx == size(species)) then
          call error_msg(proc_name, 'Species defined in "initial_conditions" not found in "species" definitions: ' &
            // trim(x_init_name))
        end if
      end do
    end do
    n_unset = count(species%x_init == unset)
    if(n_unset > 1) then
      call error_msg(proc_name, 'Number of unset species concentrations is greater than 1')
    end if
    balance_idx = first_true(species%x_init == unset)
    if(balance_idx > 0) then
      species(balance_idx)%x_init = 0.0_wp
      species(balance_idx)%x_init = 1 - sum(species%x_init)
      if(species(balance_idx)%x_init < 0.0_wp) then
        call error_msg(proc_name, 'All species initial atomic fractions have to be positive (balance species < 0)')
      end if
    end if
    if(abs(sum(species%x_init) - 1) > epsilon(1.0_wp)) then
      call error_msg(proc_name, 'Species initial atomic fractions have to sum up to 1')
    end if
    if(any(species%x_init < 0.0_wp)) then
      call error_msg(proc_name, 'All species initial atomic fractions have to be positive')
    end if

    ! FIXME: Enable free choice/location of pore phase?

    init_field%pore_phase => phases(1)

    call json%get('initial_conditions.solid.matrix_phase', phase_name)
    phase_idx = first_true(phases%name == phase_name)
    if(phase_idx == 0) call error_msg(proc_name, 'Matrix phase "' // phase_name &
      // '" defined in initial_conditions.solid not found in phase definitions')
    init_field%matrix_phase => phases(phase_idx)

    call json%get('initial_conditions.solid.precipitate_phase', phase_name)
    phase_idx = first_true(phases%name == phase_name)
    if(phase_idx == 0) call error_msg(proc_name, 'Precipitate phase "' // phase_name &
      // '"defined in initial_conditions.solid not found in phase definitions')
    init_field%precipitate_phase => phases(phase_idx)

    call json%get('initial_conditions.solid.precipitates.count', init_field%n_precipitates)

    allocate(init_field%solid_part%axes(grid%n_dims))
    init_field%solid_part%axes%shape = shape_infinite
    init_field%solid_part%axes%length = huge(0.0_wp)
    call json%info('initial_conditions.solid.finite_dimensions', n_children = n_finite_dims)
    do finite_dim = 1, n_finite_dims
      call json%get('initial_conditions.solid.finite_dimensions(' &
        // convert_to_char(finite_dim) // ').axis_label', axis_label)
      do dim = 1, grid%n_dims
        if(grid%axes(dim)%label == axis_label) then
          if(init_field%solid_part%axes(dim)%shape == shape_infinite) then
            call json%get('initial_conditions.solid.finite_dimensions(' &
              // convert_to_char(finite_dim) // ').surface', shape_string)
            select case(shape_string)
            case('plane')
              init_field%solid_part%axes(dim)%shape = shape_plane
            case('round')
              init_field%solid_part%axes(dim)%shape = shape_round
            case('default')
              call error_msg(proc_name, &
                'Unrecognized axis shape in initial_conditions.solid.shape: ' // shape_string)
            end select
            call json_get_dimensionless_value(json, 'initial_conditions.solid.finite_dimensions(' &
              // convert_to_char(finite_dim) // ').extent', constants, units, &
              init_field%solid_part%axes(dim)%length, rel_value_opt = grid%axes(dim)%length)
          else
            call error_msg(proc_name, &
              'Axis multiply defined in "initial_conditions.solid.finite_dimensions":' // axis_label)
          end if
          exit
        else if(dim == grid%n_dims) then
          call error_msg(proc_name, &
            'Axis label defined in "initial_conditions.solid.finite_dimensions" not found in grid axes: ' &
            // axis_label)
        end if
      end do
    end do

    if(count(init_field%solid_part%axes%shape == shape_round) == 1) then
      call error_msg(proc_name, 'Axis shape cannot be "round" in exactly one dimension')
    end if

    call json%info('initial_conditions.eta_shift.value', n_children = n_eta_shift_dims)
    if(n_eta_shift_dims /= grid%n_dims) then
      call error_msg(proc_name,&
        'Number of dimensions in "initial_conditions.eta_shift" is different from number of axes defined in "grid.axes"')
    end if
    allocate(init_field%eta_shift(n_eta_shift_dims))
    do dim = 1, n_eta_shift_dims
      call json_get_dimensionless_value(json, 'initial_conditions.eta_shift', &
        constants, units, init_field%eta_shift(dim), index_opt = dim)
    end do
    call json%get('initial_conditions.spinodal_noise', init_field%spinodal_noise)
    call json%get('initial_conditions.eta', type_string)
    select case(type_string)
    case('monocrystal')
      init_field%eta_type = init_monocrystal
    case('polycrystal')
      init_field%eta_type = init_polycrystal
      allocate(init_field%polycrystal)
      call read_polycrystal_settings(json, constants, units, init_field%polycrystal)
    case('single_boundary_pin')
      init_field%eta_type = init_single_boundary_pin
    case('boundary_check')
      init_field%eta_type = init_boundary_check
    case('spheres')
      init_field%eta_type = init_spheres
    case('planar_x')
      init_field%eta_type = init_planar_x
    case('planar_y')
      init_field%eta_type = init_planar_y
    case('planar_xy')
      init_field%eta_type = init_planar_xy
    case('quad_boundary_pin')
      init_field%eta_type = init_quad_boundary_pin
    case('spinodal')
      init_field%eta_type = init_spinodal
    case('spinodal_dense')
      init_field%eta_type = init_spinodal_dense
    case('sine')
      init_field%eta_type = init_sine
    case default
      call error_msg(proc_name, 'Undefined eta initial condition: ' // type_string)
    end select
  end subroutine json_read_initial_conditions

  subroutine read_polycrystal_settings(json, constants, units, polycrystal)
    type(json_file), intent(inout) :: json
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    type(polycrystal_t), intent(out) :: polycrystal
    character(*), parameter :: json_root = 'initial_conditions.solid.grain_structure'

    call json%get(json_root // '.voronoi.attempts', polycrystal%voronoi%attempts)
    call json%get(json_root // '.voronoi.free_fill_ratio', polycrystal%voronoi%free_fill_ratio)
    call json%get(json_root // '.relative_density', polycrystal%rel_density)
    call json%get(json_root // '.grain_radius_distribution.sigma', polycrystal%radius_distribution%sigma)
    call json_get_dimensionless_value(json,&
      json_root // '.grain_radius_distribution.median', constants, units, polycrystal%radius_distribution%median)
  end subroutine read_polycrystal_settings

  subroutine json_read_opencalphad_settings(json, oc_settings)
    type(json_file), intent(inout) :: json
    type(opencalphad_settings), intent(out) :: oc_settings
    character(*), parameter :: json_root = 'opencalphad'
    character(*), parameter :: proc_name = 'json_read_opencalphad_settings'
    logical :: tdb_file_exists

    call json%get(json_root // '.tdb_filename', oc_settings%tdb_filename)
    inquire(file=oc_settings%tdb_filename, exist=tdb_file_exists)
    if(.not. tdb_file_exists .and. .not. string_contains_character(oc_settings%tdb_filename, '.')) then
      oc_settings%tdb_filename = oc_settings%tdb_filename // '.tdb'
      inquire(file=oc_settings%tdb_filename, exist=tdb_file_exists)
    end if
    if(.not. tdb_file_exists) call error_msg(proc_name, 'TDB file provided not found: ' // oc_settings%tdb_filename)
  end subroutine json_read_opencalphad_settings

  subroutine json_get_dimensionless_value(json, path, constants, units, value, index_opt, rel_value_opt)
    type(json_file), intent(inout) :: json
    character(*), intent(in) :: path
    type(constants_t), intent(in) :: constants
    type(unit_factor_t), intent(in) :: units
    integer, intent(in), optional :: index_opt
    real(wp), intent(in), optional :: rel_value_opt
    real(wp), intent(out) :: value
    character(:), allocatable :: unit_name
    logical :: found
    character(*), parameter :: proc_name = 'json_get_dimensionless_value'

    if(present(index_opt)) then
      call json%get(path // '.value(' // convert_to_char(index_opt) // ')', value, found)
    else
      call json%get(path // '.value', value, found)
    end if
    if(.not. found) call error_msg(proc_name, 'Value for "' // path // '" not found')
    call json%get(path // '.unit', unit_name)
    if(unit_name == 'celsius' .or. unit_name == 'Celsius') then
      value = value + constants%zero_celsius
    end if
    if(unit_name == 'percent' .or. unit_name == 'Percent') then
      if(present(rel_value_opt)) then
        value = 0.01_wp*value*rel_value_opt
      else
        call error_msg(proc_name, 'Relative value (percent) requested, but not defined for given variable: ' // path)
      end if
    else
      if(verbosity_level >= 3 .and. this_image() == 1) then
        print list_value_format, path, value, unit_name, get_unit_expression_factor(unit_name, units), &
        'dimless', value/get_unit_expression_factor(unit_name, units)
      end if
      value = value/get_unit_expression_factor(unit_name, units)
    end if
  end subroutine json_get_dimensionless_value

  pure elemental real(wp) recursive function get_unit_expression_factor(unit_name, units) result(unit_value)
    character(*), intent(in) :: unit_name
    type(unit_factor_t), intent(in) :: units
    integer :: idx, number_start, variable_start, operation, parens_open, parens_close, parens_depth, parens_idx, stat
    logical :: operation_execute
    real(wp) :: operation_value
    character(:), allocatable :: number_string
    character(*), parameter :: proc_name = 'get_unit_expression_factor'
  
    unit_value = 0.0_wp
    if(len(unit_name) == 0) return
    idx = 1
    parens_open = 0
    parens_depth = 0
    number_start = 0
    variable_start = 0
    unit_value = 1.0_wp
    operation = op_multiply
    operation_execute = .false.
    do
      select case(unit_name(idx:idx))
      case('*')
        if(idx == 1) call error_msg(proc_name, 'Operator at initial position of unit name: ' // unit_name)
        if(variable_start > 0 .or. number_start > 0) then
          ! FIXME: exponentiation currently only works if inside parentheses
          !if(idx < len(unit_name)) then
          !    if(unit_name(idx:idx+1)) == '**' then
          !    
          !    end if
          !else
          idx = idx - 1
          operation_execute = .true.
          !end if
        else if(operation == op_unset) then
          operation = op_multiply
        else if(operation == op_multiply) then
          operation = op_power
        else
          operation = op_invalid
        end if
      case('/')
        if(idx == 1) call error_msg(proc_name, 'Operator at initial position of unit name: ' // unit_name)
        if(variable_start > 0 .or. number_start > 0) then
          idx = idx - 1
          operation_execute = .true.
        else if(operation == op_unset) then
          operation = op_divide
        else
          operation = op_invalid
        end if
      case('(')
        parens_open = idx
      case('0':'9')
        if(variable_start > 0) call error_msg(proc_name, 'Number directly after variable expression name: ' // unit_name)
        if(number_start == 0) number_start = idx
      case('a':'z','A':'Z')
        if(number_start > 0) call error_msg(proc_name, 'Variable expression name directly after number: ' // unit_name)
        if(variable_start == 0) variable_start = idx
      case default
        call error_msg(proc_name, 'Invalid character in unit name ("' // unit_name // '") : ' // unit_name(idx:idx))
      end select
      if(operation == op_invalid) then
        call error_msg(proc_name, 'Invalid operation in unit name: ' // unit_name // '(allowed: * / **)')
      end if
      if(parens_open > 0) then
        parens_depth = 1
        parens_close = 0
        do parens_idx = parens_open + 1, len(unit_name)
          select case(unit_name(parens_idx:parens_idx))
          case('(')
            parens_depth = parens_depth + 1
          case(')')
            parens_depth = parens_depth - 1
            if(parens_depth == 0) then
              parens_close = parens_idx
              exit
            end if
          end select
        end do
        if(parens_close == parens_open + 1) call error_msg(proc_name, 'Invalid parenthesis in unit given: ' // unit_name)
        if(parens_close <= parens_open) call error_msg(proc_name, 'Unclosed parenthesis in unit given: ' // unit_name)
        operation_value = get_unit_expression_factor(unit_name(parens_open + 1:parens_close - 1), units)
        unit_value = apply_operation(unit_value, operation_value, operation)
        operation = op_unset
        idx = parens_close
        parens_open = 0
      else if(operation_execute .or. idx == len(unit_name)) then
        if(number_start > 0) then
          allocate(character(idx - number_start + 1) :: number_string)
          number_string = unit_name(number_start:idx)
          read(number_string, *, iostat = stat) operation_value
          if(stat /= 0) call error_msg(proc_name, 'Failed to convert number "' // number_string // '" to type real')
          number_start = 0
          deallocate(number_string)
        else if(variable_start > 0) then
          operation_value = get_unit_value(unit_name(variable_start:idx), units)
          variable_start = 0
        else
          call error_msg(proc_name, 'No value for operation: ' // unit_name)
        end if
        unit_value = apply_operation(unit_value, operation_value, operation)
        operation = op_unset
        operation_execute = .false.
      end if
      if(idx == len(unit_name)) exit
      idx = idx + 1
    end do
  end function get_unit_expression_factor
  
  pure elemental real(wp) function apply_operation(a, b, operation) result(res)
    real(wp), intent(in) :: a, b
    integer, intent(in) :: operation
    character(*), parameter :: proc_name = 'apply_operation'
  
    res = 0.0_wp
    select case(operation)
    case(op_multiply)
      res = a*b
    case(op_divide)
      res = a/b
    case(op_power)
      res = a**b
    case(op_invalid)
      call error_msg(proc_name, 'Invalid operator value given')
    case(op_unset)
      call error_msg(proc_name, 'Operator value unset')
    case default
      call error_msg(proc_name, 'Undefined operator value given')
    end select
  end function apply_operation
  
  pure elemental real(wp) function get_unit_value(unit_name, units) result(val)
    character(*), intent(in) :: unit_name
    type(unit_factor_t), intent(in) :: units
    character(*), parameter :: proc_name = 'get_unit_value'
  
    val = 0.0_wp
    select case(unit_name)
    case('meter', 'm')
      val = units%meter
    case('second', 's')
      val = units%second
    case('mole', 'mol')
      val = units%mole
    case('joule', 'Joule', 'J')
      val = units%joule
    case('kelvin', 'Kelvin', 'celsius', 'Celsius')
      val = units%kelvin
    case('minute', 'min')
      val = units%minute
    case('hour', 'h')
      val = units%hour
    case('nodes')
      ! FIXME: this is not true for dx /= meter
      ! ideally, meter should be equal to dx
      val = 1.0_wp
    case('kilogram', 'kg')
      val = units%kilogram
    case('gram', 'g')
      val = units%gram
    case('newton', 'Newton', 'N')
      val = units%newton
    case('pascal', 'Pascal', 'Pa')
      val = units%pascal
    case default
      call error_msg(proc_name, 'Unknown unit name: ' // unit_name, module_name_opt = module_name)
    end select
  end function get_unit_value
end module mod_json
