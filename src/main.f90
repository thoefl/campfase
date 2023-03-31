! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

! Features that were implemented previously but currently not functional:
! - Non-periodic boundary conditions: Based on Fourier continuation
! - Restoring from previous simulations using HDF files as input
! 
!TODO: limit amount of label calls for remapping
!TODO: find a method to detect ringing in remapping procedure to abort early
!TODO: strided calculation of gradients (i.e. dimensions > 1) could probably be faster
!TODO: find a method to implement variable surface energies while keeping semi-implicit algorithm
!TODO: reform groups to exclude inactive order parameters from MPI collectives (?)
!FIXME: adaptive time stepping, grain remapping, etc. should use pure step timers, as others do not make much sense
!TODO: Restructure HDF file output format, should ideally resemble json file structure

program campfase

  use hdf5, only: h5open_f, h5close_f
  use mod_fftw
  use mod_hdf5

  use mod_base_functions
  use mod_globals, only: verbosity_level
  use mod_parameters
  use mod_realtype
  use mod_regular_grid
  use mod_time_dep_vars
  use mod_spectral
  use mod_fft_osc
  use mod_fcdata
  use mod_phase_set
  use mod_order_parameter_set
  use mod_phase_field_energies
  use mod_json
  use mod_mu
  use mod_pfm
  use mod_time_step
  use mod_phase_properties
  use mod_species_properties
  use mod_phase_interface
  use mod_stopwatch
  use mod_phase_properties
  use mod_initial_field_settings
  use mod_backup
  use mod_nucleation
  use mod_user_input
  use mod_unit_factors
  use mod_constants
  use mod_nucleus
  use mod_random
  use mod_files
  !use liboctq, only: gtp_equilibrium_data, GSNOACS, zero, one, tqini, tqquiet, tqtgsw, tqrpfil, tqphsts, tqgpi, tqsetc, maxc, &
  !    enter_composition_set, tqgpi, ask_default_constitution, tqlc, tqlr
  use omp_lib, only: omp_get_max_threads

  implicit none

  character(*), parameter :: proc_name = 'main'

  !!!!!!!!!!!!!!!!!!!!!!!
  !! Runtime variables !!
  !!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Date and time
  ! start_date:        Date at the beginning of the simulation, format YYYYMMDD, 
  !                    used for output filename
  ! start_time:        Time at the beginning of the simulation, format 
  !                    hhmmss[ms][ms][ms][ms], used in part for output filename
  character(8) :: start_date
  character(10) :: start_time

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! JSON file used to store simulation settings
  ! simulation_name is set to the JSON filename without extension
  character(:), allocatable :: json_filename
  type(json_file) :: json

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Simulation settings as defined per json input file
  type(regular_grid) :: grid
  type(time_step) :: dt
  type(fft_grid_max_freq_oscillations), allocatable :: fft_osc_main
  type(hdf_file), target :: hdf
  type(nucleation_data) :: nucleation
  type(grain_reindex_settings) :: grain_reindex
  type(fourier_continuation_settings) :: fourier_cont
  type(restore_file_settings) :: restore
  type(initial_field_settings) :: init_field
  type(piecewise_linear) :: temperature
  real(wp) :: pressure
  type(species_properties_t), allocatable :: species_properties(:)
  integer :: n_species, n_phases
  type(phase_properties_t), allocatable, target :: phase_properties(:)
  type(phase_interface) :: phase_inter
  type(chemical_energy_model) :: chemical_energy
  type(check_settings) :: check
  type(benchmark_settings) :: bench_fft, bench_coarray
  type(opencalphad_settings) :: oc_settings

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Simulation step and various types of "stopwatches" used for output/simulation
  ! timing
  integer(int64) :: step = 0
  type(system_stopwatch) :: wall_clock
  type(manual_stopwatch) :: sim_clock
  type(cpu_stopwatch) :: cpu_clock
  type(timer_set) :: simulation_stop, console_output, image_output, oscillation_limit, oscillation_cooldown
  type(backup_t) :: backup

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! M_eta: Mobility value for eta field
  ! mu_average: Average value of diffusion potential for species_properties
  ! randseed: Random seed storage (per image)
  ! phase_amount_oc: Phase amounts as determined by OpenCalphad
  ! oc_total_fu: Total amount of formula units calculated from phase_amount_oc

  real(wp) :: M_eta = 0.0_wp
  real(wp) :: mu_average
  integer, allocatable :: randseed(:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! General purpose variables
  ! i, j: Counters for various loops
  ! rand: Random number variable
  ! randseed :: Random seed variable
  integer :: i, j, dim, node, species, phase, n_chars
  logical :: file_exists

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Status/procedure return variables
  integer :: hdf_err, fft_err

  !TODO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Under construction
  real(wp), dimension(:), allocatable :: surf_x_node, surf_x_node_full
  real(wp), dimension(:,:), allocatable :: c_init
  logical :: use_bdfab_backup = .false.
  logical :: concentration_is_conserved
  logical :: apply_oscillation_limit = .false.
  type(nucleus) :: new_nucleus
  integer :: total_n_nucl_candidates

  ! Universal constants
  type(unit_factor_t) :: units
  type(constants_t) :: constants
























  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Scalar grid variables

  integer(int64) :: step_last_print = 0
  real(wp) :: t_last_print = 0.0_wp
  real(wp) :: steps_per_second

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Field and vector grid variables
  ! c:           Concentration field
  ! mu_c:        Chemical potential field, mu_c = f(c)
  ! x:           X coordinates
  ! y:           Y coordinates
  type(order_parameter_set) :: order_parameters
  type(phase_set) :: phases
  type(fft_field), dimension(:), allocatable :: c_field
  type(fft_field) :: mu_c, mu_eta
  type(fft_field), allocatable :: bench_fft_field
  integer :: n_fft_threads
  real(wp), allocatable, dimension(:) :: bench_coarray_data
  logical, allocatable, dimension(:) :: eta_filled[:]
  real(wp), allocatable, dimension(:) :: eta_remap[:]
  type(phase_field_energies) :: f, df_deta, df_deta_numeric
  type(phase_field_energies), allocatable, dimension(:) :: df_dc, df_dc_numeric
  type(time_dependent_real) :: df_chem_deta_common
  integer :: n_eta, n_eta_sol, eta_idx, phase_idx, species_idx, order_idx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Nucleation
  ! p_nucl:           Local nucleation chance, depending on temperature and
  !                   concentration
  ! rand_nucl:        Matrix of random variables for comparison with p_nucl to
  !                   determine if and where nucleation occurs
  ! nucl_mask:        Nucleation mask for applying the nucleating chem. potential
  ! nucl_compmask:    "Target" nucleation mask used for comparison and removal of
  !                   nucleation nodes
  ! nucl_candidate_nodes:  List of nucleation candidates to check for eligibility,
  !                   i.e. nodes where p_nucl >= rand_nucl
  integer, allocatable :: nucl_candidate_nodes(:), nucl_candidate_phase_idx(:), nucl_candidate_scramble(:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! File input/output variables
  ! hdf:      Custom type holding file name, handle, dimensions
  ! output_...:       Custom type holding associated hdf_file, data type, size, 
  !                   name for different variables
  ! input_file:       File handle used for loading
  ! input_dset:       Dataset handle used for loading
  ! input_prop:       Property list used for loading
  ! input_dspace:     Dataspace handle used for loading
  ! input_memspace:   Memory space handle used for loading
  ! hdf_err:          HDF5 error checking variable
  ! hdf%filename:    Filename of the output file used
  ! hdf%realtype:    Real type used for output, single or double precision
  !                   depending on the choice of 'wp' in mod_types.f90
  ! input_dims:       Dimensions for loading from file
  ! input_maxdims:    Maximum extent of input file
  ! input_offset:     Offset to load from input file (typically last entry)
  ! hdf5_maxdims:     Maximum extent of data storage
  ! hdf5_startdims:   Starting extent of data storage
  ! hdf5_arraychunk:  Chunk size used for array-type (matrix, vector) outputs
  ! hdf5_scalarchunk: Chunk size used for scalar outputs
  ! hdf5_xchunk:      Chunk size used for x-dependent outputs
  ! hdf5_ychunk:      Chunk size used for y-dependent outputs
  ! t_nextsave:       Next simulation time to output to file, given that
  !                   save_by_time is enabled and save interval save_dt > 0

  type(hdf_set) :: output_t_now, output_therm_now, output_t_cpu, output_t_wall, output_dt, output_step, output_kappa, &
    output_f_total_sum, output_nucl_size, &
    output_n_t_error, output_randseed, output_dt_errtol, output_nucl_dt, output_m_order, &
    output_if_nodes, output_if_energy, output_sqsum, output_eta_idx, &
    output_use_bdfab, output_stab_coeff_bdfab, output_stab_coeff_be, &
    output_surf_x_node, output_V_m, output_D_factor, output_G_m_factor, &
    output_eta, output_eta_error, output_n_eta_limit, output_eta_error_field, output_mu_eta, &
    output_c, output_c_error, output_c_error_field, output_n_c_limit, output_mu_c, output_D, &
    output_mu_eta_common, output_M_eta, output_phi
  type(hdf_set), dimension(:), allocatable :: output_n_main_osc_limit, output_main_osc_error, output_main_maxfreq
  type(hdf_set), dimension(:), allocatable :: output_axes_n_nodes, output_axes_spacing
  type(hdf_set), dimension(:), allocatable :: output_p_coeff, output_k_coeff, output_c_equi, output_omega
  type(hdf_set) :: output_f_total, output_f_grad, output_backup_step
  !integer(hid_t) :: input_file, input_dset, input_prop, input_dspace, input_memspace
  logical :: halt_program = .false.
  logical :: continue_execution = .true.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Adaptive time stepping
  ! dt%adaptive%global_err_max:       Time discretization error, calculated by Richardson
  !                 extrapolation
  ! n_saved:        Number of times saved to output file
  ! restore_backup: Restore a previously stored backup
  ! dt_reclac:      Controls flow for stages of dt error calculation
  ! c_half:         Concentration after two steps with half dt for error calc.
  real(wp) :: step_factor = 0.0_wp
  real(wp) :: total_step_factor = 0.0_wp
  integer :: n_saved = 0
  real(wp), dimension(:), allocatable :: main_var_full
  real(wp) :: t_full
  integer :: n_var_dt_limit = 0
  integer :: last_osc_step = -1E8

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Fourier continuation
  ! cont_mat_l: Left side extrapolation matrix
  ! cont_mat_r: Right side extrapolation matrix
  ! cont_bc_l:  Left boundary condition vector
  ! cont_bc_r:  Right boundary condition vector
  ! cont_bc:  Sum of left and right boundaries for concentration field
  ! type(fourier_cont_data) :: fc_data
  ! integer :: fc_dl, fc_dr

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Extended output variables
  !
  real(wp), allocatable, dimension(:) :: var_error_field

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Coarray specific
  type(team_type) :: variable_team, phase_team, core_team
  integer :: phase_team_main_image

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! OC specific
  type(gtp_equilibrium_data), pointer :: current_equi
  character(len=24), dimension(maxc) :: component_names
  real(wp), dimension(:), allocatable :: c_target
  real(wp), dimension(:), allocatable :: x_mean, c_mean, x_mean_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FFT benchmarking
  real(wp) :: t_forward_2d = 0.0_wp
  real(wp) :: t_backward_2d = 0.0_wp
  real(wp), allocatable, dimension(:) :: t_grad
  real(wp) :: t_laplacian = 0.0_wp

  real(wp), allocatable, dimension(:) :: t_co_sum, t_co_broadcast, t_co_max, t_co_min

  !!!!!!!!!!!!!!!!!!!!!!!!
  !! Executable section !!
  !!!!!!!!!!!!!!!!!!!!!!!!

  ! Get starting wall and cpu time
  if(this_image() == 1) call date_and_time(date=start_date, time=start_time)
  call co_broadcast(start_time, 1)
  call co_broadcast(start_date, 1)
  call wall_clock%init()
  call sim_clock%init()
  call cpu_clock%init()
  call sim_clock%start()

  ! Make results deterministic by fixing random seed
  ! init_random_seed is processor-dependent, use a custom portable routine for reproducibility
  !call init_random_seed(repeatable = .true., image_distinct = .false.)
  call init_random_seed(repeatable = .true., image_distinct = .false.)

  if(this_image() == 1) then
    ! Get the JSON input filename as command line argument or from user input
    ! Check if a filename was specified via command line argument
    select case(command_argument_count())
    case(1)
      call get_allocatable_command_argument(1, json_filename)
    case(2:)
      call warning_msg(proc_name, 'Command line arguments after first argument discarded')
      call get_allocatable_command_argument(1, json_filename)
    end select
    ! If yes, check if the file can be loaded
    file_exists = .false.
    if(allocated(json_filename)) then
      inquire(file=trim(json_filename), exist=file_exists)
    end if
    ! If not, ask the user to provide a filename
    if(.not. file_exists) then
      if(allocated(json_filename)) then
        call warning_msg(proc_name, 'Invalid JSON filename declared in command line')
        deallocate(json_filename)
      end if
      write(output_unit, '(a)') 'Please provide the filename containing the simulation parameters'
      do
        allocate(character(128) :: json_filename)
        read (input_unit, '(a)') json_filename
        if(trim(json_filename) == 'exit') then
          halt_program = .true.
          exit
        end if
        inquire(file=trim(json_filename), exist=file_exists)
        if(.not. file_exists) then
          if(.not. string_contains_character(json_filename, '.')) then
            json_filename = trim(json_filename) // json_extension
            inquire(file=trim(json_filename), exist=file_exists)
          end if
        end if
        if(file_exists) then
          exit
        else
          write(output_unit, '(a)') 'File does not exist, try again or type "exit"'
          deallocate(json_filename)
        end if
      end do
    end if
    n_chars = len(json_filename)
  end if
  call co_any(halt_program)
  if(halt_program) stop
  call co_broadcast(n_chars, 1)
  if(this_image() /= 1) allocate(character(n_chars) :: json_filename)
  call co_broadcast(json_filename, 1)

  call start_critical_in_order()
  ! Load JSON parameter input file
  ! Initialize
  call json%initialize(stop_on_error = .true.)
  call json%load(filename = json_filename)

  ! Read data
  call json%get('output.console.verbosity', verbosity_level)
  call json_calculate_units(json, units)
  constants = constants_t(units)
  call json_read_grid(json, constants, units,  grid)
  call json_read_timer_set_data(json, constants, units, sim_clock, wall_clock, 'simulation_stop', simulation_stop)
  call json_read_timer_set_data(json, constants, units, sim_clock, wall_clock, 'output.console', console_output)
  call json_read_timer_set_data(json, constants, units, sim_clock, wall_clock, 'output.image', image_output)
  call json_read_timer_set_data(json, constants, units, sim_clock, wall_clock, &
    'fft_oscillation_check.adaptive_dt_limit.reduced_max_factor_duration', oscillation_limit, late_allowed_opt = .true.)
  call json_read_timer_set_data(json, constants, units, sim_clock, wall_clock, &
    'fft_oscillation_check.adaptive_dt_limit.cooldown', oscillation_cooldown, late_allowed_opt = .true.)
  call json_read_hdf_output_data(json, constants, units, grid, start_date, start_time, sim_clock, wall_clock, hdf)
  call json_read_time_step(json, constants, units, sim_clock, wall_clock, dt)
  call json_read_nucleation_data(json, constants, units, sim_clock, wall_clock, nucleation)
  call json_read_grain_reindex_data(json, constants, units, sim_clock, wall_clock, grain_reindex)
  call json_read_fourier_continuation_settings(json, fourier_cont)
  call json_read_restore_file_settings(json, restore)
  call json_read_opencalphad_settings(json, oc_settings)
  call json_read_species(json, constants, units, species_properties, n_species)
  call json_read_phase_properties(json, constants, units, phase_properties, n_phases)
  call json_read_pressure_temperature_data(json, constants, units, temperature, pressure)
  call json_read_initial_conditions(json, constants, units, grid, species_properties, phase_properties, init_field)
  call json_read_interfacial_data(json, constants, units, phase_inter)
  call json_read_chemical_energy_data(json, constants, units, chemical_energy)
  call json_read_check_settings(json, check)
  call json_read_benchmark_settings(json, 'fft', bench_fft)
  call json_read_benchmark_settings(json, 'coarray', bench_coarray)

  if(n_phases + n_species > num_images()) then
    call error_msg(proc_name, 'Number of species + phases cannot be larger than number of MPI/coarray images')
  end if

  ! Close file
  call json%destroy()
  call end_critical_in_order()
  sync all

  if(this_image() == 1) then
    if(.not. simulation_stop%is_enabled()) then
      write(output_unit, '(a)', advance='no') 'This simulation will run indefinitely. '
    else
      write(output_unit, '(a)') 'This simulation will stop if any of the following conditions is fulfilled:'
      if(allocated(simulation_stop%wall_timer)) then
        write(output_unit, list_value_format, advance='no') 'wall time', simulation_stop%wall_timer%get_dt_trigger()
        write(output_unit, '(a)') ' s'
      end if
      if(allocated(simulation_stop%sim_timer)) then
        write(output_unit, list_value_format, advance='no') &
          'simulation time', simulation_stop%sim_timer%get_dt_trigger()*units%second
        write(output_unit, '(a)') ' s'
      end if
      if(allocated(simulation_stop%step_timer)) then
        print list_value_format, 'simulation step', simulation_stop%step_timer%get_dstep_trigger()
      end if
    end if
    call binary_user_choice('y', 'n', .true., continue_execution, prompt_opt = 'Is this correct?')
  end if
  call co_all(continue_execution)
  if(.not. continue_execution) stop

  fft_err = fftw_init_threads()
  n_fft_threads = omp_get_max_threads()
  if(fft_err == 0) call error_msg('main', 'Could not initialize FFTW!')
  if(verbosity_level >= 3) then
    call print_message('Initializing FFTW with ' // convert_to_char(n_fft_threads) // ' threads')
  end if
  call fftw_plan_with_nthreads(n_fft_threads)

  n_eta = num_images() - n_species
  n_eta_sol = n_eta - (n_phases - 1)

  call fft_load_wisdom(grid)
  call phases%init(grid, n_species, n_phases, species_properties%name, phase_properties, chemical_energy)
  call f%init('f', '', grid)
  ! Allocate arrays
  allocate(eta_filled(grid%n_nodes_total)[*])
  allocate(main_var_full(grid%n_nodes_total))
  allocate(c_target(n_phases))
  allocate(x_mean(n_species))
  allocate(c_mean(n_species))
  allocate(x_mean_init(n_species))
  allocate(c_field(n_species))

  ! Allocate extended output variables if selected
  if(hdf%extended) then
    allocate(var_error_field(grid%n_nodes_total))
    var_error_field = 0.0_wp
  else
    allocate(var_error_field(0))
  end if
  if(check%debug_level >= 5) then
    allocate(df_dc_numeric(n_species))
    do i = 1, n_species
      call df_dc_numeric(i)%init('df', '_dc_' // species_properties(i)%name // '_numeric', grid)
    end do
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  randseed = get_random_seed()

  if(bench_coarray%enabled) then
    allocate(t_co_broadcast(bench_coarray%count))
    allocate(t_co_sum(bench_coarray%count))
    allocate(t_co_max(bench_coarray%count))
    allocate(t_co_min(bench_coarray%count))
    t_co_broadcast = 0
    t_co_sum = 0
    t_co_max = 0
    t_co_min = 0
    call init_random_seed(repeatable = .true., image_distinct = .true.)
    allocate(bench_coarray_data(grid%n_nodes_total))
    call print_message('Starting coarray benchmark...', co_img_opt = 1)
    do i = 1, bench_coarray%count
      call random_number(bench_coarray_data)
      sync all
      call wall_clock%tic(t_co_broadcast(i))
      call co_broadcast(bench_coarray_data, mod_idx(i, num_images()))
      call wall_clock%toc(t_co_broadcast(i))
      call random_number(bench_coarray_data)
      sync all
      call wall_clock%tic(t_co_sum(i))
      call co_sum(bench_coarray_data)
      call wall_clock%toc(t_co_sum(i))
      call random_number(bench_coarray_data)
      sync all
      call wall_clock%tic(t_co_max(i))
      call co_max(bench_coarray_data)
      call wall_clock%toc(t_co_max(i))
      call random_number(bench_coarray_data)
      sync all
      call wall_clock%tic(t_co_min(i))
      call co_min(bench_coarray_data)
      call wall_clock%toc(t_co_min(i))
    end do
    if(this_image() == 1) then
      write (output_unit, '(a)', advance = 'yes') 'Coarray benchmark results: '
      write (output_unit, '(a)', advance = 'yes')&
        'wall time [s]:           co_sum     co_broadcast           co_max           co_min'
    end if
    call start_critical_in_order()
    write (output_unit, '(a)', advance = 'no') '                '
    write (output_unit, '(2(es7.2e1, " "), " ")', advance = 'no') mean(t_co_sum), maxval(t_co_sum)
    write (output_unit, '(2(es7.2e1, " "), " ")', advance = 'no') mean(t_co_broadcast), maxval(t_co_broadcast)
    write (output_unit, '(2(es7.2e1, " "), " ")', advance = 'no') mean(t_co_max), maxval(t_co_max)
    write (output_unit, '(2(es7.2e1, " "), " ")', advance = 'yes') mean(t_co_min), maxval(t_co_min)
    call end_critical_in_order()
    call millisleep(50)
    sync all
    deallocate(bench_coarray_data)
  end if

  if(bench_fft%enabled) then
    call wall_clock%start()
    ! Make sure everything is allocated and ready
    allocate(bench_fft_field)
    allocate(t_grad(grid%n_dims))
    t_grad = 0
    call bench_fft_field%init(grid, 'fft_benchmark', .false.)
    call bench_fft_field%fft_allocate_check()
    call random_number(bench_fft_field%real%vals)
    bench_fft_field%real%time = 0
    call bench_fft_field%calc_grad()
    call print_message('Starting FFT benchmark...', co_img_opt = 1)
    ! Start the benchmark
    do i = 1, bench_fft%count
      call bench_fft_field%reset
      call random_number(bench_fft_field%real%vals)
      bench_fft_field%real%time = 0
      call wall_clock%tic(t_forward_2d)
      call bench_fft_field%forward_alldim()
      call wall_clock%toc(t_forward_2d)
      call bench_fft_field%real%invalidate()
      call wall_clock%tic(t_backward_2d)
      call bench_fft_field%backward_alldim()
      call wall_clock%toc(t_backward_2d)
      do dim = 1, grid%n_dims
        call wall_clock%tic(t_grad(dim))
        call bench_fft_field%calc_grad(dim)
        call wall_clock%toc(t_grad(dim))
      end do
      call bench_fft_field%fft%invalidate()
      call wall_clock%tic(t_laplacian)
      call bench_fft_field%calc_real_laplacian()
      call wall_clock%toc(t_laplacian)
    end do
    sync all
    if(this_image() == 1) then
      write (output_unit, '(a)', advance = 'yes') 'FFT benchmark results: '
      write (output_unit, '(a)', advance = 'no') 'wall time [s]:    forward all dim   backward all dim'
      do dim = 1, grid%n_dims
        write (output_unit, '(a)', advance = 'no') '    real gradient '
        write (output_unit, '(i1)', advance = 'no') dim
      end do
      write (output_unit, '(a)', advance = 'yes') '     real laplacian'
    end if
    call start_critical_in_order()
    write (output_unit, '(a)', advance = 'no') '                 '
    write (output_unit, '(es16.9)', advance = 'no') t_forward_2d
    write (output_unit, '(a)', advance = 'no') '   '
    write (output_unit, '(es16.9)', advance = 'no') t_backward_2d
    write (output_unit, '(a)', advance = 'no') '   '
    do dim = 1, grid%n_dims
      write (output_unit, '(es16.9)', advance = 'no') t_grad(dim)
      write (output_unit, '(a)', advance = 'no') '   '
    end do
    write (output_unit, '(es16.9)', advance = 'yes') t_laplacian
    call millisleep(10)
    call end_critical_in_order()
    sync all
    deallocate(bench_fft_field)
    call wall_clock%stop()
    if(this_image() == 1) then
      call binary_user_choice('y', 'n', .true., continue_execution, prompt_opt = 'Continue?')
    end if
    call co_all(continue_execution)
    if(.not. continue_execution) stop
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  x_mean = species_properties%x_init
  call opencalphad_setup(oc_settings, species_properties, phase_properties, pressure, units, current_equi, component_names)

  if(verbosity_level >= 5) call print_message('Forming core teams...')
  if(this_image() <= n_species + n_phases .or. check%debug_level >= 5) then
    form team(primary_team_num, core_team)
    do i = 1, n_species
      call c_field(i)%init(grid, c_name // '_' // species_properties(i)%name, dt%bdfab_enabled)
    end do
  else
    form team(secondary_team_num, core_team)
  end if

  phase_idx = max(min(this_image() - n_species, n_phases), 0)

  ! FIXME: Setting phase indices could be prettier...
  do i = 1, n_phases
    phase_properties(i)%idx%phase = i
    phase_properties(i)%idx%image = n_species + i
    phase_properties(i)%idx%eta%first = i
    if(i == n_phases) then
      phase_properties(i)%n_order = 1 + num_images() - n_species - n_phases
    else
      phase_properties(i)%n_order = 1
    end if
    phase_properties(i)%idx%eta%last = phase_properties(i)%idx%eta%first + phase_properties(i)%n_order - 1
  end do

  if(phase_idx > 0) then
    species_idx = 0
    eta_idx = this_image() - n_species
    if(verbosity_level >= 5) call print_message('Forming phase_team ' // convert_to_char(phase_idx))
    form team(phase_idx, phase_team)
    if(verbosity_level >= 5) call print_message('Forming variable_team ' // convert_to_char(eta_team_num))
    form team(eta_team_num, variable_team)
  else
    eta_idx = 0
    order_idx = 0
    species_idx = this_image()
    if(verbosity_level >= 5) call print_message('Forming phase_team ' // convert_to_char(c_team_num))
    form team(c_team_num, phase_team)
    if(verbosity_level >= 5) call print_message('Forming variable_team ' // convert_to_char(species_idx))
    form team(species_idx, variable_team)
    call mu_c%init(grid, mu_name // '_' // c_name // '_' // trim(species_properties(species_idx)%name), dt%bdfab_enabled)
    ! Setup matrices for Fourier continuation
    if(.not. all(c_field(species_idx)%real%grid%axes%bound_cond%is_periodic)) then
      call error_msg(proc_name, 'Non-periodic boundary conditions currently not implemented')
      !if(c_field(species_idx)%bc%type_max == bc_dirichlet) then
      !    fc_dr = fc_d_dirichlet
      !else if(c_field(species_idx)%bc%type_max == bc_neumann) then
      !    fc_dr = fc_d_reflect
      !else
      !    stop 'ERROR: Boundary condition type not recognized for south'&
      !         'Fourier continuation.'
      !end if
      !if(c_field(species_idx)%bc%type_min == bc_dirichlet) then
      !    fc_dl = fc_d_dirichlet
      !else if(c_field(species_idx)%bc%type_min == bc_neumann) then
      !    fc_dl = fc_d_reflect
      !else
      !    stop 'ERROR: Boundary condition type not recognized for north'&
      !         'Fourier continuation.'
      !end if
      !call set_fcdata(c_field(species_idx)%bc, fc_nodes, fc_dl, fc_dr, fc_w, fc_data)
    end if
  end if
  if(verbosity_level >= 5) call print_message('phase_idx = '&
    // convert_to_char(phase_idx) // ' species_properties idx = ' // convert_to_char(species_idx))

  change team(phase_team)
  if(phase_idx > 0) then
    order_idx = this_image()
    if(order_idx == 1) then
      phase_team_main_image = this_image()
    else
      phase_team_main_image = -1
    end if
    call co_max(phase_team_main_image)
    if(phase_team_main_image <= 0 .or. phase_team_main_image > num_images()) then
      error stop 'ERROR: Could not determine phase team main image for broadcasting!'
    end if
  end if
  end team

  ! TODO: Best way to calculate initial equilibrium?
  ! call opencalphad_equilib_phase(273.15_wp, pressure, x_mean, phase_idx, units, constants, &
  !    current_equi, component_names, c_field, phases, gridmin_opt = .true.)
  ! Disable automatic generation of composition sets
  call tqtgsw(GSNOACS)
  call opencalphad_equilib_phase(temperature%get_value_at(sim_clock%get_time()), pressure, x_mean, phase_idx, units, constants, &
    current_equi, component_names, c_field, phases, gridmin_opt = .true.)

  if(check%equilibrium) then
    if(this_image() == 1) then
      call opencalphad_equilib_phase(temperature%get_value_at(sim_clock%get_time()), pressure, x_mean, phase_idx, &
        units, constants, current_equi, component_names, c_field, phases)
      ! FIXME: allow setting number of checks
      do step = 1, 10000
        call opencalphad_equilib_phase(temperature%get_value_at(sim_clock%get_time()), pressure, x_mean, phase_idx, &
          units, constants, current_equi, component_names, c_field, phases)
        call tqlc(output_unit, current_equi)
        call tqlr(output_unit, current_equi)
      end do
    end if
    sync all
    stop
  end if

  !FIXME: eta_labels does not allow partial presence of phases
  change team(variable_team)
  if(team_number() == eta_team_num) then
    call order_parameters%init(grid, phase_properties, sim_clock%get_time())
    call order_parameters%set_if_param(phase_inter)
    call mu_eta%init(grid, &
      mu_eta_name // '_' // trim(phase_properties(phase_idx)%name) // '_' // convert_to_char(order_idx), &
      dt%bdfab_enabled)
    call df_chem_deta_common%init(grid, mu_eta_name // '_' // trim(phase_properties(phase_idx)%name) // '_' // common_name)
    call order_parameters%apply_initial_conditions(init_field)
    call nucleation%init(grid, order_parameters%eta)
    call df_deta%init('df', '_d' // order_parameters%eta%real%var_name, grid)
    if(check%debug_level >= 5) call df_deta_numeric%init('df', '_d' // order_parameters%eta%real%var_name // '_numeric', grid)
  else
    ! Since the concentration fields are idle otherwise, precalculate the FFT plans here
    call c_field(species_idx)%fft_allocate_check()
  end if
  end team

  ! FIXME: This is a rather crude workaround to define n_eta on non-eta images
  call co_broadcast(order_parameters%n_eta, n_species + 1)

  allocate(df_dc(n_species))
  do i = 1, n_species
    call df_dc(i)%init('df', '_dc_' // trim(species_properties(i)%name), grid)
  end do

  if(dt%adaptive%enabled) then
    call start_critical_in_order()
    call json%initialize(stop_on_error = .true.)
    call json%load(filename = json_filename)
    if(phase_idx > 0) then
      call json_read_fft_oscillation_data(json, order_parameters%eta%real%grid, order_parameters%eta%real%var_name, fft_osc_main)
    else
      call json_read_fft_oscillation_data(json, c_field(species_idx)%real%grid, c_field(species_idx)%real%var_name, fft_osc_main)
    end if
    call json%destroy()
    call end_critical_in_order()
    sync all

    dt%adaptive%max_factor_osc = dt%adaptive%max_factor**(2*dt%adaptive%timer%step_timer%get_dstep_trigger() &
      /real(oscillation_limit%step_timer%get_dstep_trigger(), wp))
  end if

  change team(variable_team)
  if(team_number() == eta_team_num) then
    call order_parameters%calc_phi(phases)
  end if
  end team
  call phases%phi_broadcast()
  ! TODO: should generalize surf_x_node to all dimensions that are non-periodic
  allocate(surf_x_node(grid%axes(1)%n_nodes))
  allocate(surf_x_node_full(grid%axes(1)%n_nodes))
  surf_x_node = 0.0_wp
  if(.not. all(grid%axes%bound_cond%is_periodic)) then
    call find_surface(grid, phases%phi(1)%vals, surf_x_node)
  end if

  allocate(c_init(n_phases, n_species))
  change team(core_team)
  if(team_number() == primary_team_num) then
    if(init_field%n_precipitates > 0) then
      c_init = get_c_init(phases)
    else
      call tqgpi(phase, init_field%matrix_phase%name, current_equi)
      call tqgnp(i, j, current_equi)
      call tqphsts(-1, -3, j, zero, current_equi)
      call tqphsts(phase, 0, j, one, current_equi)
      call opencalphad_equilib_phase(temperature%get_value_at(sim_clock%get_time()), &
        pressure, x_mean, phase_idx, units, constants, current_equi, component_names, c_field, phases)
      c_init = get_c_init(phases, &
        species_properties%x_init*init_field%matrix_phase%mol_per_formula_unit/init_field%matrix_phase%V_m)
      do i = 1, size(phase_properties)
        call tqgpi(phase, phase_properties(i)%name, current_equi)
        if(verbosity_level >= 5) call print_message('Opencalphad: Setting phase '&
          // phase_properties(i)%name // ' at phase tuple idx ' // convert_to_char(phase) // ' to status entered.')
        call tqphsts(phase, 0, j, one, current_equi)
      end do
      call opencalphad_equilib_phase(temperature%get_value_at(sim_clock%get_time()), &
        pressure, x_mean, phase_idx, units, constants, current_equi, component_names, c_field, phases)
    end if
    do concurrent(node = 1:grid%n_nodes_total, species = 1:n_species)
      c_field(species)%real%vals(node) = phases%phi_mix(c_init(:,species), node)
    end do
    c_field%real%time = phases%phi(1)%time
  end if
  end team

  !do i = 1,n_species
  !    if(.not. c_field(i)%bc%is_periodic) then
  !        call fourier_clean_var_lb(c_field(i)%real%vals, surf_x_node)
  !        call prediffuse_interface(c_field(i), surf_x_node, 6.0_wp, 50)
  !    end if
  !end do

  ! FIXME: This should go in a subroutine, preferrably type-bound to hdf
  ! Setup HDF5 if saving or loading was requested
  if(hdf%output%is_enabled() .or. restore%enabled) then
    ! FIXME
    if(restore%enabled) call error_msg(proc_name, 'Restoring from HDF file currently not implemented')
    call h5open_f(hdf_err)
  end if

  !if(species_idx > 0) then
  !    if(nucl_init_random) then
  !        ! Select random nodes for the initial seeds
  !        do while (nucl_sites < nucl_initseeds)
  !            i = random_int(1, grid%n_nodes_total)
  !            !if(eligible_nucleus(c_c, c_ni, c_matrix, i, nucleation%r_nodes/dx, phase_inter%get_width(), nucl_node, &
  !            !    field_bc, sim_mode, fc_nodes)) then
  !            !    call nucl_add(i, nucl_node, nucl_sites, nucl_counter)
  !            !    if(verbosity_level >= 1) call print_message('Nucleating...', step)
  !            !end if
  !        end do
  !    else
  !        do i = 1, nucl_initseeds
  !            call nucl_add((nx + 1)*ny/2, nucl_node, nucl_sites, nucl_counter)
  !                if(verbosity_level >= 1) call print_message('Nucleating...', step)
  !        end do
  !    end if
  !    call nucl_setmask(&
  !        nucl_node, nucleation%r_nodes/dx, nucl_mask, phase_inter%get_width(), field_bc)
  !    c_prev_nucl_node = c_c(nucl_node)
  !end if
  !call co_broadcast(nucl_mask, c_c_idx)
  !call co_broadcast(nucl_sites, c_c_idx)


  if(verbosity_level >= 5) call print_message('Saving initial data to hdf5...', step)

  call cpu_clock%start()
  call wall_clock%start()

  ! Setup HDF5 file for saving data if requested
  if(hdf%output%is_enabled()) then
    if(this_image() == 1) then
      call hdf%touch()
      call copy_file(json_filename, hdf%filename // '.json')
    end if

    call start_critical_in_order()
    call hdf%open()
    if(this_image() == 1) then
      allocate(output_axes_n_nodes(grid%n_dims))
      allocate(output_axes_spacing(grid%n_dims))
      ! FIXME: grid parameter output should be prettier
      do dim = 1, grid%n_dims
        call output_axes_n_nodes(dim)%init(hdf, grid%axes(dim)%n_nodes, 'n_nodes_' // grid%axes(dim)%label)
        call output_axes_spacing(dim)%init(hdf, grid%axes(dim)%spacing*units%meter, 'node_spacing_' // grid%axes(dim)%label)
      end do
      call output_t_now%init(hdf, sim_clock%get_time(), t_now_name)
      call output_therm_now%init(hdf, temperature%get_value_at(sim_clock%get_time()), temp_now_name)
      call output_t_cpu%init(hdf, cpu_clock%get_time(), t_cpu_name)
      call output_t_wall%init(hdf, wall_clock%get_time(), t_wall_name)
      call output_dt%init(hdf, dt%val*units%t, dt_name)
      call output_V_m%init(hdf, units%V_m, V_m_name)
      call output_D_factor%init(hdf, units%D, D_factor_name)
      call output_G_m_factor%init(hdf, units%G_m, G_m_factor_name)
      call output_step%init(hdf, step, step_name)
      call output_backup_step%init(hdf, step, backup_step_name)
      call output_n_t_error%init(hdf, dt%adaptive%n_min_factor, n_t_error_name)
      !call output_nucl_k1%init(hdf, nucleation%k_1, nucl_k1_name)
      !call output_nucl_k2%init(hdf, nucleation%k_2, nucl_k2_name)
      call output_dt_errtol%init(hdf, dt%adaptive%rel_err_tol, dt_errtol_name)
      call output_nucl_size%init(hdf, nucleation%r_nodes, nucl_size_name)
      !call output_nucl_node%init(hdf, nucl_node, nucl_node_name)
      call output_randseed%init(hdf, [size(randseed)], randseed, randseed_name)
      call output_nucl_dt%init(hdf, sim_clock%get_time(), nucl_dt_name)
      call output_f_total_sum%init(hdf, f%total_sum, f_total_sum_name)
      call output_stab_coeff_bdfab%init(hdf, stab_coeff_bdfab, stab_coeff_bdfab_name)
      call output_stab_coeff_be%init(hdf, stab_coeff_be, stab_coeff_be_name)
      call output_use_bdfab%init(hdf, dt%use_bdfab, use_bdfab_name)
      call output_surf_x_node%init(hdf, [grid%axes(1)%n_nodes], surf_x_node, surf_x_node_name)
      call output_f_total%init(hdf, grid%axes%n_nodes, f%total%vals, f_total_name)
      call output_f_grad%init(hdf, grid%axes%n_nodes, f%gradient%vals, f_grad_name)

      ! Save constant values
      do dim = 1, grid%n_dims
        call output_axes_n_nodes(dim)%append(grid%axes(dim)%n_nodes)
        call output_axes_spacing(dim)%append(grid%axes(dim)%spacing*units%meter)
      end do
      ! FIXME: V_m factor should not be called 'V_m'
      call output_V_m%append(units%V_m)
      call output_D_factor%append(units%D)
      call output_G_m_factor%append(units%G_m)
      !call output_nucl_k1%append(nucleation%k_1)
      !call output_nucl_k2%append(nucleation%k_2)
      call output_nucl_size%append(nucleation%r_nodes)
      call output_dt_errtol%append(dt%adaptive%rel_err_tol)
      call output_stab_coeff_bdfab%append(stab_coeff_bdfab)
      call output_stab_coeff_be%append(stab_coeff_be)
    end if
    if(allocated(fft_osc_main)) then
      allocate(output_n_main_osc_limit(fft_osc_main%n_freqs))
      allocate(output_main_osc_error(fft_osc_main%n_freqs))
      allocate(output_main_maxfreq(fft_osc_main%n_freqs))
      do i = 1, fft_osc_main%n_freqs
        call output_n_main_osc_limit(i)%init(hdf, fft_osc_main%max_freqs(i)%n_limit, &
          n_osc_limit_name // '_' // fft_osc_main%var_name // '_' // fft_osc_main%max_freqs(i)%name)
        call output_main_osc_error(i)%init(hdf, [fft_osc_main%max_freqs(i)%error_nsave], &
          fft_osc_main%max_freqs(i)%get_error(), &
          osc_error_name // '_' // fft_osc_main%var_name // '_' // fft_osc_main%max_freqs(i)%name)
        call output_main_maxfreq(i)%init(hdf, fft_osc_main%max_freqs(i)%prev_maxfreq, &
          osc_maxfreq_name // '_' // fft_osc_main%var_name // '_' // fft_osc_main%max_freqs(i)%name)
      end do
    end if
    if(species_idx > 0) then
      call output_c%init(hdf,&
        grid%axes%n_nodes, c_field(species_idx)%real%vals*units%c,&
        trim(c_field(species_idx)%real%var_name))
      call output_D%init(hdf, &
        species_properties(species_idx)%diffusion_coefficient%get_value(temperature%get_value_at(sim_clock%get_time())), &
        D_name // '_' // trim(species_properties(species_idx)%name))
      if(dt%adaptive%enabled) then
        call output_c_error%init(hdf, dt%adaptive%local_err, &
          trim(c_field(species_idx)%real%var_name) // '_' // error_name)
        call output_n_c_limit%init(hdf, n_var_dt_limit, &
          trim(c_field(species_idx)%real%var_name) // '_' // n_limit_name)
        if(hdf%extended) then
          call output_c_error_field%init(hdf, grid%axes%n_nodes, var_error_field*units%c, &
            trim(c_field(species_idx)%real%var_name) // '_' // error_field_name)
        end if
      end if
      if(hdf%extended) then
        call output_mu_c%init(hdf, grid%axes%n_nodes, c_field(species_idx)%real%vals, &
          mu_name // '_' // trim(c_field(species_idx)%real%var_name))
      end if
    end if
    if(eta_idx > 0) then
      if(eta_idx == 1) then
        call output_if_nodes%init(hdf, phase_inter%get_width(), if_nodes_name)
        call output_if_nodes%append(phase_inter%get_width())
        call output_if_energy%init(hdf, phase_inter%get_energy(), if_energy_name)
        call output_kappa%init(hdf, phase_inter%get_kappa(), kappa_name)
        call output_m_order%init(hdf, phase_inter%get_m_order(), m_order_name)
        call output_sqsum%init(hdf, grid%axes%n_nodes, order_parameters%sqsum%vals, sqsum_name)
        call output_eta_idx%init(hdf, grid%axes%n_nodes, order_parameters%eta_labels%vals, eta_idx_name)
      end if
      call output_eta%init(hdf, grid%axes%n_nodes, order_parameters%eta%real%vals, &
        trim(order_parameters%eta%real%var_name))
      if(dt%adaptive%enabled) then
        call output_eta_error%init(hdf, dt%adaptive%local_err, &
          trim(order_parameters%eta%real%var_name) // '_' // error_name)
        call output_n_eta_limit%init(hdf, n_var_dt_limit, trim(order_parameters%eta%real%var_name) // '_' // n_limit_name)
        if(hdf%extended) then
          call output_eta_error_field%init(hdf, grid%axes%n_nodes, var_error_field, &
            trim(order_parameters%eta%real%var_name) // '_' // error_field_name)
        end if
      end if
      if(hdf%extended) then
        call output_mu_eta%init(hdf, grid%axes%n_nodes, order_parameters%eta%real%vals, &
          mu_name // '_' // trim(order_parameters%eta%real%var_name))
      end if
      if(order_idx == 1) then
        call output_M_eta%init(hdf, M_eta, M_eta_name // '_' // trim(phase_properties(phase_idx)%name))
        if(phase_idx < n_phases) then
          call output_phi%init(hdf, grid%axes%n_nodes, phases%phi(phase_idx)%vals, phi_name // '_' &
            // trim(phase_properties(phase_idx)%name))
        end if
        allocate(output_p_coeff(n_species))
        allocate(output_k_coeff(n_species))
        allocate(output_c_equi(n_species))
        allocate(output_omega(n_species))
        do i = 1, n_species
          call output_c_equi(i)%init(hdf, phases%phi(phase_idx)%equi_val(i)*units%c, &
            c_equi_name // '_' // trim(species_properties(i)%name) // '_' // trim(phase_properties(phase_idx)%name))
          call output_p_coeff(i)%init(hdf, phases%phi(phase_idx)%p_coeff(i), &
            p_coeff_name // '_' // trim(species_properties(i)%name) // '_' // trim(phase_properties(phase_idx)%name))
          call output_k_coeff(i)%init(hdf, phases%phi(phase_idx)%k_coeff(i), &
            k_coeff_name // '_' // trim(species_properties(i)%name) // '_' // trim(phase_properties(phase_idx)%name))
          call output_omega(i)%init(hdf, phases%phi(phase_idx)%omega(i), &
            omega_name // '_' // trim(species_properties(i)%name) // '_' // trim(phase_properties(phase_idx)%name))
        end do
        if(hdf%extended) then
          call output_mu_eta_common%init(hdf, grid%axes%n_nodes, df_chem_deta_common%vals, &
            df_chem_deta_common%var_name)
        end if
      end if
    end if
    call hdf%close()
    call end_critical_in_order()
  end if

  ! Restore from a previous file if requested
  !if(restore_file) then
  !    call h5fopen_f(hdf5_infilename, H5F_ACC_RDONLY_F, input_file, hdf_err)
  !    
  !    call h5dopen_f(input_file, c_c_name, input_dset, hdf_err)
  !    call h5dget_create_plist_f(input_dset, input_prop, hdf_err)
  !    call h5dget_space_f(input_dset, input_dspace, hdf_err)
  !    call h5sget_simple_extent_dims_f(&
  !        input_dspace, input_dims, input_maxdims, hdf_err)
  !    if(restore_idx < 0) then
  !        input_offset = [integer(hsize_t) :: 0, 0, input_dims(3) + restore_idx]
  !    else
  !        input_offset = [integer(hsize_t) :: 0, 0, restore_idx]
  !    end if
  !    call h5screate_simple_f (3, hdf5_arraychunk, input_memspace, hdf_err)
  !    call h5sselect_hyperslab_f(&
  !        input_dspace, H5S_SELECT_SET_F, input_offset, hdf5_arraychunk, hdf_err)
  !    call h5dread_f(input_dset, hdf%realtype, c_c, hdf5_arraychunk, hdf_err, &
  !        input_memspace, file_space_id=input_dspace)
  !    
  !    call h5dopen_f(input_file, phi_nbc_name, input_dset, hdf_err)
  !    call h5dget_create_plist_f(input_dset, input_prop, hdf_err)
  !    call h5dget_space_f(input_dset, input_dspace, hdf_err)
  !    call h5sget_simple_extent_dims_f(&
  !        input_dspace, input_dims, input_maxdims, hdf_err)
  !    if(restore_idx < 0) then
  !        input_offset = [integer(hsize_t) :: 0, 0, input_dims(3) + restore_idx]
  !    else
  !        input_offset = [integer(hsize_t) :: 0, 0, restore_idx]
  !    end if
  !    call h5screate_simple_f (3, hdf5_arraychunk, input_memspace, hdf_err)
  !    call h5sselect_hyperslab_f(&
  !        input_dspace, H5S_SELECT_SET_F, input_offset, hdf5_arraychunk, hdf_err)
  !    call h5dread_f(input_dset, hdf%realtype, phi_nbc, hdf5_arraychunk, hdf_err, &
  !        input_memspace, file_space_id=input_dspace)
  !    
  !    if (restore_c_only == .false.) then    
  !        call h5dopen_f(input_file, t_now_name, input_dset, hdf_err)
  !        call h5dget_create_plist_f(input_dset, input_prop, hdf_err)
  !        call h5dget_space_f(input_dset, input_dspace, hdf_err)
  !        call h5sget_simple_extent_dims_f(&
  !            input_dspace, input_dims, input_maxdims, hdf_err)
  !        input_offset = [integer(hsize_t) :: 0, 0, input_dims(3) - 1]
  !        call h5screate_simple_f (3, hdf5_scalarchunk, input_memspace, hdf_err)
  !        call h5sselect_hyperslab_f(input_dspace, H5S_SELECT_SET_F, &
  !            input_offset, hdf5_scalarchunk, hdf_err)
  !        call h5dread_f(input_dset, hdf%realtype, sim_clock%get_time(), hdf5_scalarchunk, &
  !            hdf_err, input_memspace, file_space_id=input_dspace)
  !    
  !        deallocate(nucl_node)
  !        allocate(nucl_node(nx))
  !        call h5dopen_f(input_file, nucl_node_name, input_dset, hdf_err)
  !        call h5dget_create_plist_f(input_dset, input_prop, hdf_err)
  !        call h5dget_space_f(input_dset, input_dspace, hdf_err)
  !        call h5sget_simple_extent_dims_f(&
  !            input_dspace, input_dims, input_maxdims, hdf_err)
  !        input_offset = [integer(hsize_t) :: 0, 0, input_dims(3) - 1]
  !        call h5screate_simple_f (3, hdf5_xchunk, input_memspace, hdf_err)
  !        call h5sselect_hyperslab_f(&
  !            input_dspace, H5S_SELECT_SET_F, input_offset, hdf5_xchunk, hdf_err)
  !        call h5dread_f(input_dset, H5T_NATIVE_INTEGER, nucl_node, hdf5_xchunk, &
  !            hdf_err, input_memspace, file_space_id=input_dspace)
  !        nucl_node = pack(nucl_node, nucl_node /= 0)
  !        print *, nucl_node
  !        nucl_sites = size(nucl_node)
  !        print *, nucl_sites
  !        call nucl_setmask(nucl_node, nucleation%r_nodes/dx, nucl_mask, phase_inter%get_width(), &
  !            bc_c)
  !    
  !        nucl_last_t = sim_clock%get_time()
  !        t_nextsave = sim_clock%get_time()
  !        temperature%get_value_at(sim_clock%get_time()) = get_therm(sim_clock%get_time(), temp_start, temp_rate, temp_hold, &
  !            temp_inithold)
  !    end if
  !    
  !    call h5fclose_f(input_file, hdf_err)
  !end if

  change team(core_team)
  if(team_number() == primary_team_num) then
    ! Ensure that concentrations are synchronous at the start, no matter what
    if(verbosity_level >= 5) call print_message('Broadcasting initial concentration fields...')
    call fft_field_broadcast(c_field, c_name, step)
    call c_average(c_field, species_idx, c_mean, x_mean)
    x_mean_init = x_mean
  end if
  end team


  if(this_image() == 1) then
    concentration_is_conserved = .true.
    do species = 1, n_species
      if(.not. all(c_field(species)%real%grid%axes%bound_cond%is_periodic)) then
        concentration_is_conserved = .false.
        exit
      end if
    end do
  end if
  call co_broadcast(concentration_is_conserved, 1)

  !call VTinit(fft_err)
  !print *, 'VTINIT: ', fft_err
  !call VTtraceon()

  if(.not. restore%enabled) nucleation%last_t = sim_clock%get_time()
  call simulation_stop%new_schedule(step)
  call simulation_stop%schedule(step)
  call console_output%new_schedule(step)
  call image_output%new_schedule(step)
  call hdf%output%new_schedule(step)
  call hdf%eta_output%new_schedule(step)
  call grain_reindex%timer%new_schedule(step)
  call dt%adaptive%timer%new_schedule(step)
  call nucleation%timer%new_schedule(step)
  call nucleation%timer%schedule(step)
  oscillation_limit%force = .not. oscillation_limit%is_enabled()
  oscillation_cooldown%force = .not. oscillation_cooldown%is_enabled()
  call oscillation_limit%new_schedule(step)
  call oscillation_cooldown%new_schedule(step)

  if(allocated(simulation_stop%sim_timer)) then
    total_step_factor = get_total_step_factor(sim_clock%get_time(), simulation_stop%sim_timer%get_t_trigger(), &
      species_properties, temperature, M_eta)
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Main calculation loop
  do
    ! This ensures correct entry into the adaptive time step calculation substeps
    if(step == 0 .and. dt%adaptive%enabled) then
      dt%adaptive%timer%force = .true.
    end if

    if(verbosity_level >= 5) call print_message('Starting main loop...', step)

    if(check%debug_level >= 4) then
      call c_average(c_field, species_idx, c_mean, x_mean)
      if(size(failed_images()) > 0 .or. size(stopped_images()) > 0) then
        error stop 'ERROR: Image(s) failed or stopped!'
      end if
      if(.not. co_all_equal(sim_clock%get_time())) then
        call img_printhead(step)
        print list_value_format, ': sim_time', sim_clock%get_time()*units%t, &
          'dt', dt%val*units%t, 'last_osc_step', last_osc_step, 'use_bdfab', dt%use_bdfab
        error stop 'ERROR: sim_clock%get_time() different across images!'
      else if(.not. co_all_equal(temperature%get_value_at(sim_clock%get_time()))) then
        error stop 'ERROR: temperature%get_value_at(sim_clock%get_time()) different across images!'
      else if(.not. co_all_equal(step)) then
        error stop 'ERROR: step different across images!'
      else if(.not. all(co_all_equal(x_mean))) then
        error stop 'ERROR: base phase concentration different across images!'
      end if

      change team(variable_team)
      if(team_number() == eta_team_num) then
        if(.not. co_all_equal(order_parameters%eta%interface_data%get_kappa())) &
          error stop 'ERROR: kappa different across images!'
        if(.not. co_all_equal(order_parameters%eta%interface_data%get_m_order())) &
          error stop 'ERROR: m_order different across images!'
        if(any(order_parameters%eta%real%vals < eta_min) .or. any(order_parameters%eta%real%vals > eta_max)) then
          error stop 'ERROR: eta_' // trim(phase_properties(phase_idx)%name) // '_' // convert_to_char(order_idx) //&
            ' outside valid value range (' // convert_to_char(eta_min) // ' - ' // convert_to_char(eta_max) // ')!'
        end if
        if(.not. all(ieee_is_finite(order_parameters%eta%real%vals))) error stop 'ERROR: eta contains NaN value(s)'
      else if(species_idx > 0) then
        if (any(c_field(species_idx)%real%vals < c_min) .or. any(c_field(species_idx)%real%vals > c_max)) then
          error stop 'ERROR: c_' // trim(species_properties(species_idx)%name) // ' outside valid value range ('&
            // convert_to_char(c_min) // ' - ' // convert_to_char(c_max) // ')!'
        else if(concentration_is_conserved .and. abs(x_mean(species_idx) - x_mean_init(species_idx)) > 1E-4_wp) then
          error stop 'ERROR: x_mean for ' // trim(species_properties(species_idx)%name) //&
            ' has changed significantly from its initial value! (initial: ' //&
          convert_to_char(x_mean_init(species_idx)) //&
            ', current: ' // convert_to_char(x_mean(species_idx))
        end if
        if(.not. all(ieee_is_finite(c_field(species_idx)%real%vals))) then
          error stop 'ERROR: concentration contains NaN value(s)'
        end if
      end if
      end team
    end if

    if(grain_reindex%timer%is_due(step)) then
      allocate(eta_remap(grid%n_nodes_total)[*])
      change team(phase_team)
      if(team_number() == n_phases) then
        eta_remap = order_parameters%eta%real%vals
        call eta_reassign(order_parameters%eta%real%grid, eta_remap, &
          order_parameters%eta%real%var_name, grain_reindex, eta_filled)
        order_parameters%eta%real%vals = eta_remap
        call order_parameters%eta%fft%invalidate()
      end if
      end team
      deallocate(eta_remap)
      call grain_reindex%timer%schedule(step)
      dt%bdfab_stage = 0
    end if

    if(nucleation%timer%is_due(step)) then
      change team(variable_team)
      if(team_number() == eta_team_num) then
        call order_parameters%eta_to_eta_labels()
        call nucleation%remove_expired_nuclei(step)
        if(phase_idx > 0 .and. order_idx == 1) then
          nucl_candidate_nodes = nucleation%get_nucl_candidates(phases, species_properties, phase_idx, &
            temperature%get_value_at(sim_clock%get_time()), sim_clock%get_time(), phase_inter%get_energy(), &
            df_chem_deta_common%vals, constants)
          allocate(nucl_candidate_phase_idx(size(nucl_candidate_nodes)))
          nucl_candidate_phase_idx = phase_idx
        else
          allocate(nucl_candidate_nodes(0))
          allocate(nucl_candidate_phase_idx(0))
        end if
        call concatenate_image_arrays(nucl_candidate_nodes)
        call concatenate_image_int_arrays(nucl_candidate_phase_idx)
        total_n_nucl_candidates = size(nucl_candidate_nodes)
        if(total_n_nucl_candidates > 0) then
          allocate(nucl_candidate_scramble(total_n_nucl_candidates))
          if(this_image() == 1) then
            nucl_candidate_scramble = [(i, i=1, total_n_nucl_candidates)]
            call scramble(nucl_candidate_scramble, 10)
          end if
          call co_broadcast(nucl_candidate_scramble, 1)
          nucl_candidate_phase_idx = nucl_candidate_phase_idx(nucl_candidate_scramble)
          nucl_candidate_nodes = nucl_candidate_nodes(nucl_candidate_scramble)
          call nucleation%expiration_template%new_schedule(step)
          call nucleation%expiration_template%schedule(step)
          do i = 1, size(nucl_candidate_nodes)
            node = nucl_candidate_nodes(i)
            call new_nucleus%nucleus_init(grid, nucl_candidate_phase_idx(i), nucleation%r_nodes, &
              grid%linear_to_dim_idx(node), nucleation%expiration_template, step)
            if(new_nucleus%fits_in_mask(nucleation%eta_labels)) then
              call nucleation%enqueue_nucleus(new_nucleus)
            end if
          end do
          call nucleation%calc_eta_mask(phase_idx)
          if(order_idx == 1) then
            call output_image(nucleation%eta_labels, &
              'nucl_eta_labels_' // trim(phases%phi(phase_idx)%properties%name), &
              grid%axes(1)%n_nodes)
            call output_image(nucleation%eta_mask%real%vals, &
              'nucl_eta_mask_' // trim(phases%phi(phase_idx)%properties%name), grid%axes(1)%n_nodes)
          end if
          order_parameters%eta%real%vals = min(abs(order_parameters%eta%real%vals), 1 + nucleation%eta_mask%real%vals)
          dt%bdfab_stage = 0
          call order_parameters%eta%fft%invalidate()
          call order_parameters%sqsum%invalidate()
          call order_parameters%quadsum%invalidate()
          deallocate(nucl_candidate_scramble)
        end if
        deallocate(nucl_candidate_nodes)
        deallocate(nucl_candidate_phase_idx)
      end if
      end team
      nucleation%last_t = sim_clock%get_time()
      call nucleation%timer%schedule(step)
    end if

    if(verbosity_level >= 5) call print_message('Checking for adaptive time stepping...', step)
    if(dt%adaptive%timer%is_due(step)) then
      dt%adaptive%timer%force = .false.
      if(verbosity_level >= 4) call print_message('Entering adaptive time step calculation', step)
      if(dt%adaptive%substep == 0) then
        call backup%store(step, sim_clock, c_field, phases, order_parameters, nucleation, randseed)
        step = -1
        dt%adaptive%substep = 1
      else if(dt%adaptive%substep == 1) then
        if(species_idx > 0) then
          main_var_full = c_field(species_idx)%real%vals
          surf_x_node_full = surf_x_node
        else
          main_var_full = order_parameters%eta%real%vals
        end if
        t_full = sim_clock%get_time()
        call backup%restore(step, sim_clock, c_field, phases, order_parameters, nucleation, randseed)
        ! Reset mu, i.e. make it invalid
        if(species_idx > 0) call mu_c%reset
        if(eta_idx > 0) call mu_eta%reset
        call f%reset()
        step = -dt%adaptive%n_steps_richardson
        dt%val = dt%val/dt%adaptive%n_steps_richardson
        use_bdfab_backup = dt%use_bdfab
        dt%use_bdfab = dt%adaptive%use_bdfab
        dt%bdfab_stage = 0
        dt%adaptive%substep = 2
      else if(dt%adaptive%substep == 2) then
        dt%val = dt%val*dt%adaptive%n_steps_richardson
        if(abs(t_full - sim_clock%get_time()) > 1E3*spacing(sim_clock%get_time())) then
          call img_printhead(step)
          print list_value_format, 'Time after full step t_full', t_full, 't_partial - t_full', sim_clock%get_time() - t_full
          error stop 'ERROR: Inconsistent full/partial step times in adaptive time stepping routine!'
        end if
        dt%use_bdfab = use_bdfab_backup
        if(species_idx > 0) then
          ! FIXME: reimplement non-periodic BCs
          !if(.not. c_field(species_idx)%bc%is_periodic) then
          !FIXME: more efficient cleaning?
          !FIXME: restore fourier continuation
          !call fourier_clean_var_lb(main_var_full, surf_x_node)
          !call fourier_clean_var_lb(main_var_full, surf_x_node_full)
          !call fourier_clean_var_lb(c_field(species_idx)%real%vals, surf_x_node)
          !call fourier_clean_var_lb(c_field(species_idx)%real%vals, surf_x_node_full)
          !end if
          if(species_properties(species_idx)%instant_diffusion) then
            dt%adaptive%local_err = 0
          else
            dt%adaptive%local_err = norm2_error(main_var_full, c_field(species_idx)%real%vals)
          end if
        else
          dt%adaptive%local_err = norm2_error(main_var_full, order_parameters%eta%real%vals)
        end if
        if(verbosity_level >= 5) call print_message('dt error value: ' // convert_to_char(dt%adaptive%local_err), step)
        dt%adaptive%global_err_max = max(dt%adaptive%local_err, tiny(0.0_wp))
        if(verbosity_level >= 5) call print_message('Waiting for adt error synchronization...', step)
        call co_max(dt%adaptive%global_err_max)
        if(dt%adaptive%global_err_max == dt%adaptive%local_err) call inc(dt%adaptive%n_local_limiting)
        call backup%restore(step, sim_clock, c_field, phases, order_parameters, nucleation, randseed)
        dt%adaptive%substep = 0
        dt%bdfab_stage = 0
        if(species_idx > 0) call mu_c%reset
        if(eta_idx > 0) call mu_eta%reset
        call f%reset()
        if(hdf%extended) then
          if(species_idx > 0) then
            var_error_field = c_field(species_idx)%real%vals - main_var_full
          else
            var_error_field = order_parameters%eta%real%vals - main_var_full
          end if
        end if
        if(dt%adaptive%global_err_max <= dt%adaptive%rel_err_max) then
          if(oscillation_limit%is_due(step)) then
            if(dt%use_bdfab) then
              dt%val = dt%val*min(sqrt(dt%adaptive%rel_err_tol/dt%adaptive%global_err_max), dt%adaptive%max_factor)
            else
              dt%val = dt%val*min(dt%adaptive%rel_err_tol/dt%adaptive%global_err_max, dt%adaptive%max_factor)
            end if
          else
            if(dt%use_bdfab) then
              dt%val = dt%val*min(sqrt(dt%adaptive%rel_err_tol/dt%adaptive%global_err_max), dt%adaptive%max_factor_osc)
            else
              dt%val = dt%val*min(dt%adaptive%rel_err_tol/dt%adaptive%global_err_max, dt%adaptive%max_factor_osc)
            end if
          end if
          if(allocated(hdf%output%sim_timer)) then
            if(hdf%output%sim_timer%elapsed) then
              dt%val = nearest_divisor_less_or_equal(dt%val, hdf%output%sim_timer%get_dt_trigger()) &
                + spacing(sim_clock%get_time())
            else
              dt%val = nearest_divisor_less_or_equal(dt%val, hdf%output%sim_timer%get_dt_until_trigger()) &
                + spacing(sim_clock%get_time())
            end if
          end if
          call dt%adaptive%timer%schedule(step)
        else
          dt%adaptive%n_min_factor = dt%adaptive%n_min_factor + 1
          if(verbosity_level >= 2) then
            call print_message('Maximum time step error exceeded, reducing time step from dt = ' //&
              convert_to_char(dt%val*units%t) // ' to ' // convert_to_char(dt%val*dt%adaptive%min_factor*units%t) //&
              ' (reduction count: ' // convert_to_char(dt%adaptive%n_min_factor) // ')', backup%step, 1)
          end if
          if(hdf%save_on_dt_min_factor) then
            hdf%output%force = .true.
            if(verbosity_level >= 2) call print_message('Forcing save for time stepping error!', backup%step, 1)
          end if
          dt%val = dt%val*dt%adaptive%min_factor
          cycle
        end if
        if(verbosity_level >= 5) call print_message('dt set to ' // convert_to_char(dt%val), step)
        if(dt%val < spacing(sim_clock%get_time())) then
          call error_msg(proc_name, 'time step value smaller than smallest possible simulation time increment')
        end if
        step_factor = get_step_factor(dt%val, species_properties, temperature, sim_clock%get_time(), M_eta)
        call dt%set_step_factor(step, step_factor)
        ! FIXME: this switching algorithm probably needs some work
        if(dt%bdfab_enabled .and. .not. dt%force_bdfab) then
          if(.not. dt%use_bdfab&
            .and. oscillation_limit%is_due(step)&
            .and. dt%adaptive%global_err_max >= 0.95_wp*dt%adaptive%rel_err_tol&
            .and. dt%step_factor%val < 0.9_wp*dt%step_factor%bdfab_max) then

            if(verbosity_level >= 3) call print_message('Switching bdfab on!', backup%step, 1)
            dt%use_bdfab = .true.
            dt%val = dt%val*dt%adaptive%min_factor
            dt%bdfab_stage = 0
          else if(dt%use_bdfab .and. dt%adaptive%global_err_max <= dt%adaptive%bdfab_switch_rel_error_max &
            .and. .not. oscillation_limit%is_due(step)) then
            if(verbosity_level >= 3) call print_message('Switching bdfab off', backup%step, 1)
            dt%use_bdfab = .false.
          end if
        end if
        !FIXME: rather crude - is there a better method to prevent negative eta?
        if(eta_idx > 0) then
          order_parameters%eta%real%vals = abs(order_parameters%eta%real%vals)
          call order_parameters%eta%fft%invalidate()
        end if
        if(verbosity_level >= 3 .and. this_image() == 1) then
          print list_value_format, 'Step', step, 'dt error', dt%adaptive%global_err_max, &
            'dt', dt%val*units%t, 'time', sim_clock%get_time()*units%t
        end if
      end if
    end if

    call simulation_stop%poll()
    call console_output%poll()
    call image_output%poll()
    call hdf%output%poll()
    call hdf%eta_output%poll()

    call simulation_stop%sync()
    call console_output%sync()
    call image_output%sync()
    call hdf%output%sync()
    call hdf%eta_output%sync()

    !FIXME: Probably unnecessary to recalculate the equilibrium after every step - check if x_mean changed significantly?
    ! Recalculating after every step may actually contribute to oscillations
    if(step == 0 .or. abs(temperature%get_slope_at(sim_clock%get_time())) > epsilon(0.0_wp) .or. .not. &
      concentration_is_conserved .or. any(phase_properties%calphad_use_phase_conc)) then
      if(verbosity_level >= 5) call print_message('Recalculating thermodynamic data...', step)
      call opencalphad_equilib_phase(temperature%get_value_at(sim_clock%get_time()), pressure, x_mean, phase_idx, &
        units, constants, current_equi, component_names, c_field, phases)
      if(phase_idx > 0) then
        M_eta = phases%phi(phase_idx)%properties%mobility_factor&
          *minval(species_properties%diffusion_coefficient%get_value(temperature%get_value_at(sim_clock%get_time())))
      end if
      if(verbosity_level >= 5) call print_message('Thermodynamic data calculated!', step)
    end if

    if(dt%use_bdfab .and. dt%bdfab_stage == 2 .and. step /= 0) then
      if(verbosity_level >= 5) call print_message('Doubling dt for BDF/AB', step)
      dt%val = 2*dt%val
    end if

    if(step >= 0 .and. dt%adaptive%enabled .and. allocated(fft_osc_main)) then
      if(oscillation_cooldown%is_due(step)) then
        apply_oscillation_limit = any(fft_osc_main%max_freqs%is_limiting())
      end if
    end if

    if(verbosity_level >= 5) call print_message('Waiting for oscillation data sync...', step)
    call co_any(apply_oscillation_limit)
    if(apply_oscillation_limit) then
      dt%val = dt%val/dt%adaptive%max_factor
      if(allocated(hdf%output%sim_timer)) then
        if(hdf%output%sim_timer%elapsed) then
          dt%val = nearest_divisor_less_or_equal(dt%val, hdf%output%sim_timer%get_dt_trigger()) &
            + spacing(sim_clock%get_time())
        else
          dt%val = nearest_divisor_less_or_equal(dt%val, hdf%output%sim_timer%get_dt_until_trigger()) &
            + spacing(sim_clock%get_time())
        end if
      end if
      dt%bdfab_stage = 0
      call oscillation_cooldown%schedule(step)
      call oscillation_limit%new_schedule(step)
      call oscillation_limit%schedule(step)
      apply_oscillation_limit = .false.
    end if
    if(verbosity_level >= 5) call print_message('Oscillation data synced!', step, 1)

    if(dt%use_bdfab .and. dt%bdfab_stage == 0 .and. dt%adaptive%substep /= 1) then
      if(verbosity_level >= 5) call print_message('Halving dt for BDF/AB', step)
      dt%val = 0.5_wp*dt%val
    end if

    change team(variable_team)
    if(team_number() == eta_team_num) then
      call nucleation%remove_expired_nuclei(step)
        call nucleation%calc_eta_mask(phase_idx)
        if(.not. nucleation%queue_is_empty()) then
          if(order_idx == 1) then
            order_parameters%eta%real%vals = max(abs(order_parameters%eta%real%vals), nucleation%eta_mask%real%vals)
            call order_parameters%eta%fft%invalidate()
          end if
          call order_parameters%sqsum%invalidate()
          call order_parameters%quadsum%invalidate()
        end if
      call order_parameters%calc_phi(phases)
    else
      if(.not. species_properties(species_idx)%instant_diffusion) then
        call c_field(species_idx)%forward_alldim()
      end if

      !FIXME: restore non-periodic BCs
      !if(.not. c_field(species_idx)%bc%is_periodic .and. lazy_boundary) then
      !    call fourier_cont_var_lb(c_field(species_idx), fc_data, surf_x_node)
      !end if
      !if(c_field(species_idx)%bc%is_periodic .or. lazy_boundary) then
      !    call c_field(species_idx)%forward_alldim()
      !end if
      !if(.not. c_field(species_idx)%bc%is_periodic .and. dt%prestep_c .and. .not. dt%use_bdfab) then
      !    !FIXME: balancing necessary?
      !    !call fourier_balance_var_lb(c_field(species_idx), phases%phi, surf_x_node)
      !    !FIXME: necessary?
      !    !call fourier_clean_var_lb(c_field(species_idx)%real%vals, surf_x_node)
      !end if
    end if
    end team

    if(check%debug_level >= 5) then
      call phases%phi_broadcast()
      call fft_field_broadcast(c_field, c_name, step)
      call calc_numeric_energy_derivatives(order_parameters, phases, c_field, species_idx, &
        variable_team, df_dc_numeric, df_deta_numeric)
      call phases%calc_f0_mix()
      if(species_idx > 0) then
        call phases%calc_c_mix(species_idx)
        call phases%calc_p_mix(species_idx)
        call calc_df_chem_dc(c_field, phases, species_idx, df_dc(species_idx)%chem)
      end if
    end if

    if(hdf%output%is_due(step)) then
      if(verbosity_level >= 2) call print_message('Performing save (pre-step)', step, 1)
      change team(variable_team)
      if(team_number() == eta_team_num) then
        call order_parameters%eta_to_eta_labels()
      end if
      end team
      call f%calc_total(order_parameters, phases, c_field, species_idx, variable_team)
      randseed = get_random_seed()

      sync all
      if(verbosity_level >= 5) call print_message('Save data synced!', step, 1)

      call start_critical_in_order()
      call hdf%open()
      if(this_image() == 1) then
        call output_t_now%append(sim_clock%get_time()*units%t)
        call output_therm_now%append(temperature%get_value_at(sim_clock%get_time()))
        call output_t_cpu%append(cpu_clock%get_time())
        call output_t_wall%append(wall_clock%get_time())
        call output_dt%append(dt%val*units%t)
        call output_step%append(step)
        call output_n_t_error%append(dt%adaptive%n_min_factor)
        call output_nucl_dt%append((sim_clock%get_time() - nucleation%last_t)*units%t)
        call output_f_total_sum%append(f%total_sum*units%joule)
        call output_use_bdfab%append(dt%use_bdfab)
        call output_surf_x_node%append(surf_x_node)
        call output_f_total%append(f%total%vals*units%joule/units%meter**3)
        call output_f_grad%append(f%gradient%vals*units%joule/units%meter**3)
        call output_randseed%append(randseed)
        if(hdf%extended) then
          call phases%calc_f0_mix()
        end if
      end if
      if(allocated(fft_osc_main)) then
        do i = 1, fft_osc_main%n_freqs
          call output_n_main_osc_limit(i)%append(fft_osc_main%max_freqs(i)%n_limit)
          call output_main_osc_error(i)%append(fft_osc_main%max_freqs(i)%get_error())
          call output_main_maxfreq(i)%append(fft_osc_main%max_freqs(i)%prev_maxfreq)
        end do
      end if
      if(species_idx > 0) then
        call output_c%append(c_field(species_idx)%real%vals*units%c)
        call output_c_error%append(dt%adaptive%local_err)
        call output_n_c_limit%append(n_var_dt_limit)
        call output_D%append(species_properties(species_idx)%diffusion_coefficient% &
          get_value(temperature%get_value_at(sim_clock%get_time()))*units%D)
        if(hdf%extended) call output_c_error_field%append(var_error_field*units%c)
      end if
      if(eta_idx > 0) then
        if(eta_idx == 1) then
          call output_if_energy%append(phase_inter%get_energy()*units%if_energy)
          call output_kappa%append(order_parameters%eta%interface_data%get_kappa()*units%kappa)
          call output_m_order%append(order_parameters%eta%interface_data%get_m_order()*units%m_order)
          call output_eta_idx%append(order_parameters%eta_labels%vals)
          call output_sqsum%append(order_parameters%sqsum%vals)
        end if
        call output_eta_error%append(dt%adaptive%local_err)
        call output_n_eta_limit%append(n_var_dt_limit)
        if(order_idx == 1) then
          call output_M_eta%append(M_eta*units%M_eta)
          if(phase_idx < n_phases) call output_phi%append(phases%phi(phase_idx)%vals)
          do i = 1, n_species
            call output_c_equi(i)%append(phases%phi(phase_idx)%equi_val(i)*units%c)
            call output_p_coeff(i)%append(phases%phi(phase_idx)%p_coeff(i)*units%p_coeff)
            call output_k_coeff(i)%append(phases%phi(phase_idx)%k_coeff(i)*units%k_coeff)
            call output_omega(i)%append(phases%phi(phase_idx)%omega(i)*units%omega)
          end do
        end if
        if(hdf%extended) then
          call output_eta%append(order_parameters%eta%real%vals)
          call output_eta_error_field%append(var_error_field)
        else if(hdf%eta_output%is_due(step)) then
          call print_message('Saving etas...', step, 1)
          call hdf%eta_output%schedule(step)
          call output_eta%append(order_parameters%eta%real%vals)
          if(this_image() == 1) call output_backup_step%append(step)
        end if
      end if
      call hdf%close()
      call millisleep(50)
      call end_critical_in_order()
      sync all
    end if

    if(dt%prestep_c .and. .not. dt%use_bdfab) then
      change team(core_team)
      if(team_number() == primary_team_num) then
        call phases%phi_broadcast()
      end if
      if(species_idx > 0) then
        !call set_boundaries(c_field(species_idx), phases%phi, fc_data, main_fft_now, fft, n_c_fft, surf_x_node, t_coa)
        if(species_properties(species_idx)%instant_diffusion) then
          call phases%calc_c_mix(species_idx)
          call phases%calc_p_mix(species_idx)
          mu_average = (c_mean(species_idx)&
            - mean(phases%c_mix(species_idx)%real%vals))/mean(phases%p_mix(species_idx)%real%vals)
          do concurrent(phase = 1:n_phases)
            c_target(phase) = mu_average/phases%phi(phase)%p_coeff(species_idx)&
              + phases%phi(phase)%equi_val(species_idx)
          end do
          do concurrent(node = 1:grid%n_nodes_total)
            c_field(species_idx)%real%vals(node) = phases%phi_mix(c_target, node)
          end do
          c_field(species_idx)%real%time = c_field(species_idx)%real%time + dt%val
        else
          call calc_mu_c(c_field(species_idx), phases, species_idx, mu_c, dt%use_bdfab, dt%bdfab_stage)
          call spectral_step_c(c_field(species_idx), mu_c, species_properties(species_idx)% &
            diffusion_coefficient%get_value(temperature%get_value_at(sim_clock%get_time())), &
            dt%val, dt%use_bdfab, dt%bdfab_stage, phases%equal_p_coeff(species_idx))
        end if
        !FIXME: optimize cleaning positions
        !if(.not. c_field(species_idx)%bc%is_periodic) call fourier_clean_var_lb(c_field(species_idx)%real%vals, surf_x_node)
      else if(phase_idx > 0) then
        if(order_parameters%eta%is_active()) then
          call order_parameters%eta%forward_alldim()
          call calc_df_order_deta(order_parameters, df_deta%order)
          if(order_idx == 1) then
            if(.not. all(grid%axes%bound_cond%is_periodic)) then
              call find_surface(grid, phases%phi(1)%vals, surf_x_node)
            end if
            call phases%calc_f0_mix()
            do concurrent(species = 1:n_species)
              call phases%calc_c_mix(species)
              call phases%calc_p_mix(species)
            end do
          end if
        end if
      end if
      if(team_number() == primary_team_num) then
        call fft_field_broadcast(c_field, c_name, step)
        if(.not. concentration_is_conserved) then
          call c_average(c_field, species_idx, c_mean, x_mean)
        end if
      end if
      end team
      change team(phase_team)
      if(team_number() /= c_team_num) then
        !FIXME: prestep_c argument should use a named variable/constant
        if(this_image() == phase_team_main_image) then
          call calc_df_chem_deta_common(&
            c_field, order_parameters, species_properties, phases, phase_idx, .true., dt%val, df_chem_deta_common)
        end if
        if(num_images() > 1) then
          call real_broadcast(df_chem_deta_common%vals, mu_name, step)
          call co_broadcast(df_chem_deta_common%time, phase_team_main_image)
          call real_broadcast(x_mean, 'x_mean', step)
        end if
        if(order_parameters%eta%is_active()) then
          call calc_df_chem_deta(order_parameters, df_chem_deta_common, df_deta%chem)
          call spectral_step_eta(order_parameters, df_deta, mu_eta, M_eta, dt%val, dt%use_bdfab, dt%bdfab_stage)
        else
          order_parameters%eta%real%time = order_parameters%eta%real%time + dt%val
        end if
      end if
      end team
    else
      change team(core_team)
      if(team_number() == primary_team_num) then
        call phases%phi_broadcast()
        call fft_field_broadcast(c_field, c_name, step)
        if(.not. concentration_is_conserved) then
          call c_average(c_field, species_idx, c_mean, x_mean)
        end if
      else if(team_number() == secondary_team_num .and. order_parameters%eta%is_active()) then
        call order_parameters%eta%forward_alldim()
        call calc_df_order_deta(order_parameters, df_deta%order)
      end if
      end team
      !FIXME: There is likely room for optimization here
      change team(phase_team)
      if(team_number() /= c_team_num) then
        if(this_image() == phase_team_main_image) then
          call phases%calc_f0_mix()
          do concurrent(species = 1:n_species)
            call phases%calc_c_mix(species)
            call phases%calc_p_mix(species)
          end do
          if(.not. all(grid%axes%bound_cond%is_periodic)) then
            call find_surface(grid, phases%phi(1)%vals, surf_x_node)
          end if
          !FIXME: dt%prestep_c argument should use a named variable/constant
          call calc_df_chem_deta_common(&
            c_field, order_parameters, species_properties, phases, phase_idx, .false., dt%val, df_chem_deta_common)
        end if
        if(num_images() > 1) then
          call real_broadcast(df_chem_deta_common%vals, mu_name, step)
          call co_broadcast(df_chem_deta_common%time, phase_team_main_image)
          call real_broadcast(x_mean, 'x_mean', step)
        end if
        if(this_image() == phase_team_main_image) then
          call order_parameters%eta%forward_alldim()
          call calc_df_order_deta(order_parameters, df_deta%order)
        end if
        if(order_parameters%eta%is_active()) then
          call calc_df_chem_deta(order_parameters, df_chem_deta_common, df_deta%chem)
          call spectral_step_eta(order_parameters, df_deta, mu_eta, M_eta, dt%val, dt%use_bdfab, dt%bdfab_stage)
        else
          order_parameters%eta%real%time = order_parameters%eta%real%time + dt%val
        end if
      else if(species_idx > 0) then
        !call set_boundaries(c_field(species_idx), phases%phi, fc_data, main_fft_now, fft, n_c_fft, surf_x_node, t_coa)
        if(species_properties(species_idx)%instant_diffusion) then
          call phases%calc_c_mix(species_idx)
          call phases%calc_p_mix(species_idx)
          mu_average = (c_mean(species_idx)&
            - mean(phases%c_mix(species_idx)%real%vals))/mean(phases%p_mix(species_idx)%real%vals)
          do concurrent(phase = 1:n_phases)
            c_target(phase) = mu_average/phases%phi(phase)%p_coeff(species_idx)&
              + phases%phi(phase)%equi_val(species_idx)
          end do
          do concurrent(node = 1:grid%n_nodes_total)
            c_field(species_idx)%real%vals(node) = phases%phi_mix(c_target, node)
          end do
          c_field(species_idx)%real%time = c_field(species_idx)%real%time + dt%val
        else
          call calc_mu_c(c_field(species_idx), phases, species_idx, mu_c, dt%use_bdfab, dt%bdfab_stage)
          call spectral_step_c(c_field(species_idx), mu_c, species_properties(species_idx)% &
            diffusion_coefficient%get_value(temperature%get_value_at(sim_clock%get_time())), &
            dt%val, dt%use_bdfab, dt%bdfab_stage, phases%equal_p_coeff(species_idx))
        end if
        !FIXME: optimize cleaning positions
        !if(.not. c_field(species_idx)%bc%is_periodic) call fourier_clean_var_lb(c_field(species_idx)%real%vals, surf_x_node)
      end if
      end team
    end if

    if(check%debug_level >= 5) then
      if(species_idx > 0) then
        call df_dc(species_idx)%chem%output_image()
        call df_dc_numeric(species_idx)%chem%output_image()
        print list_value_format, 'df_chem_dc error for ' // species_properties(species_idx)%name, &
          norm2_error(df_dc(species_idx)%chem%vals, df_dc_numeric(species_idx)%chem%vals)
      else if(phase_idx > 0) then
        if(order_parameters%eta%is_active()) then
          print list_value_format, 'df_order_deta error for eta ' // convert_to_char(eta_idx), &
            norm2_error(df_deta%order%vals, df_deta_numeric%order%vals)
          print list_value_format, 'df_chem_deta error for eta ' // convert_to_char(eta_idx), &
            norm2_error(df_deta%chem%vals, df_deta_numeric%chem%vals)
          call df_deta%order%output_image()
          call df_deta_numeric%order%output_image()
          call df_deta%chem%output_image()
          call df_deta_numeric%chem%output_image()
          if(order_idx == 1) call df_chem_deta_common%output_image()
        end if
      end if
    end if

    if(hdf%output%is_due(step)) then
      call hdf%output%schedule(step)
      if(verbosity_level >= 2) call print_message('Performing save (post-step)', step, 1)

      if(hdf%extended) then
        call start_critical_in_order()
        call hdf%open()
        if(species_idx > 0) then
          if(.not. species_properties(species_idx)%instant_diffusion) call mu_c%backward_alldim()
          call output_mu_c%append(mu_c%real%vals*units%mu_c)
        end if
        if(eta_idx > 0) then
          if(order_idx == 1) call output_mu_eta_common%append(df_chem_deta_common%vals*units%df_chem_deta)
          call output_mu_eta%append(mu_eta%real%vals*units%df_chem_deta)
        end if
        call hdf%close()
        call millisleep(50)
        call end_critical_in_order()
      end if

      n_saved = n_saved + 1
    end if

    if(this_image() == 1 .and. verbosity_level >= 4) then
      print '(a)', 'Current parameter values (scaled):'
      print list_value_format, kappa_name, phase_inter%get_kappa(), m_order_name, phase_inter%get_m_order()
      do j = 1, n_species
        print '(a)', species_properties(j)%name // ':'
        print list_value_format, D_name, species_properties(j)%diffusion_coefficient%get_value(&
          temperature%get_value_at(sim_clock%get_time()))
      end do
      do i = 1, n_phases
        do j = 1, n_species
          print '(a)', species_properties(j)%name // ' in ' // phase_properties(i)%name // ':'
          print list_value_format, c_equi_name, phases%phi(i)%equi_val(j), p_coeff_name, phases%phi(i)%p_coeff(j)
        end do
      end do
    end if
    if(verbosity_level >= 5) then
      call img_printhead(step)
      print list_value_format, 'adding dt', dt%val*units%t, 'to sim_clock%get_time()', sim_clock%get_time()*units%t
    end if
    call sim_clock%increment_time(dt%val)

    if(allocated(fft_osc_main) .and. dt%adaptive%substep == 0) then
      if(species_idx > 0) then
        call fft_osc_main%update_all(c_field(species_idx)%fft%vals)
      end if
      if(eta_idx > 0) then
        call fft_osc_main%update_all(order_parameters%eta%fft%vals)
      end if
    end if

    if(console_output%is_due(step)) then
      call console_output%schedule(step)
      if(verbosity_level >= 5) call print_message('Entering printing routine...', step)
      if(.not. co_all_equal(sim_clock%get_time())) then
        print *, this_image(), sim_clock%get_time()*units%t
        error stop 'ERROR: sim_clock%get_time() different across images!'
      else if(.not. co_all_equal(temperature%get_value_at(sim_clock%get_time()))) then
        error stop 'ERROR: temperature%get_value_at(sim_clock%get_time()) different across images!'
      else if(.not. co_all_equal(step)) then
        error stop 'ERROR: step different across images!'
      end if
      steps_per_second = (step - step_last_print)&
        /max((wall_clock%get_time() - t_last_print), epsilon(0.0_wp))
      t_last_print = wall_clock%get_time()
      step_last_print = step
      call f%calc_total(order_parameters, phases, c_field, species_idx, variable_team)
      if(verbosity_level >= 1) then
        call print_info(step, cpu_clock%get_time(), sim_clock%get_time(), wall_clock%get_time(), &
          temperature%get_value_at(sim_clock%get_time()), steps_per_second, simulation_stop, dt, f, order_parameters, &
          total_step_factor, variable_team, units)
      end if
      if(verbosity_level >= 2) call print_extra_info(dt, c_field, order_parameters, fft_osc_main, species_idx, phase_idx)
      !FIXME: restore checking for increasing total free energy
      !if(abs(temperature%get_value_at(sim_clock%get_time()) - temp_hold) < delta_num .and. f%total_sum > f%total_sum_prev) then
      !    call print_message('WARNING: Total free energy increased! (difference: ' //&
      !        convert_to_char(f%total_sum - f%total_sum_prev) // ')', step, 1)
      !end if
    end if

    if(image_output%is_due(step)) then
      call image_output%schedule(step)
      if(verbosity_level >= 5) call print_message('Entering image output routine...', step)
      call f%calc_total(order_parameters, phases, c_field, species_idx, variable_team)
      change team(variable_team)
      if(team_number() == eta_team_num) then
        call order_parameters%eta_to_eta_labels()
        if(eta_idx == 1) then
          call f%total%output_image()
          call order_parameters%sqsum%output_image()
          call order_parameters%eta_labels%output_image()
          call phases%output_composite_image(order_parameters%sqsum)
        end if
        call order_parameters%eta%real%output_image(a_min_opt = 0.0_wp, a_max_opt = 1.0_wp)
        if(order_idx == 1 .and. phase_idx < n_phases) then
          call phases%phi(phase_idx)%output_image(a_min_opt = 0.0_wp, a_max_opt = 1.0_wp)
        end if
        if(hdf%extended .and. order_parameters%eta%is_active()) then
          call mu_eta%real%output_image()
        end if
      else
        call c_field(species_idx)%real%output_image()
        if(hdf%extended) then
          if(.not. species_properties(species_idx)%instant_diffusion) then
            call mu_c%backward_alldim()
            call mu_c%real%output_image()
          end if
        end if
      end if
      end team
    end if

    if(simulation_stop%is_enabled()) then
      if(simulation_stop%is_due(step)) exit
    end if
    if(dt%bdfab_stage /= 1 .or. dt%adaptive%substep == 1 .or. .not. dt%use_bdfab) step = step + 1
    if(dt%use_bdfab .and. dt%bdfab_stage < 3) dt%bdfab_stage = dt%bdfab_stage + 1
  end do
  ! end of main calculation loop
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(this_image() == 1 .and. hdf%output%is_enabled()) then
    call h5close_f(hdf_err)
  end if

  !call VTtraceoff()
  !call VTfini(fft_err)

  print *, 'Image', this_image(), 'stopping.'
  sync all
  stop
end program campfase
