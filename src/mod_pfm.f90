! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_pfm

  use mod_base_functions
  use mod_globals, only: verbosity_level
  use mod_parameters
  use mod_regular_grid
  use mod_graphics
  use mod_time_dep_vars
  use mod_particle
  use mod_fcdata
  use mod_normal_dist
  use mod_spectral
  use mod_fft_osc
  use mod_mu
  use mod_phase_set
  use mod_order_parameter_set
  use mod_phase_field_energies
  use mod_time_step
  use mod_species_properties
  use mod_stopwatch
  use mod_unit_factors
  use mod_remap_data
  implicit none
  private

  public :: print_info, print_extra_info, eta_reassign, find_surface, c_average, calc_mu_c, spectral_step_c, spectral_step_eta, &
    real_broadcast, get_step_factor, get_total_step_factor

contains

  subroutine print_info(step, t_cpu, t_sim, t_wall, temp_now, steps_per_second, simulation_stop, dt, f, order_parameters, &
      total_step_factor, variable_team, units)
    integer(int64), intent(in) :: step
    type(team_type), intent(in) :: variable_team
    real(wp), intent(in) :: t_cpu, t_sim, t_wall, temp_now, steps_per_second, total_step_factor
    type(order_parameter_set), intent(in) :: order_parameters
    type(timer_set), intent(in) :: simulation_stop
    type(time_step), intent(in) :: dt
    type(phase_field_energies), intent(in) :: f
    type(unit_factor_t), intent(in) :: units
    integer :: n_eta_active
    integer, parameter :: t_limitfactor = 5
    integer, parameter :: min_s = 60
    integer, parameter :: h_s = 60*min_s
    integer, parameter :: d_s = 24*h_s
    integer, parameter :: w_s = 7*d_s
    real(wp) :: t_remain
    character(3) :: t_unit_char

    if(this_image() == 1) then
      write(output_unit, '(a,i0)', advance = 'no') 'Step number ', step
      if(allocated(simulation_stop%step_timer)) then
        write(output_unit, '(a,i0)', advance = 'no') ' of ', simulation_stop%step_timer%get_step_trigger()
      end if
      call newline()
      write(output_unit, '(a)', advance = 'no') 'Time stepping: '
      if(dt%use_bdfab) then
        write(output_unit, '(a)') '2nd order Backward differentiation formula/Adams-Bashforth'
      else
        write(output_unit, '(a)') '1st order semi-implicit Backward Euler'
      end if
      write(output_unit, '(2(a,es10.3),a)') 'Current time step size: ', dt%val*units%t, ', ', steps_per_second, ' steps per second'
      write(output_unit, '(a,es10.3,a,f7.2,a)') 'Simulation time: ', t_sim*units%t, ' s, Temperature: ', temp_now, ' K'
      write(output_unit, '(2(a,es10.3),a)') 'CPU time: ', t_cpu, ' s, wall time: ', t_wall, ' s'
    end if
    sync all
    call millisleep(10)

    change team(variable_team)
    if(team_number() == eta_team_num) then
      n_eta_active = order_parameters%eta%get_n_eta_active()
      if(order_parameters%eta%idx == 1) then
        write(output_unit, '(a,i0,a,es16.9)') 'Active eta fields: ', n_eta_active, ', total free energy: ', f%total_sum
      end if
    end if
    end team
    sync all
    call millisleep(10)

    if(this_image() == 1 .and. step > 0 .and. simulation_stop%is_enabled()) then
      t_remain = huge(0.0_wp)
      if(allocated(simulation_stop%wall_timer)) then
        t_remain = min(t_remain, simulation_stop%wall_timer%get_dt_until_trigger())
      end if
      if(allocated(simulation_stop%step_timer)) then
        t_remain = min(t_remain, t_wall*(real(simulation_stop%step_timer%get_step_trigger(),wp)/step - 1))
      end if
      if(allocated(simulation_stop%sim_timer)) then
        t_remain = min(t_remain, t_wall*(total_step_factor/dt%step_factor%mean/step - 1))
      end if
      if(t_remain <= t_limitfactor*min_s) then
        t_unit_char = 's'
      else if(t_remain <= t_limitfactor*h_s) then
        t_remain = t_remain/min_s
        t_unit_char = 'min'
      else if(t_remain <= t_limitfactor*d_s) then
        t_remain = t_remain/h_s
        t_unit_char = 'h'
      else if(t_remain <= t_limitfactor*w_s) then
        t_remain = t_remain/d_s
        t_unit_char = 'd'
      else if(t_remain < 1E4*w_s) then
        t_remain = t_remain/w_s
        t_unit_char = 'w'
      else
        t_remain = 9999.9_wp
        t_unit_char = 'w'
      end if
      write(output_unit, '(a,f6.1,a,a)') 'Estimated remaining time: ', t_remain, ' ', t_unit_char
      call newline()
    end if
  end subroutine print_info

  subroutine print_extra_info(dt, c_field, order_parameters, fft_osc, species_idx, phase_idx)

    type(time_step), intent(in) :: dt
    type(fft_field), intent(in) :: c_field(:)
    type(order_parameter_set), intent(in) :: order_parameters
    type(fft_grid_max_freq_oscillations), allocatable, intent(in) :: fft_osc
    integer, intent(in) :: species_idx, phase_idx
    character(:), allocatable :: var_name, var_name_format
    integer :: i, max_var_name_len, osc_count_chars, osc_limit_chars, osc_count_pad, osc_limit_pad, osc_image
    integer(int64) :: n_fft_forward, n_fft_backward
    character(*), parameter :: osc_separator = '/'
    character(*), parameter :: proc_name = 'print_extra_info'
    integer, parameter :: osc_count_digits = 2
    integer, parameter :: osc_limit_digits = 5
    integer, parameter :: osc_separator_len = len(osc_separator)
    logical :: fft_osc_enabled, any_fft_osc_enabled

    if(.not. dt%adaptive%enabled) return
    if(species_idx > 0) then
      var_name = trim(c_field(species_idx)%real%var_name)
      n_fft_forward = c_field(species_idx)%n_forward_alldim
      n_fft_backward = c_field(species_idx)%n_backward_alldim
    else if(phase_idx > 0) then
      var_name = trim(order_parameters%eta%real%var_name)
      n_fft_forward = order_parameters%eta%n_forward_alldim
      n_fft_backward = order_parameters%eta%n_backward_alldim
    else
      allocate(character(0) :: var_name)
    end if
    max_var_name_len = len(var_name)
    call co_max(max_var_name_len)
    var_name_format = '(a' // convert_to_char(max_var_name_len) // ')'

    osc_image = 0
    fft_osc_enabled = allocated(fft_osc)
    if(fft_osc_enabled) then
      osc_image = this_image()
      osc_count_chars = (osc_separator_len + osc_count_digits)*fft_osc%n_freqs - osc_separator_len
      osc_limit_chars = (osc_separator_len + osc_limit_digits)*fft_osc%n_freqs - osc_separator_len
      osc_count_pad = osc_count_chars - 5
      osc_limit_pad = osc_limit_chars - 10
    end if
    any_fft_osc_enabled = fft_osc_enabled
    call co_any(any_fft_osc_enabled)
    if(any_fft_osc_enabled) then
      call co_min(osc_image)
      call co_broadcast(osc_count_chars, osc_image)
      call co_broadcast(osc_count_pad, osc_image)
      call co_broadcast(osc_limit_chars, osc_image)
      call co_broadcast(osc_limit_pad, osc_image)
    end if

    if(.not. co_all_equal(dt%adaptive%n_min_factor)) then
      call error_msg(proc_name, 'Number of t reductions not equal across images!')
    end if

    if(this_image() == 1) then
      write(output_unit, '(4(a,es10.3),a)') 'Stepping factor: Current:', dt%step_factor%val, ', mean: ', dt%step_factor%mean, &
        ', max: ', dt%step_factor%max, ' (BDFAB: ', dt%step_factor%bdfab_max, ')'
      write(output_unit, '(a)') 'Adaptive timestep calculation:'
      write(output_unit, '(a,i0)') 'Number of timestep reductions: ', dt%adaptive%n_min_factor
      write (output_unit, '(a,a)', advance = 'no') repeat(' ', max_var_name_len), &
        '     dt_error  n_fft_forward  n_fft_backward  n_dt_limit'
      if(any_fft_osc_enabled) then
        write (output_unit, '(a)', advance = 'no') '  osc: '
        if(osc_count_pad > 0) write (output_unit, '(a)', advance = 'no') repeat(' ', osc_count_pad)
        write (output_unit, '(a)', advance = 'no') 'count  '
        if(osc_limit_pad > 0) write (output_unit, '(a)', advance = 'no') repeat(' ', osc_limit_pad)
        write (output_unit, '(a)', advance = 'no') 'n_limiting'
      end if
      call newline()
    end if
    call start_critical_in_order()
    write (output_unit, var_name_format, advance = 'no') var_name
    write (output_unit, '(a,es12.5,a,i13,a,i14,a,i10)', advance = 'no') &
      ' ', dt%adaptive%local_err, '  ', n_fft_forward, '  ', n_fft_backward, '  ', dt%adaptive%n_local_limiting
    if(fft_osc_enabled) then
      write (output_unit, '(a)', advance = 'no') '       '
      if(osc_count_pad < 0) write (output_unit, '(a)', advance = 'no') repeat(' ', -osc_count_pad)
      do i = 1, fft_osc%n_freqs
        write (output_unit, '(i' // convert_to_char(osc_count_digits) // ')', advance = 'no') &
          abs(fft_osc%max_freqs(i)%n_osc)
        if(i < fft_osc%n_freqs) write (output_unit, '(a)', advance = 'no') osc_separator
      end do
      write (output_unit, '(a)', advance = 'no') '  '
      if(osc_limit_pad < 0) write (output_unit, '(a)', advance = 'no') repeat(' ', -osc_limit_pad)
      do i = 1, fft_osc%n_freqs
        write (output_unit, '(i' // convert_to_char(osc_limit_digits) // ')', advance = 'no') &
          fft_osc%max_freqs(i)%n_limit
        if(i < fft_osc%n_freqs) write (output_unit, '(a)', advance = 'no') osc_separator
      end do
    end if
    call newline()
    call millisleep(5)
    call end_critical_in_order()
    sync all
    if(any_fft_osc_enabled .and. this_image() == osc_image) then
      write (output_unit, '(a)', advance = 'no') '(Oscillation frequencies: '
      do i = 1, fft_osc%n_freqs
        write (output_unit, '(a)', advance = 'no') fft_osc%max_freqs(i)%name
        if(i < fft_osc%n_freqs) write (output_unit, '(a)', advance = 'no') osc_separator
      end do
      write (output_unit, '(a)', advance = 'yes') ')'
    end if

    call millisleep(5)
    sync all

    if(this_image() == num_images()) then
      call newline()
      call newline()
    end if
  end subroutine print_extra_info

  !FIXME: eta remapping currently assumes eta_idx == this_image() - this could be dangerous
  subroutine get_remap_grains(grid, eta, eta_name, settings, eta_remap_data, allow_performance_remap)
    type(regular_grid), intent(in) :: grid
    real(wp), intent(in) :: eta(:)
    character(*), intent(in) :: eta_name
    type(grain_reindex_settings), intent(in) :: settings
    logical, intent(in) :: allow_performance_remap
    type(remap_data), intent(inout) :: eta_remap_data
    integer :: i, ngrains, nlabels, n_remap, remap_idx, label_calls,&
      n_remap_merge, perf_idx, nlabels_target, grain_1, grain_2, grain_select
    integer, allocatable, dimension(:) :: merge_grain_1, merge_grain_2,&
      core_nodes, grain_overlaps, grain_overlaps_lower, grain_nnodes, labels
    real(wp), allocatable, dimension(:) :: merge_seps
    real(wp) :: minsep, minsep_upper, minsep_lower, a_fill_ratio
    logical, allocatable, dimension(:) :: fill_mask, core_mask, grain_mask
    logical :: performance_remap
    character(*), parameter :: proc_name = 'get_remap_grains'

    if(verbosity_level >= 3) call print_message('Entering get_remap_grains', co_img_opt = 1)
    allocate(labels(grid%n_nodes_total))

    label_calls = 0
    remap_idx = 0
    core_mask = eta > settings%separators%intra
    grain_mask = eta > settings%separators%inter
    fill_mask = eta > settings%separators%fill
    a_fill_ratio = real(count(fill_mask),wp)/grid%n_nodes_total

    ! Find the core parts of all grains - these are used both to count the
    ! number of actual grains and determine if a coherent region is part of
    ! a grain - labels not containing any core parts are discarded.
    ! The intra-grain separator defines indisputably coherent regions - should
    ! two merging grains have a larger separator, they can no longer be
    ! distinguished at this point.
    call label_areas(grid, core_mask, labels, ngrains)
    label_calls = label_calls + 1
    if(verbosity_level >= 1) then
      call output_image(labels, 'remap_' // eta_name // '_core_labels', grid%axes(1)%n_nodes)
      call output_image(core_mask, 'remap_' // eta_name // '_core_grains', grid%axes(1)%n_nodes)
    end if

    ! Find and store the first nodes of all grains for later identification of labels
    allocate(core_nodes(ngrains))
    allocate(grain_overlaps(ngrains))
    allocate(grain_overlaps_lower(ngrains))
    allocate(grain_nnodes(ngrains))
    grain_overlaps = 0
    grain_overlaps_lower = 0
    grain_nnodes = 0
    do i = 1, ngrains
      core_nodes(i) = first_true(labels == i)
      grain_nnodes(i) = count(labels == i)
    end do

    ! Find the grains as determined by the inter-grain separator
    call label_nodes(grid, grain_mask, core_nodes, labels, nlabels, grain_overlaps)
    label_calls = label_calls + 1
    if(verbosity_level >= 1) then
      call output_image(labels, 'remap_' // eta_name // '_grain_labels', grid%axes(1)%n_nodes)
      call output_image(grain_mask, 'remap_' // eta_name // '_grain_grains', grid%axes(1)%n_nodes)
    end if

    ! Grains with risk of merging are identified by a seperator between the
    ! intra- and inter-grain separators - therefore, the number of merges is
    ! equal to the difference between the number of labels when using the two
    ! separators.
    n_remap_merge = max(ngrains - nlabels, 0)
    perf_idx = n_remap_merge

    ! Add all grains for performance remapping if fill percentage is below threshold
    ! Performance remapping is only an option if no grains with merging risk exist
    performance_remap = allow_performance_remap &
      .and. a_fill_ratio <= settings%for_performance%fill_ratio &
      .and. n_remap_merge == 0

    ! If merging grains have been identified, find their respective separators
    if(n_remap_merge > 0) then
      ! Allocate for pairs of grains which are candidates for remapping due
      ! to merging risk. The actual grains to be remapped will be selected
      ! according to criteria defined below.
      allocate(merge_seps(n_remap_merge))
      allocate(merge_grain_1(n_remap_merge))
      allocate(merge_grain_2(n_remap_merge))
      merge_seps = 0.0_wp
      merge_grain_1 = 0
      merge_grain_2 = 0
      minsep_upper = settings%separators%intra
      ! Iteratively search for the minimum separators between merging grains
      ! up to a defined relative tolerance.
      do remap_idx = 1, n_remap_merge
        minsep_lower = settings%separators%inter
        nlabels_target = ngrains - remap_idx + 1
        grain_overlaps_lower = 0
        do
          minsep = (minsep_upper + minsep_lower)/2
          call label_nodes(grid, eta > minsep, core_nodes, labels, nlabels, grain_overlaps)
          label_calls = label_calls + 1
          if(nlabels == nlabels_target) then
            minsep_upper = minsep
            if((1 - minsep_lower/minsep_upper) < settings%separators%rel_tol .and. any(grain_overlaps_lower > 0)) exit
          else if(nlabels < nlabels_target) then
            grain_overlaps_lower = grain_overlaps
            minsep_lower = minsep
          else
            ! In rare cases (when minimum separators for two different
            ! sets of merging grains are very similar), nlabels will
            ! be greater than nlabels_target, the target number of
            ! labels
            minsep_upper = minsep
          end if
        end do
        ! Store the determined minimum separator
        merge_seps(remap_idx) = minsep
        ! Identify the grains that merge near the given separator
        merge_grain_2(remap_idx) = first_true(grain_overlaps_lower > 0 .and. grain_overlaps == 0)
        merge_grain_1(remap_idx) = grain_overlaps_lower(merge_grain_2(remap_idx))
        if(verbosity_level >= 1) then
          call img_printhead()
          write(output_unit, '(3(a,i0),2(a,es12.5))') &
            'Found merge ', remap_idx, ' out of ', n_remap_merge, ' after ', label_calls, &
            ' label calls, separator = ', minsep_lower, ' - ', minsep_upper
        end if
      end do
      if(verbosity_level >= 1) then
        call img_printhead()
        write(output_unit, '(a)', advance = 'no') 'Merge detection completed.'
        call newline()
      end if
      ! Grain remapping priority:
      !   - only choose grains which have no merges at higher separators
      !   - in merging pairs, prefer grain with the smaller amount of nodes
      !   - multiply overlapping grains should be rare, so treat them equally

      ! Exclude grains occuring in multiple merging pairs so that only the 
      ! ones with the largest separator remains as a valid option
      do i = 1, n_remap_merge - 1
        grain_1 = merge_grain_1(i)
        grain_2 = merge_grain_2(i)
        where(merge_grain_1(i+1:n_remap_merge) == grain_1 .or. merge_grain_1(i+1:n_remap_merge) == grain_2)
          merge_grain_1(i+1:n_remap_merge) = -1
        end where
        where(merge_grain_2(i+1:n_remap_merge) == grain_1 .or. merge_grain_2(i+1:n_remap_merge) == grain_2)
          merge_grain_2(i+1:n_remap_merge) = -1
        end where
      end do
      ! Some merging pairs might be too complex to resolve according to
      ! the criteria above, so discard them and hope for options later on
      n_remap_merge = n_remap_merge - count(merge_grain_1 == -1 .and. merge_grain_2 == -1)
      perf_idx = n_remap_merge
      if(performance_remap) then
        n_remap = ngrains
      else
        n_remap = n_remap_merge
      end if
      if(verbosity_level >= 1) then
        call img_printhead()
        call newline()
        write(output_unit, '(a)') 'Merge grain 1:   '
        do i = 1, size(merge_grain_1)
          write (output_unit, '(i8)', advance = 'no') merge_grain_1(i)
          write(output_unit, '(a)', advance = 'no') ' '
        end do
        call newline()
        write(output_unit, '(a)') 'Merge grain 2:   '
        do i = 1, size(merge_grain_2)
          write (output_unit, '(i8)', advance = 'no') merge_grain_2(i)
          write(output_unit, '(a)', advance = 'no') ' '
        end do
        call newline()
        write(output_unit, '(a)') 'Merge nodes 1:'
        do i = 1, size(merge_grain_1)
          if(merge_grain_1(i) > 0 .and. merge_grain_1(i) <= size(core_nodes)) then
            write (output_unit, '(i8)', advance = 'no') core_nodes(merge_grain_1(i))
          else
            write(output_unit, '(a)', advance = 'no') 'INVALID '
          end if
          write(output_unit, '(a)', advance = 'no') ' '
        end do
        call newline()
        write(output_unit, '(a)') 'Merge nodes 2:'
        do i = 1, size(merge_grain_2)
          if(merge_grain_2(i) > 0 .and. merge_grain_2(i) <= size(core_nodes)) then
            write (output_unit, '(i8)', advance = 'no') core_nodes(merge_grain_2(i))
          else
            write(output_unit, '(a)', advance = 'no') 'INVALID '
          end if
          write(output_unit, '(a)', advance = 'no') ' '
        end do
        call newline()
      end if
      eta_remap_data = remap_data(n_remap, num_images())
      eta_remap_data%grains(1:n_remap_merge)%merge_risk = .true.
      eta_remap_data%grains(1:n_remap_merge)%separator = settings%separators%inter
      remap_idx = 1
      ! Choose grains for remapping from the pairs determined above
      ! If only one of the grains is a valid option, choose it
      ! Else, select the grain with the smaller size, i.e. number of nodes
      do i = 1, size(merge_grain_1)
        grain_1 = merge_grain_1(i)
        grain_2 = merge_grain_2(i)
        if(grain_1 == -1) then
          if(grain_2 == -1) cycle
          grain_select = grain_2
        else if(merge_grain_2(i) == -1) then
          grain_select = grain_1
        else if(grain_nnodes(grain_1) >= grain_nnodes(grain_2)) then
          grain_select = grain_1
        else
          grain_select = grain_2
        end if
        eta_remap_data%grains(remap_idx)%node = core_nodes(grain_select)
        eta_remap_data%grains(remap_idx)%separator = merge_seps(i)
        remap_idx = remap_idx + 1
      end do
      ! Check if the remapping selection is valid
      if(any(eta_remap_data%grains(1:n_remap_merge)%node <= 0)) then
        call img_printhead()
        call newline()
        write(output_unit, '(a)') 'Merge grain 1:   '
        do i = 1, size(merge_grain_1)
          write (output_unit, '(i8)', advance = 'no') merge_grain_1(i)
          write(output_unit, '(a)', advance = 'no') ' '
        end do
        call newline()
        write(output_unit, '(a)') 'Merge grain 2:   '
        do i = 1, size(merge_grain_2)
          write (output_unit, '(i8)', advance = 'no') merge_grain_2(i)
          write(output_unit, '(a)', advance = 'no') ' '
        end do
        call newline()
        write(output_unit, '(a)') 'Merge nodes 1:'
        do i = 1, size(merge_grain_1)
          if(merge_grain_1(i) > 0 .and. merge_grain_1(i) <= size(core_nodes)) then
            write (output_unit, '(i8)', advance = 'no') core_nodes(merge_grain_1(i))
          else
            write(output_unit, '(a)', advance = 'no') 'INVALID '
          end if
          write(output_unit, '(a)', advance = 'no') ' '
        end do
        call newline()
        write(output_unit, '(a)') 'Merge nodes 2:'
        do i = 1, size(merge_grain_2)
          if(merge_grain_2(i) > 0 .and. merge_grain_2(i) <= size(core_nodes)) then
            write (output_unit, '(i8)', advance = 'no') core_nodes(merge_grain_2(i))
          else
            write(output_unit, '(a)', advance = 'no') 'INVALID '
          end if
          write(output_unit, '(a)', advance = 'no') ' '
        end do
        call newline()
        write(output_unit, '(a)') 'Remap nodes:  '
        do i = 1, n_remap
          write (output_unit, '(i8)', advance = 'no') eta_remap_data%grains(i)%node
          write(output_unit, '(a)', advance = 'no') ' '
        end do
        write(output_unit, '(a,i0,a)') ' (last ', n_remap - n_remap_merge, ' for performance)'
        call warning_msg(proc_name, 'Faulty choice of remap nodes, skipping remapping...')
        n_remap = 0
        eta_remap_data = remap_data(n_remap, num_images())
        performance_remap = .false.
      end if
    else if(nlabels > ngrains) then
      call warning_msg(proc_name, 'Number of labels increased for lower separator, skipping remapping...')
      n_remap = 0
      eta_remap_data = remap_data(n_remap, num_images())
      performance_remap = .false.
    else
      if(performance_remap) then
        n_remap = ngrains
      else
        n_remap = 0
      end if
      eta_remap_data = remap_data(n_remap, num_images())
      eta_remap_data%grains%merge_risk = .false.
      eta_remap_data%grains%separator = settings%separators%inter
    end if
    ! If performance remapping is enabled, add all grains with the default
    ! inter-grain separator
    if(performance_remap) then
      do concurrent(i = 1:ngrains)
        eta_remap_data%grains(i)%node = core_nodes(i)
      end do
    end if
    if(verbosity_level >= 1) then
      call img_printhead()
      write(output_unit, '(a)', advance = 'no') 'Waiting for sync...'
      call newline()
      sync all
      if(this_image() == 1) write(output_unit, '(a)') 'get_remap_grains results:'
      call start_critical_in_order()
      call millisleep(5)
      call img_printhead()
      write(output_unit, '(2(i0,a))') label_calls, ' calls to label_areas, ', n_remap, ' remap nodes found.'
      call end_critical_in_order()
      sync all
    end if
  end subroutine get_remap_grains

  ! Due to a bug with allocatable coarrays in the Intel Fortran Compiler,
  ! this subroutine uses coarrays as input/output arguments that should
  ! be allocated inside the subroutine instead (since they are unused outside
  ! the subroutine), particularly 'eta_filled' and 'eta_remap_data'.
  ! The bug is reported on the Intel forums here:
  ! https://community.intel.com/t5/Intel-Fortran-Compiler/coarray-coarray-this-image-for-coarrays-allocated-inside-teams/m-p/1367474#M160524
  subroutine eta_reassign(grid, eta, eta_name, settings, eta_filled)
    type(regular_grid), intent(in) :: grid
    real(wp), intent(inout) :: eta(:)[*]
    character(*), intent(in) :: eta_name
    type(grain_reindex_settings), intent(in) :: settings
    logical, contiguous, intent(inout) :: eta_filled(:)[*]
    integer :: i, j, eta_target, node, stop_idx, n_remap_nodes, n_overlap,&
      my_grain_idx, ext_grain_idx, img, n_total_conflicts, eta_idx_perf, max_eta_idx
    logical, dimension(grid%n_nodes_total) :: eta_reassign_mask(grid%n_nodes_total)
    logical :: remap_found, fill_changed, find_candidates, candidates_changed
    integer, allocatable, dimension(:) :: overlap_eta_idx,&
      overlap_ext_grains, overlap_my_grains, node_coordinates
    logical, allocatable, dimension(:) :: grain_targets
    type(int_vector), allocatable, dimension(:) :: remap_overlaps
    type(remap_data), allocatable, dimension(:) :: eta_remap_data
    real(wp), allocatable, dimension(:) :: eta_remap
    real(wp) :: last_sep
    character(*), parameter :: proc_name = 'eta_reassign'
    !dir$ vector aligned

    allocate(overlap_eta_idx(0))
    allocate(overlap_ext_grains(0))
    allocate(overlap_my_grains(0))
    allocate(grain_targets(num_images()))
    allocate(eta_remap_data(num_images()))
    allocate(eta_remap(grid%n_nodes_total))
    grain_targets = .false.
    find_candidates = .false.
    last_sep = 0.0_wp

    max_eta_idx = get_highest_n_eta_active(eta)
    eta_idx_perf = max(max_eta_idx - settings%for_performance%n_eta, settings%for_performance%n_eta_min_active)
    call get_remap_grains(grid, eta, eta_name, settings, eta_remap_data(this_image()), this_image() >= eta_idx_perf)
    n_remap_nodes = eta_remap_data(this_image())%get_n()
    if(any(eta_remap_data(this_image())%grains%node > 0 .and. .not. eta_remap_data(this_image())%grains%merge_risk)) then
      ! if performance remapping is active, all grains could potentially be remapped
      eta_filled = .false.
    else
      ! else, use the default intergrain separator for the fill map
      eta_filled = eta > settings%separators%fill
    end if
    if(verbosity_level >= 1) then
      call output_image(eta, 'remap_' // eta_name // '_before', grid%axes(1)%n_nodes)
    end if
    allocate(remap_overlaps(n_remap_nodes))
    ! Get and store node coordinates scheduled for remapping
    do i = 1,n_remap_nodes
      node = eta_remap_data(this_image())%grains(i)%node
      node_coordinates = grid%linear_to_dim_idx(node)
      eta_reassign_mask = .false.
      call flood_fill_4(grid, eta > eta_remap_data(this_image())%grains(i)%separator, node_coordinates, eta_reassign_mask)
      eta_remap_data(this_image())%grains(i)%sections = mask_to_sections(eta_reassign_mask)
    end do
    do i = 1, n_remap_nodes - 1
      do j = i + 1, n_remap_nodes
        if(sections_overlap(eta_remap_data(this_image())%grains(i)%sections, &
          eta_remap_data(this_image())%grains(j)%sections)) then
          call output_image(eta_remap_data(this_image())%get_grain_map(grid%n_nodes_total), &
            'remap_' // eta_name // '_overlap_error', grid%axes(1)%n_nodes)
          call error_msg(proc_name, 'Grain areas for remapping overlap, check separators!')
        end if
      end do
    end do
    ! Coarray derived type is very buggy in gfortran (and also in ifort), so just keep local copies
    if(verbosity_level >= 1) then
      call output_image(eta_filled, 'remap_' // eta_name // '_filled_ideal', grid%axes(1)%n_nodes)
      call output_image(eta_remap_data(this_image())%get_grain_map(grid%n_nodes_total), &
        'remap_' // eta_name // '_schedule_ideal', grid%axes(1)%n_nodes)
    end if


    call broadcast_remap_data(eta_remap_data)
    ! Exclude grains that are impossible to remap
    do
      fill_changed = .false.
      do i = 1,n_remap_nodes
        if(eta_remap_data(this_image())%grains(i)%node == -1) cycle
        remap_found = .false.
        node = eta_remap_data(this_image())%grains(i)%node
        if(eta_remap_data(this_image())%grains(i)%merge_risk) then
          stop_idx = num_images()
        else
          stop_idx = this_image() - 1
        end if
        do eta_target = 1, stop_idx
          if(eta_target == this_image()) cycle
          if(verbosity_level >= 1) then
            call img_printhead()
            write(output_unit, '(2(a,i0))') 'Checking if remap possible for grain ', i, ' to eta ', eta_target
          end if
          if(.not. section_mask_overlap(eta_remap_data(this_image())%grains(i)%sections, eta_filled(:)[eta_target])) then
            remap_found = .true.
            if(find_candidates) then
              eta_remap_data(this_image())%grains(i)%targets(eta_target)%eligible = .true.
            else
              exit
            end if
          end if
        end do
        if(.not. remap_found) then
          if(find_candidates) call error_msg(proc_name, 'Found no remap destination in candidate checking phase!')
          if(verbosity_level >= 1) then
            call img_printhead()
            write(output_unit, '(a,i0,a)') 'Grain ', i, ' scheduled for remap has no possible destination!'
          end if
          fill_changed = .true.
          exit
        end if
        if(find_candidates .and. verbosity_level >= 1) then
          call img_printhead()
          write(output_unit, '(a)', advance = 'no') 'Found remap candidate targets: '
          print *, eta_remap_data(this_image())%grains(i)%targets%eligible
        end if
      end do
      !sync all
      if(fill_changed) then
        eta_remap_data(this_image())%grains(i)%node = -1
        eta_filled = eta_filled .or. sections_to_mask(eta_remap_data(this_image())%grains(i)%sections, grid%n_nodes_total)
      end if
      if(find_candidates) exit
      call co_any(fill_changed)
      find_candidates = .not. fill_changed
    end do
    call broadcast_remap_data(eta_remap_data)
    if(verbosity_level >= 1) then
      call output_image(eta_filled, 'remap_' // eta_name // '_filled_possible', grid%axes(1)%n_nodes)
      call output_image(eta_remap_data(this_image())%get_grain_map(grid%n_nodes_total), &
        'remap_' // eta_name // '_schedule_possible', grid%axes(1)%n_nodes)
    end if
    ! Find overlaps between remaining grains scheduled for remapping
    do img = 1, num_images()
      if(img == this_image()) cycle
      do my_grain_idx = 1, size(eta_remap_data(this_image())%grains)
        do ext_grain_idx = 1, size(eta_remap_data(img)%grains)
          associate(my_grain => eta_remap_data(this_image())%grains(my_grain_idx), &
              ext_grain => eta_remap_data(img)%grains(ext_grain_idx))
            if(my_grain%node <= 0 .or. ext_grain%node <= 0) cycle
            if(sections_overlap(my_grain%sections, ext_grain%sections)) then
              if(.not. any(overlap_my_grains == my_grain_idx&
                .and. overlap_ext_grains == ext_grain_idx&
                .and. overlap_eta_idx == img)) then

                overlap_my_grains = [overlap_my_grains, my_grain_idx]
                overlap_ext_grains = [overlap_ext_grains, ext_grain_idx]
                overlap_eta_idx = [overlap_eta_idx, img]
              end if
            end if
          end associate
        end do
      end do
    end do
    n_overlap = size(overlap_eta_idx)
    if(verbosity_level >= 1) then
      call start_critical_in_order()
      call millisleep(5)
      call img_printhead()
      write(output_unit, '(a)', advance = 'no') 'Found overlaps with other grains scheduled for remap: '
      call newline()
      do i = 1, n_overlap
        write(output_unit, '(3(a,i0))', advance = 'no') &
          'Local grain ', overlap_my_grains(i), ' with grain ', overlap_ext_grains(i), ' of img ', overlap_eta_idx(i)
      end do
      call end_critical_in_order()
    end if
    ! Lock in targets that result in no conflicts
    sync all
    do
      candidates_changed = .false.
      do concurrent(my_grain_idx = 1:n_remap_nodes,&
          eta_remap_data(this_image())%grains(my_grain_idx)%node > 0 .and. any(overlap_my_grains == my_grain_idx))
        grain_targets = eta_remap_data(this_image())%grains(my_grain_idx)%targets%eligible
        do concurrent(i = 1:n_overlap, overlap_my_grains(i) == my_grain_idx)
          ext_grain_idx = overlap_ext_grains(i)
          img = overlap_eta_idx(i)
          grain_targets = grain_targets .and. .not. eta_remap_data(img)%grains(ext_grain_idx)%targets%eligible
        end do
        if(any(grain_targets) .and. &
          any(eta_remap_data(this_image())%grains(my_grain_idx)%targets%eligible .neqv. grain_targets)) then
          candidates_changed = .true.
          eta_remap_data(this_image())%grains(my_grain_idx)%targets%eligible = grain_targets
        end if
      end do
      call co_any(candidates_changed)
      if(.not. candidates_changed) exit
      call broadcast_remap_data(eta_remap_data)
    end do
    ! Select targets so as to keep the amount of conflicts at a minimum
    do
      candidates_changed = .false.
      ! Calculate target badness, i.e. relative reduction of options
      do concurrent(i = 1:size(eta_remap_data(this_image())%grains))
        eta_remap_data(this_image())%grains(i)%targets%badness = 0.0_wp
      end do
      do i = 1, n_overlap
        my_grain_idx = overlap_my_grains(i)
        ext_grain_idx = overlap_ext_grains(i)
        img = overlap_eta_idx(i)
        where(eta_remap_data(this_image())%grains(my_grain_idx)%targets%eligible &
            .and. eta_remap_data(img)%grains(ext_grain_idx)%targets%eligible)
          eta_remap_data(this_image())%grains(my_grain_idx)%targets%badness = &
            eta_remap_data(this_image())%grains(my_grain_idx)%targets%badness &
            + 1.0_wp/count(eta_remap_data(img)%grains(ext_grain_idx)%targets%eligible)
        end where
      end do
      call broadcast_remap_data(eta_remap_data)
      ! Try selecting eta targets with minimum badness
      call eta_remap_data(this_image())%grains%select_target()
      ! Calculate the amount of resulting conflicts
      do i = 1, n_overlap
        my_grain_idx = overlap_my_grains(i)
        ext_grain_idx = overlap_ext_grains(i)
        img = overlap_eta_idx(i)
        if(eta_remap_data(this_image())%grains(my_grain_idx)%target_idx ==&
          eta_remap_data(img)%grains(ext_grain_idx)%target_idx) then
          eta_remap_data(this_image())%grains(my_grain_idx)%n_conflicts =&
            eta_remap_data(this_image())%grains(my_grain_idx)%n_conflicts + 1
        end if
      end do
      ! Lock in non-conflicting options
      do i = 1, eta_remap_data(this_image())%get_n()
        if(eta_remap_data(this_image())%grains(i)%n_conflicts == 0&
          .and. eta_remap_data(this_image())%grains(i)%target_idx > 0&
          .and. count(eta_remap_data(this_image())%grains(i)%targets%eligible) > 1) then

          candidates_changed = .true.
          eta_remap_data(this_image())%grains(i)%targets%eligible = .false.
          eta_remap_data(this_image())%grains(i)%targets(eta_remap_data(this_image())%grains(i)%target_idx)%eligible = .true.
        end if
      end do
      n_total_conflicts = sum(eta_remap_data(this_image())%grains%n_conflicts, eta_remap_data(this_image())%grains%node > 0)
      call co_sum(n_total_conflicts)
      call co_any(candidates_changed)
      if(.not. candidates_changed .or. n_total_conflicts == 0) exit
    end do
    ! If there are any conflicts left, resolve them based on simple priority rules:
    ! - Grains with merging risk before performance-oriented remapping
    ! - Grains with higher amount of conflicts are lower priority
    ! - Grains from higher index eta fields before lower index eta fields
    if(n_total_conflicts > 0) then
      do i = 1, n_overlap
        my_grain_idx = overlap_my_grains(i)
        if(eta_remap_data(this_image())%grains(my_grain_idx)%n_conflicts == 0) cycle
        ext_grain_idx = overlap_ext_grains(i)
        img = overlap_eta_idx(i)
        if(eta_remap_data(this_image())%grains(my_grain_idx)%target_idx ==&
          eta_remap_data(img)%grains(ext_grain_idx)%target_idx) then
          if(.not. eta_remap_data(this_image())%grains(my_grain_idx)%merge_risk .and.&
            eta_remap_data(img)%grains(ext_grain_idx)%merge_risk .or.&
            (eta_remap_data(this_image())%grains(my_grain_idx)%merge_risk .eqv.&
            eta_remap_data(img)%grains(ext_grain_idx)%merge_risk&
            .and. this_image() < img)) then

            print *, this_image(), 'setting grain', my_grain_idx, 'target to 0'
            eta_remap_data(this_image())%grains(my_grain_idx)%target_idx = 0
            eta_remap_data(this_image())%grains(my_grain_idx)%node = -1
            eta_filled = eta_filled .or. &
              sections_to_mask(eta_remap_data(this_image())%grains(my_grain_idx)%sections, grid%n_nodes_total)
            !do concurrent(j = 1:size(remap_coords(my_grain_idx)%v))
            !eta_filled(remap_coords(my_grain_idx)%v(j)) = .true.
            !eta_remap_data(this_image())%grain_map(remap_coords(my_grain_idx)%v(j)) = 0
            !end do
          end if
        end if
      end do
      eta_remap_data(this_image())%grains%n_conflicts = 0
      call broadcast_remap_data(eta_remap_data)
      do i = 1, n_overlap
        my_grain_idx = overlap_my_grains(i)
        ext_grain_idx = overlap_ext_grains(i)
        img = overlap_eta_idx(i)
        if(eta_remap_data(this_image())%grains(my_grain_idx)%target_idx > 0&
          .and. eta_remap_data(this_image())%grains(my_grain_idx)%target_idx ==&
          eta_remap_data(img)%grains(ext_grain_idx)%target_idx) then
          eta_remap_data(this_image())%grains(my_grain_idx)%n_conflicts =&
            eta_remap_data(this_image())%grains(my_grain_idx)%n_conflicts + 1
        end if
      end do
      n_total_conflicts = sum(eta_remap_data(this_image())%grains%n_conflicts)
      call co_sum(n_total_conflicts)
      if(n_total_conflicts > 0) then
        call start_critical_in_order()
        call millisleep(5)
        call img_printhead()
        write(output_unit, '(a)', advance = 'no') 'Got the following remapping targets: '
        call newline()
        do i = 1, eta_remap_data(this_image())%get_n()
          write(output_unit, '(3(a,i0))') 'remap node ', i, ', selected: ', &
            eta_remap_data(this_image())%grains(i)%target_idx, &
            ', conflicts: ', eta_remap_data(this_image())%grains(i)%n_conflicts
        end do
        call end_critical_in_order()
        sync all
        call millisleep(100)
        call error_msg(proc_name, 'could not resolve remapping conflicts!')
      end if
    end if
    ! Check if remapping targets are still available, especially if space made
    ! available by remapping of other grains is still unoccupied
    sync all
    do
      candidates_changed = .false.
      do i = 1, n_remap_nodes
        img = eta_remap_data(this_image())%grains(i)%target_idx
        if(img <= 0) cycle
        !if(any_idx(eta_filled(:)[img], remap_coords(i)%v)) then
        if(section_mask_overlap(eta_remap_data(this_image())%grains(i)%sections, eta_filled(:)[img])) then
          call img_printhead()
          write(output_unit, '(2(a,i0),a)') 'Target ', img, ' no longer available for grain ', i, ', setting to 0'
          eta_remap_data(this_image())%grains(i)%target_idx = 0
          candidates_changed = .true.
          exit
        end if
      end do
      sync all
      if(candidates_changed) then
        eta_filled = eta_filled .or. sections_to_mask(eta_remap_data(this_image())%grains(i)%sections, grid%n_nodes_total)
        !do concurrent(j = 1:size(remap_coords(i)%v))
        !    eta_filled(remap_coords(i)%v(j)) = .true.
        !eta_remap_data(this_image())%grain_map(remap_coords(i)%v(j)) = 0
        !end do
      end if
      call co_any(candidates_changed)
      if(.not. candidates_changed) exit
    end do
    sync all
    if(verbosity_level >= 1) then
      call output_image(eta_filled, 'remap_' // eta_name // '_filled_final', grid%axes(1)%n_nodes)
      call output_image(eta_remap_data(this_image())%get_grain_map(grid%n_nodes_total), &
        'remap_' // eta_name // '_schedule_final', grid%axes(1)%n_nodes)
    end if
    if(verbosity_level >= 1) then
      call start_critical_in_order()
      call millisleep(5)
      call img_printhead()
      write(output_unit, '(a)', advance = 'no') 'Got the following remapping targets: '
      call newline()
      do i = 1, size(eta_remap_data(this_image())%grains)
        if(eta_remap_data(this_image())%grains(i)%node <= 0) cycle
        write(output_unit, '(2(a,i0))') 'remap node ', i, ', selected: ', eta_remap_data(this_image())%grains(i)%target_idx
      end do
      call end_critical_in_order()
    end if

    eta_remap = merge(eta, 0.0_wp, eta_remap_data(this_image())%get_grain_map(grid%n_nodes_total) > 0)
    where(eta_remap_data(this_image())%get_grain_map(grid%n_nodes_total) > 0) 
      eta = settings%separators%inter
    end where
    do eta_target = 1, num_images()
      sync all
      if(eta_target == this_image()) then
        sync images(left_image())
      else
        if(this_image() /= next_image(eta_target)) sync images(left_image())
        call img_printhead()
        write(output_unit, '(a,i0)') 'Writing data to eta ', eta_target
        do i = 1, n_remap_nodes
          if(eta_remap_data(this_image())%grains(i)%target_idx /= eta_target) cycle
          eta_remap_data(this_image())%grains(i)%node = 0
          do j = 1, size(eta_remap_data(this_image())%grains(i)%sections)
            associate(first_idx => eta_remap_data(this_image())%grains(i)%sections(j)%first, &
                last_idx => eta_remap_data(this_image())%grains(i)%sections(j)%last)
              eta(first_idx:last_idx)[eta_target] = eta_remap(first_idx:last_idx)
            end associate
          end do
          !do concurrent(j = 1:size(remap_coords(i)%v))
          !node = remap_coords(i)%v(j)
          !eta(node)[eta_target] = eta_remap(node)
          !end do
        end do
        sync images(right_image())
      end if
    end do

    sync all
    if(verbosity_level >= 1) then
      call output_image(eta, 'remap_' // eta_name // '_after', grid%axes(1)%n_nodes)
      call start_critical_in_order()
      call millisleep(5)
      call img_printhead()
      write(output_unit, '(3(i0,a))') count(eta_remap_data(this_image())%grains%node /= 0), &
        ' out of ', size(eta_remap_data(this_image())%grains%node), ' left (', &
        count(eta_remap_data(this_image())%grains%node /= 0 .and. eta_remap_data(this_image())%grains%merge_risk), &
        ' with risk of merging)'
      call end_critical_in_order()
      sync all
      call millisleep(5)
      if(this_image() == 1) write(output_unit, '(a)') 'Grain reassignment finished!'
    end if
  end subroutine eta_reassign

  !FIXME: this implicitly assumes that this_image() == eta_idx
  integer function get_highest_n_eta_active(eta_vals) result(highest_n_eta_active)
    real(wp), intent(in) :: eta_vals(:)

    if(any(abs(eta_vals) > epsilon(1.0_wp))) then
      highest_n_eta_active = this_image()
    else
      highest_n_eta_active = 0
    end if
    call co_max(highest_n_eta_active)
  end function get_highest_n_eta_active

  pure function true_indices(a)
    logical, intent(in) :: a(:)
    integer :: true_indices(count(a))
    integer :: true_count, i

    true_count = 0
    true_indices = -1
    ! Parallelizing this loop leads to segmentation faults
    ! Not thread-safe as order is important
    !DEC$ NOPARALLEL
    do i = 1, size(a)
      if(a(i)) then
        true_count = true_count + 1
        true_indices(true_count) = i
      end if
    end do
    if(true_count /= count(a)) error stop &
      'ERROR: Not all true indices found in function true_indices!'
    if(any(true_indices < 0)) error stop &
      'ERROR: Not all true indices correctly assigned in function true_indices!'
  end function true_indices

  pure logical function any_idx(a, idx)
    logical, intent(in) :: a(:)
    integer, intent(in) :: idx(:)
    integer :: i

    any_idx = .false.
    do concurrent(i = 1:size(idx))
      if(a(idx(i))) any_idx = .true.
    end do
  end function any_idx

  pure subroutine find_surface(grid, phi_pore, surf_x_node)
    type(regular_grid), intent(in) :: grid
    real(wp), intent(in) :: phi_pore(:)
    real(wp), intent(out) :: surf_x_node(:)
    real(wp) :: phi_1, phi_2
    real(wp), parameter :: surf_phi = 0.5_wp
    integer, allocatable :: p(:)
    integer :: linear_idx, sub_node, x_idx
    type(regular_grid) :: sub_grid
    character(*), parameter :: proc_name = 'find_surface'

    surf_x_node = 0.0_wp
    allocate(p(grid%n_dims))
    p = 1
    call sub_grid%init(grid%axes(2:))
    do concurrent(sub_node = 1:size(surf_x_node))
      p(2:) = sub_grid%linear_to_dim_idx(sub_node)
      linear_idx = grid%dim_to_linear_idx(p)
      phi_1 = phi_pore(linear_idx)
      if(phi_1 < surf_phi) call error_msg(proc_name, 'First node in x-dimension is not pore phase')
      do x_idx = 2, grid%axes(1)%n_nodes
        phi_2 = phi_pore(linear_idx + x_idx - 1)
        if(phi_1 >= surf_phi .and. phi_2 < surf_phi) then
          surf_x_node(sub_node) = x_idx - 1 + (surf_phi - phi_1)/(phi_2 - phi_1)
          exit
        end if
        phi_1 = phi_2
        if(x_idx == grid%axes(1)%n_nodes) call error_msg(proc_name, 'Could not find surface')
      end do
    end do
  end subroutine find_surface

  ! TODO: this needs to be extended to n dimensions
  !pure function find_surface(grid, phi_pore) result(surf_x_node)
  !    type(regular_grid), intent(in) :: grid
  !    real(wp), intent(in) :: phi_pore(:)
  !    real(wp) :: surf_x_node(product(grid%axes(2:)%n_nodes))
  !    real(wp) :: phi_1, phi_2
  !    integer :: node, node_x, node_y
  !    real(wp), parameter :: surf_phi = 0.5_wp
  !    
  !    surf_x_node = 0.0_wp
  !    do concurrent(node_y = 1:grid%axes(2)%n_nodes)
  !        phi_1 = phi_pore(get_node(1, node_y))
  !        do node_x = 2, grid%axes(1)%n_nodes
  !            node = get_node(node_x, node_y)
  !            phi_2 = phi_pore(node)
  !            if(phi_1 >= surf_phi .and. phi_2 < surf_phi) then
  !                surf_x_node(node_y) = node_x - 1 + (surf_phi - phi_1)/(phi_2 - phi_1)
  !                exit
  !            end if
  !            phi_1 = phi_2
  !        end do
  !    end do
  !end function find_surface

  !pure subroutine fourier_clean_var_lb(c_vals, surf_x_node)
  !    real(wp), intent(inout) :: c_vals(:)
  !    real(wp), contiguous, intent(in) :: surf_x_node(:)
  !    integer :: fc_xstart, fc_xend, node_y
  !    
  !    do concurrent(node_y = 1:ny)
  !        fc_xend = (node_y - 1)*nx + nint(surf_x_node(node_y)) - 1
  !        fc_xstart = fc_xend - fc_nodes + 1
  !        c_vals(fc_xstart:fc_xend) = 0.0_wp
  !    end do
  !end subroutine fourier_clean_var_lb

  !pure subroutine fourier_balance_var_lb(c_field, phi_field, surf_x_node)
  !    type(fft_field), intent(inout) :: c_field
  !    type(fft_field), intent(in) :: phi_field(:)
  !    real(wp), contiguous, intent(in) :: surf_x_node(:)
  !    integer :: fc_xstart, fc_xend
  !    
  !    do concurrent(i = 1:ny)
  !        fc_xend = (i - 1)*nx + nint(surf_x_node(i)) - 1
  !        fc_xstart = fc_xend - fc_nodes + 1
  !        !FIXME: Maybe setting to 0 is sufficient? What is the influence on eta?
  !        c_field%vals(fc_xstart:fc_xend) = mix_kks_local(c_field, phi_field, fc_xstart, fc_xend)
  !    end do
  !end subroutine fourier_balance_var_lb

  !TODO: reimplement non-periodic boundary conditions
  !subroutine set_boundaries(c_field, phi_field, fc_data, c_fft, fft, n_c_fft, surf_x_node)
  !    type(fft_field), intent(inout) :: c_field
  !    type(fft_field), intent(inout) :: phi_field(:)
  !    type(fourier_cont_data), intent(in) :: fc_data
  !    complex(wp), dimension(grid%n_nodes_total_fft), intent(out) :: c_fft
  !    type(fft_field), intent(in) :: fft
  !    integer(int64), intent(inout) :: n_c_fft
  !    integer :: phi_center_idx
  !    real(wp), intent(out) :: surf_x_node(ny)
  !    integer, parameter :: phi_nodes = fc_nodes + 5*nint(if_nodes)
  !    real(wp) :: smooth_x(phi_nodes), phi_offset
  !    
  !    !if(.not. c_field%bc%is_periodic) then
  !    !    surf_x_node = find_surface(phi_field(1)%vals)
  !    !    if(.not. lazy_boundary) then
  !    !        call fourier_cont_var_lb(c_field, fc_data, surf_x_node)
  !    !        call fft%forward_alldim(c_field%vals, c_fft)
  !    !        n_c_fft = n_c_fft + 1
  !    !    end if
  !    !    ! FIXME: Can this be improved?
  !    !    ! Phi seems to cause oscillations, broadening the interface zone helps somewhat.
  !    !    ! There should be a better function to handle this.
  !    !    !do concurrent(i = 1:phi_nodes)
  !    !    !    smooth_x(i) = 2*real(i - 1, wp)/(fc_nodes - 1) - 1
  !    !    !end do
  !    !    do concurrent(i = 1:ny)
  !    !        phi_field(1)%vals((i - 1)*nx + 1:(i - 1)*nx + buffer_nodes_l + phi_nodes) = 0.0_wp
  !    !        !phi_offset = 2*(surf_x_node(i) - nint(surf_x_node(i)))/(fc_nodes - 1)
  !    !        !!surf_x_nearest_node = nint(surf_x_node(i))
  !    !        !!bc_offset = surf_x_node(i) - surf_x_nearest_node
  !    !        !!fc_xend = (i - 1)*nx + surf_x_nearest_node - 1
  !    !        !!fc_xstart = fc_xend - fc_nodes + 1
  !    !        !phi_field(1)%vals((i - 1)*nx + nint(surf_x_node(i)) - fc_nodes:(i - 1)*nx + nint(surf_x_node(i)) - fc_nodes + phi_nodes - 1) =&
  !    !        !    smoothing_func(smooth_x - phi_offset, fc_w_phi)
  !    !        !!phi_field(1)%vals((i - 1)*nx + 1:(i - 1)*nx + nint(surf_x_node(i)) + phi_center) =&
  !    !        !!    equi_phi(x(1:nint(surf_x_node(i)) + phi_center) - surf_x_node(i) + phi_center, 2*if_nodes)
  !    !    end do
  !    !end if
  !end subroutine set_boundaries

  subroutine c_average(c_field, species_idx, c_mean, x_mean)
    type(fft_field), intent(in) :: c_field(:)
    integer, intent(in) :: species_idx
    !real(wp), intent(in) :: surf_x_node(:)
    real(wp), dimension(:), intent(out) :: x_mean, c_mean
    real(wp), dimension(size(c_field)) :: n
    real(wp) :: n_total
    ! integer :: fc_xend, fc_xstart
    integer :: i

    if(species_idx > 0) then
      !FIXME: restore fourier continuation feature
      !if(c_field(species_idx)%real%grid%axes%bound_cond%is_periodic) then
      n(species_idx) = sum(c_field(species_idx)%real%vals)
      !else
      !n(species_idx) = 0.0_wp
      !do concurrent(i = 1:ny)
      !    fc_xend = (i - 1)*nx + nint(surf_x_node(i)) - 1
      !    fc_xstart = fc_xend - fc_nodes + 1
      !    n(species_idx) = n(species_idx) + sum(c_field(species_idx)%real%vals((i - 1)*nx + 1:fc_xstart - 1)) + sum(c_field(species_idx)%real%vals(fc_xend + 1:i*nx))
      !end do
      !end if
      n_total = n(species_idx)
      call co_sum(n_total)
      x_mean(species_idx) = n(species_idx)/n_total
      c_mean(species_idx) = n(species_idx)/c_field(species_idx)%real%grid%n_nodes_total
    else
      n_total = 0
      call co_sum(n_total)
    end if
    do i = 1, size(c_mean)
      call co_broadcast(x_mean(i), i)
      call co_broadcast(c_mean(i), i)
    end do
  end subroutine c_average

  subroutine calc_mu_c(c_field, phases, species_idx, mu_c, use_bdfab, bdfab_stage)
    type(fft_field), intent(inout) :: c_field
    type(fft_field), intent(inout) :: mu_c
    type(phase_set), intent(inout) :: phases
    logical, intent(in) :: use_bdfab
    integer, intent(in) :: species_idx, bdfab_stage
    integer :: i
    character(*), parameter :: proc_name = 'calc_mu_c'

    do i = 1, phases%n_phases - 1
      call c_field%real%check_time_is_equal(phases%phi(i), proc_name)
    end do
    if(mu_c%fft%time_is_equal(c_field%real)) return
    if(use_bdfab .and. bdfab_stage > 2) call mu_c%swap_fft()
    call phases%calc_c_mix(species_idx)
    call phases%c_mix(species_idx)%calc_fft_laplacian()
    if(phases%equal_p_coeff(species_idx)) then
      call mu_c%fft_allocate_check()
      mu_c%fft%vals = -phases%c_mix(species_idx)%fft_laplacian%vals
      mu_c%fft%time = phases%c_mix(species_idx)%fft_laplacian%time
    else
      call phases%calc_p_mix(species_idx)
      mu_c%real%vals = (c_field%real%vals - phases%c_mix(species_idx)%real%vals)/phases%p_mix(species_idx)%real%vals
      mu_c%real%time = c_field%real%time
      call mu_c%calc_real_laplacian()
      call phases%p_mix(species_idx)%calc_real_laplacian()
      mu_c%real%vals = phases%p_mix(species_idx)%real%vals*mu_c%real_laplacian%vals &
        - mu_c%real%vals*phases%p_mix(species_idx)%real_laplacian%vals
      call mu_c%fft%invalidate()
      call mu_c%fft_laplacian%invalidate()
      call mu_c%real_laplacian%invalidate()
      call mu_c%forward_alldim()
      mu_c%fft%vals = 0.5_wp*(mu_c%fft%vals - phases%c_mix(species_idx)%fft_laplacian%vals)
!      call phases%calc_log_p_mix(species_idx)
!      mu_c%real%vals = c_field%real%vals - phases%c_mix(species_idx)%real%vals
!      mu_c%real%time = c_field%real%time
!      call phases%log_p_mix(species_idx)%calc_real_laplacian()
!      call phases%log_p_mix(species_idx)%calc_grad()
!      call mu_c%calc_grad()
!      do concurrent(i = 1:mu_c%real%grid%n_nodes_total)
!        mu_c%real%vals(i) = -(mu_c%real%vals(i)*phases%log_p_mix(species_idx)%real_laplacian%vals(i)&
!          + fft_nabla_product(mu_c, phases%log_p_mix(species_idx), i))
!      end do
!      do dim = 1, mu_c%real%grid%n_dims
!        call mu_c%real_grad(dim)%invalidate()
!      end do
!      call mu_c%forward_alldim()
!      do concurrent(node = 1:mu_c%fft%grid%n_nodes_total)
!        mu_c%fft%vals(node) = mu_c%fft%vals(node) - phases%c_mix(species_idx)%fft_laplacian%vals(node)
!      end do
    end if
  end subroutine calc_mu_c

  subroutine spectral_step_c(c_field, mu_c, D, dt, use_bdfab, bdfab_stage, equal_p_coeff)
    type(fft_field), intent(inout) :: c_field
    type(fft_field), intent(inout) :: mu_c
    real(wp), intent(in) :: D, dt
    logical, intent(in) :: use_bdfab, equal_p_coeff
    integer, intent(in) :: bdfab_stage
    character(*), parameter :: proc_name = 'spectral_step_c'
    !dir$ vector aligned

    !TODO: include variable D by expanding nabla(D*nabla(c)) to D*laplacian(c) + nabla(D)*nabla(c)
    if(verbosity_level >= 5) call print_message('Entering spectral_step_c...')
    if(equal_p_coeff) then
      call spectral_solve(c_field, mu_c, D, dt, use_bdfab, bdfab_stage)
    else
      call spectral_solve(c_field, mu_c, D, dt, use_bdfab, bdfab_stage, a_coeff_opt = 0.5_wp)
    end if
    ! call spectral_solve(c_field, mu_c, D, dt, use_bdfab, bdfab_stage)
    call c_field%backward_alldim()
  end subroutine spectral_step_c

  subroutine spectral_step_eta(order_parameters, df_deta, mu_eta, M_eta, dt, use_bdfab, bdfab_stage)
    type(order_parameter_set), intent(inout) :: order_parameters
    type(phase_field_energies), intent(in) :: df_deta
    type(fft_field), intent(inout) :: mu_eta
    real(wp), intent(in) :: M_eta, dt
    logical, intent(in) :: use_bdfab
    integer, intent(in) :: bdfab_stage
    character(*), parameter :: proc_name = 'spectral_step_eta'
    !dir$ vector aligned

    if(verbosity_level >= 5) call print_message('Entering spectral_step_eta...')
    if(use_bdfab .and. bdfab_stage > 2) call mu_eta%swap_fft()
    call df_deta%chem%check_time_is_equal(df_deta%order, proc_name)
    mu_eta%real%vals = -(df_deta%chem%vals + df_deta%order%vals)
    mu_eta%real%time = df_deta%chem%time
    call mu_eta%forward_alldim()
    call spectral_solve(order_parameters%eta, mu_eta, M_eta, dt, use_bdfab, bdfab_stage,&
      a_coeff_opt = order_parameters%eta%interface_data%get_kappa())
    call order_parameters%eta%backward_alldim()
  end subroutine spectral_step_eta

  subroutine real_broadcast(var, var_name, step, img_opt)
    real(wp), intent(inout) :: var(:)
    character(*), intent(in) :: var_name
    integer(int64), intent(in) :: step
    integer, intent(in), optional :: img_opt
    integer :: img

    if(present(img_opt)) then
      img = img_opt
    else
      img = 1
    end if
    if(verbosity_level >= 5) call print_message('Waiting for ' // var_name // ' sync...', step)
    call co_broadcast(var, img)
    if(verbosity_level >= 1 .and. verbosity_level >= 5) call print_message(var_name // ' synced!', step, 1)
  end subroutine real_broadcast

  !pure subroutine prediffuse_interface(a, surf_x_node, erf_target, n_diff_nodes)
  !    type(fft_field), intent(inout) :: a
  !    real(wp), contiguous, intent(in) :: surf_x_node(:)
  !    real(wp), intent(in) :: erf_target
  !    integer, intent(in) :: n_diff_nodes
  !    real(wp) :: bc_offset, c_inner
  !    integer :: diff_xstart, diff_xend, surf_x_nearest_node
  !    
  !    do concurrent(i = 1:ny)
  !        surf_x_nearest_node = nint(surf_x_node(i))
  !        bc_offset = (surf_x_nearest_node - surf_x_node(i))*erf_target/(n_diff_nodes - 1)
  !        diff_xstart = (i - 1)*nx + surf_x_nearest_node
  !        diff_xend = diff_xstart + n_diff_nodes - 1
  !        c_inner = a%vals(diff_xend + 1)
  !        a%vals(diff_xstart:diff_xend) = erf(linspace(0.0_wp, erf_target, n_diff_nodes) + bc_offset)*c_inner
  !    end do
  !end subroutine prediffuse_interface

  real(wp) function get_step_factor(dt_val, species_properties, temperature, t_now, M_eta) result(step_factor)
    real(wp), intent(in) :: dt_val, M_eta, t_now
    type(species_properties_t), intent(in) :: species_properties(:)
    type(piecewise_linear), intent(in) :: temperature

    step_factor = dt_val*max(&
      maxval(species_properties%diffusion_coefficient%get_value(temperature%get_value_at(t_now)), &
      .not. species_properties%instant_diffusion), M_eta)
    call co_max(step_factor)
  end function get_step_factor

  real(wp) function get_total_step_factor(t_start, t_stop, species_properties, temperature, M_eta) result(total_step_factor)
    real(wp), intent(in) :: M_eta, t_start, t_stop
    type(species_properties_t), intent(in) :: species_properties(:)
    type(piecewise_linear), intent(in) :: temperature
    real(wp) :: t_now, dt_val
    real(wp), allocatable :: profile_times(:)
    logical :: isotherm
    integer, parameter :: n_partition = 1E6

    total_step_factor = 0.0_wp
    t_now = t_start
    allocate(profile_times(temperature%get_n_val()))
    profile_times = temperature%get_x_known()
    dt_val = (profile_times(temperature%get_low_idx(t_now) + 1) - t_now)/n_partition
    do
      isotherm = abs(temperature%get_slope_at(t_now)) < epsilon(1.0_wp)
      if(isotherm) dt_val = profile_times(temperature%get_low_idx(t_now) + 1) - t_now
      dt_val = min(dt_val, t_stop - t_now)
      total_step_factor = total_step_factor + get_step_factor(dt_val, species_properties, temperature, t_now, M_eta)
      t_now = t_now + dt_val
      if(t_now >= t_stop) exit
      if(isotherm) then
        dt_val = (profile_times(temperature%get_low_idx(t_now) + 1) - t_now)/n_partition
        isotherm = .false.
      end if
    end do
  end function get_total_step_factor
end module mod_pfm
