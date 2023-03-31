! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_mu

  use mod_base_functions
  use mod_globals, only: verbosity_level
  use liboctq, only: gtp_equilibrium_data, GSNOACS, zero, one, tqini, tqquiet, tqtgsw, tqrpfil, tqphsts, tqgpi, tqsetc, maxc, &
    enter_composition_set, tqgnp, ask_default_constitution, set_default_constitution, tqlc, tqlr, gx, tqce, tqgetv, tqgpi2, &
    get_phase_variance, get_sublattice_number, get_state_var_value, tqgphc1, tqcceq, tqselceq
  use mod_th, only: tqgdmat_th
  use mod_parameters, only: eta_team_num, opencalphad_settings
  use mod_species_properties
  use mod_phase_properties
  use mod_regular_grid
  use mod_time_dep_vars
  use mod_spectral
  use mod_phase_set
  use mod_order_parameter_set
  use mod_phase_field_energies
  use mod_constants
  use mod_unit_factors
  implicit none
  private

  public :: opencalphad_equilib_phase, get_c_init, calc_df_chem_dc, calc_df_chem_deta_common, calc_df_chem_deta, &
    calc_df_order_deta, calc_numeric_energy_derivatives, gtp_equilibrium_data, GSNOACS, zero, one, tqini, tqquiet, &
    tqtgsw, tqrpfil, tqphsts, tqsetc, maxc, enter_composition_set, tqgpi, ask_default_constitution, tqlc, tqlr, &
    opencalphad_setup, tqcceq, tqgnp, tqselceq

contains

  subroutine opencalphad_setup(oc_settings, species_properties, phase_properties, pressure, units, current_equi, component_names)
    type(opencalphad_settings), intent(in) :: oc_settings
    type(species_properties_t), intent(in) :: species_properties(:)
    type(phase_properties_t), intent(in) :: phase_properties(:)
    real(wp), intent(in) :: pressure
    type(unit_factor_t), intent(in) :: units
    type(gtp_equilibrium_data), pointer, intent(out) :: current_equi
    character(len=24), dimension(maxc), intent(out) :: component_names
    integer :: i, j, phase_idx, compset_idx, phase_tuple_idx, n_phase_tuples, n_phases_oc
    character(16) :: compset_phase_name, compset_suffix
    type(string), dimension(:), allocatable :: phase_name_parts

    ! Initialize the OpenCalphad application software interface (OCASI)
    if(verbosity_level >= 5) call print_message('OpenCalphad: tqini')
    call tqini(i, current_equi)
    ! Silence output
    if(verbosity_level <= 3) call tqquiet(.true.)
    ! Load the selected species_properties from the tdb file given in oc_settings%tdb_filename
    if(verbosity_level >= 5) call print_message('OpenCalphad: tqrpfil')
    call tqrpfil(oc_settings%tdb_filename, i, n_phase_tuples, species_properties%name, current_equi, component_names)
    ! Suspend all phases
    if(verbosity_level >= 5) call print_message('OpenCalphad: Selecting phases')
    call tqphsts(-1, -3, n_phase_tuples, zero, current_equi)
    ! Enter selected phases
    do i = 1, size(phase_properties)
      if(trim(phase_properties(i)%name) == 'pore') cycle
      if(phase_properties(i)%is_compset) then
        call split_string(phase_properties(i)%name, '_', phase_name_parts)
        compset_phase_name = join_string(phase_name_parts(1:size(phase_name_parts) - 1), '_')
        compset_suffix = trim(phase_name_parts(size(phase_name_parts))%c)
        call tqgpi2(phase_idx, compset_idx, compset_phase_name, current_equi)
        call enter_composition_set(phase_idx, '    ', compset_suffix, compset_idx)
      else
        call tqgpi2(phase_idx, compset_idx, phase_properties(i)%name, current_equi)
      end if
      if(len_trim(phase_properties(i)%default_comp) > 0) then
        ! variable 'last' in ask_default_constitution refers to the current index in the command line string
        j = 1
        call ask_default_constitution(phase_properties(i)%default_comp, j, phase_idx, compset_idx, current_equi)
        call set_default_constitution(phase_idx, compset_idx, current_equi)
      end if
      if(phase_idx <= 0) then
        error stop 'ERROR: Phase ' // trim(phase_properties(i)%name) // ' not found in selected database!'
      end if
    end do
    call tqgnp(n_phases_oc, n_phase_tuples, current_equi)
    do i = 1, size(phase_properties)
      call tqgpi(phase_tuple_idx, phase_properties(i)%name, current_equi)
      if(verbosity_level >= 5) call print_message('Opencalphad: Setting phase '&
        // phase_properties(i)%name // ' at phase tuple idx ' // convert_to_char(phase_tuple_idx) // ' to status entered.')
      call tqphsts(phase_tuple_idx, 0, n_phase_tuples, one, current_equi)
    end do
    call tqsetc('P ', 0, 0, pressure*units%pascal, component_names, current_equi)
    call tqsetc('N ', 0, 0, 1.0_wp, component_names, current_equi)
    if(verbosity_level >= 5) call print_message('OpenCalphad: initialized!')
  end subroutine opencalphad_setup

  ! FIXME: This probably only works if species names are in the same order as in calphad equilib data, i.e. in alphabetical order (?)
  subroutine opencalphad_equilib_phase(temp, pressure, x_base, phase_idx, units, constants, current_equi, component_names, &
      c_field, phases, gridmin_opt)
    real(wp), intent(in) :: temp, pressure, x_base(:)
    integer, intent(in) :: phase_idx
    type(gtp_equilibrium_data), pointer, intent(inout) :: current_equi
    type(fft_field), intent(inout) :: c_field(:)
    type(phase_set), intent(inout) :: phases
    real(wp), dimension(phases%n_species, phases%n_phases) :: mole_equi, conc_equi, delta_mu, omega
    real(wp), dimension(phases%n_phases) :: f0
    logical, intent(in), optional :: gridmin_opt
    type(unit_factor_t), intent(in) :: units
    type(constants_t), intent(in) :: constants
    real(wp) :: fetch_var(1), n_phase_total
    real(wp), dimension(phases%n_species) :: x_eq, cpot
    real(wp), dimension(:), allocatable :: mugrad
    real(wp), dimension(2) :: Tp
    character(len=24), dimension(maxc), intent(in) :: component_names
    character(len=24), dimension(maxc) :: consnames 
    character(len=128) :: encoded
    real(wp), allocatable, dimension(:,:) :: ysave
    real(wp), allocatable, dimension(:) :: nsites, mugrad_diag
    real(wp) :: yfrac(maxc)
    real(wp), dimension(phases%n_species) :: extra(2), mu_global_eq
    integer :: i, species_idx, phase_tuple_idx, oc_phase_idx, compset_idx, n_phase_tuples, n3, n1, nsub, spix(maxc),&
      gridmin_int, node, pore_idx
    integer, dimension(:), allocatable :: nkl
    logical :: phase_stoichiometric

    if(present(gridmin_opt)) then
      if(gridmin_opt) then
        gridmin_int = 0
      else
        gridmin_int = -1
      end if
    else
      gridmin_int = -1
    end if

    consnames = ''
    omega = -huge(0.0_wp)
    mole_equi = -huge(0.0_wp)
    delta_mu = -huge(0.0_wp)
    f0 = -huge(0.0_wp)
    cpot = 0
    n3 = 0
    Tp = [temp*units%kelvin, pressure*units%pascal]

    pore_idx = 0
    if(phase_idx > 0) then
      if(trim(phases%phi(phase_idx)%properties%name) == 'pore') pore_idx = phase_idx
    end if
    call co_max(pore_idx)

    ! Set conditions
    if(verbosity_level >= 5) call print_message('Setting T = ' // convert_to_char(temp))
    call tqsetc('T ', 0, 0, temp*units%kelvin, component_names, current_equi)
    do species_idx = 1, phases%n_species - 1
      if(verbosity_level >= 5) then
        call print_message('Setting X(' // phases%species_name(species_idx) // ') = ' // convert_to_char(x_base(species_idx)))
      end if
      call tqsetc('X ', species_idx, 0, x_base(species_idx), component_names, current_equi)
      if(gx%bmperr /= 0) error stop 'ERROR: opencalphad_equilib_phase: OC could not set conditions!'    
    enddo
    if(verbosity_level >= 5) call tqlc(output_unit, current_equi)
    ! Calculate global equilibrium
    call tqce(' ', gridmin_int, 0, n_phase_tuples, zero, current_equi, ysave)
    if(gx%bmperr /= 0) then
      call tqlr(output_unit, current_equi)
      error stop 'ERROR: opencalphad_equilib_phase: OC calculation error!'    
    end if
    ! Print results
    if(verbosity_level >= 5) call tqlr(output_unit, current_equi)
    n3 = phases%n_species
    call tqgetv('mu', -1, 0, n3, mu_global_eq, current_equi, component_names)
    if(gx%bmperr /= 0) error stop
    ! Pore phase
    if(pore_idx > 0) then
      mole_equi(:,pore_idx) = 0.0_wp
      delta_mu(:,pore_idx) = 0.0_wp
      omega(:,pore_idx) = sqrt(huge(0.0_wp))
      f0(pore_idx) = 0.0_wp
    end if
    ! Calculate darken matrix
    if(phase_idx > 0 .and. this_image() <= phases%grain_img_idx) then
      if(phases%phi(phase_idx)%properties%calphad_use_phase_conc) then
        phases%phi(phase_idx)%amount = 0
        phases%phi(phase_idx)%x_mean = 0
        ! FIXME: there should be a better way to check if fields have been initialized
        !if(any(c_field%vals(1) > 0) .or. any(phi_field(1)%vals > 0)) then
        if(any(phases%phi(phase_idx)%p_coeff > 0)) then
          do i = 1, phases%n_species
            call phases%calc_c_mix(i)
            call phases%calc_p_mix(i)
          end do
          do concurrent(node = 1:phases%grid%n_nodes_total)
            phases%phi(phase_idx)%amount = phases%phi(phase_idx)%amount + phases%get_phi(phase_idx, node)
          end do
          if(phases%phi(phase_idx)%amount > delta_num) then
            do concurrent(i = 1:phases%n_species, node = 1:phases%grid%n_nodes_total)
              ! FIXME: what is the error by disallowing negative amounts?
              phases%phi(phase_idx)%x_mean(i) = phases%phi(phase_idx)%x_mean(i)&
                + phases%get_phi(phase_idx, node)*max(get_c_phase(phases, c_field(i), i, phase_idx, node), 0.0_wp)
            end do
          else
            phases%phi(phase_idx)%x_mean = 0
          end if
        else
          phases%phi(phase_idx)%amount = 0
          phases%phi(phase_idx)%x_mean = 0
        end if
      end if
      if(phase_idx /= pore_idx) then
        call tqgpi(phase_tuple_idx, phases%phi(phase_idx)%properties%name, current_equi)
        call tqgpi2(oc_phase_idx, compset_idx, phases%phi(phase_idx)%properties%name, current_equi)

        if(phases%phi(phase_idx)%amount <= 1.0_wp .or. .not. phases%phi(phase_idx)%properties%calphad_use_phase_conc) then
          call get_state_var_value('DG(' // phases%phi(phase_idx)%properties%name // ')', f0(phase_idx), encoded, current_equi)
          f0(phase_idx) = -f0(phase_idx)*constants%R_gas*temp
          call tqgetv('X', phase_tuple_idx, -1, n3, x_eq, current_equi, component_names)
          phases%phi(phase_idx)%x_mean = x_eq
        else
          f0(phase_idx) = 0.0_wp
          phases%phi(phase_idx)%x_mean = phases%phi(phase_idx)%x_mean/max(sum(phases%phi(phase_idx)%x_mean), delta_num)
        end if
        phases%phi(phase_idx)%x_mean = min(max(phases%phi(phase_idx)%x_mean, delta_num), 1.0_wp - delta_num)
        if(verbosity_level >= 5) call print_message('Phase ' // phases%phi(phase_idx)%properties%name // ' with amount ' &
          // convert_to_char(phases%phi(phase_idx)%amount) //' has the mean atomic fractions ' &
          // convert_to_char(phases%phi(phase_idx)%x_mean))

        allocate(mugrad(phases%n_species**2))
        allocate(mugrad_diag(phases%n_species))

        call get_phase_variance(oc_phase_idx, n1)
        phase_stoichiometric = n1 == 0
        if(.not. phase_stoichiometric) then
          if(verbosity_level >= 5) then
            call print_message('Calculating phase equilibrium for ' // phases%phi(phase_idx)%properties%name)
          end if
          call tqgdmat_th(phase_tuple_idx, temp*units%kelvin, pressure*units%pascal, phases%phi(phase_idx)%x_mean, cpot, &
            mugrad, current_equi, verbosity_level < 5)
          if(gx%bmperr /= 0) then
            print *, 'ERROR: opencalphad_equilib_phase: OC could not calculate Darken matrix!'
            print *, 'Conditions: x = ', phases%phi(phase_idx)%x_mean, &
              ', temp = ', temp, ', cpot = ', cpot, ', mugrad = ', mugrad
            ! FIXME: restore energy parameter output
            !do i = 1, n_species
            !    print *, 'Energy parameters ' // species_name(i) // ': c_equi = ', c_field(i)%equi_val,&
            !        ', p_coeff = ', c_field(i)%p_coeff, ', k_coeff = ', c_field(i)%k_coeff
            !end do
            print *, units%p_coeff, units%k_coeff, units%c
            error stop 
          end if
        end if
        call get_sublattice_number(phase_idx, nsub, current_equi)
        allocate(nkl(nsub))
        allocate(nsites(nsub))
        call tqgphc1(phase_tuple_idx, nsub, nkl, spix, yfrac, nsites, extra, current_equi)
        phases%phi(phase_idx)%properties%mol_per_formula_unit = extra(1)
        call tqgetv('NP', phase_tuple_idx, -1, n3, fetch_var, current_equi, component_names)
        phases%phi(phase_idx)%properties%equi_fraction = fetch_var(1)/phases%phi(phase_idx)%properties%mol_per_formula_unit &
          *phases%phi(phase_idx)%properties%V_m
        if(.not. phase_stoichiometric) then
          omega(:, phase_idx) = mugrad(::phases%n_species+1)/units%omega
          delta_mu(:, phase_idx) = (cpot - mu_global_eq)/units%chem_pot
        else
          if(verbosity_level >= 5) call print_message('Phase ' // phases%phi(phase_idx)%properties%name //&
            ' has fixed stoichiometry, setting p = p_max')
          omega(:, phase_idx) = sqrt(huge(0.0_wp))
          delta_mu(:, phase_idx) = 0.0_wp
        end if
        mole_equi(:, phase_idx) = phases%phi(phase_idx)%x_mean*phases%phi(phase_idx)%properties%mol_per_formula_unit  ![mol/mol_FU]
        deallocate(mugrad)
        deallocate(mugrad_diag)
        deallocate(nkl)
        deallocate(nsites)
      end if
    end if

    if(verbosity_level >= 5) call print_message('Awaiting OC data synchronization...')
    call co_max(omega)
    call co_max(mole_equi)
    call co_max(delta_mu)
    call co_max(f0)
    if(phase_idx > 0) then
      n_phase_total = phases%phi(phase_idx)%properties%equi_fraction
    else
      n_phase_total = 0
    end if
    call co_sum(n_phase_total)
    if(phase_idx > 0) phases%phi(phase_idx)%properties%equi_fraction = phases%phi(phase_idx)%properties%equi_fraction/n_phase_total
    do i = 1, phases%n_phases
      call co_broadcast(phases%phi(i)%amount, phases%phi(i)%properties%idx%image)
      call co_broadcast(phases%phi(i)%properties%equi_fraction, phases%phi(i)%properties%idx%image)
      ! FIXME: this should be called x_approximation (as the point at which approx is developed) rather than x_mean
      call co_broadcast(phases%phi(i)%x_mean, phases%phi(i)%properties%idx%image)
      call co_broadcast(phases%phi(i)%properties%mol_per_formula_unit, phases%phi(i)%properties%idx%image)
    end do
    if(verbosity_level >= 5) call print_message('OC data synchronized!')

    mu_global_eq = mu_global_eq/units%chem_pot
    !FIXME: does V_m factor into mu and/or omega?
    do concurrent(i = 1:phases%n_phases)
      if(phases%phi(i)%properties%V_m == 0) then
        conc_equi(:,i) = 0.0_wp
      else
        conc_equi(:,i) = mole_equi(:,i)/phases%phi(i)%properties%V_m
      end if
    end do
    call phases%set_energy_params(f0, mu_global_eq, delta_mu, omega, conc_equi)
  end subroutine opencalphad_equilib_phase

  pure function get_c_init(phases, c_gamma_opt) result(c_init)
    type(phase_set), intent(in) :: phases
    real(wp), intent(in), optional :: c_gamma_opt(phases%n_species)
    real(wp) :: c_init(phases%n_phases, phases%n_species)
    integer :: i, j

    do concurrent(i = 1:phases%n_species, j = 1:phases%n_phases)
      c_init(j,i) = phases%phi(j)%equi_val(i)
    end do
    if(present(c_gamma_opt)) then
      c_init(phases%n_phases,:) = c_gamma_opt
    end if
  end function get_c_init

  pure elemental real(wp) function get_df_order_deta(order_parameters, node) result(df_order_deta)
    type(order_parameter_set), intent(in) :: order_parameters
    integer, intent(in) :: node

    ! 'm_order*eta' could be included in implicit term, but doesn't seem to affect stability either way
    associate(m_order => order_parameters%eta%interface_data%get_m_order(), &
        eta => order_parameters%eta%real%vals(node), &
        sqsum => order_parameters%sqsum%vals(node), &
        gamma => order_parameters%eta%interface_data%get_gamma())
      df_order_deta = m_order*eta*(eta**2*(1 - 2*gamma) + 2*gamma*sqsum - 1)
    end associate
  end function get_df_order_deta

  pure real(wp) function get_df_chem_dc(c_field, phases, species_idx, node) result(mu_c)
    type(phase_set), intent(in) :: phases
    type(fft_field), intent(in) :: c_field(:)
    integer, intent(in) :: species_idx, node

    associate(equi_mu => phases%equi_mu(species_idx), &
        c => c_field(species_idx)%real%vals(node), &
        c_mix => phases%c_mix(species_idx)%real%vals(node), &
        p_mix => phases%p_mix(species_idx)%real%vals(node), &
        p_coeff => phases%phi(1)%p_coeff(species_idx))
      if(phases%equal_p_coeff(species_idx)) then
        mu_c = equi_mu + 2*(c - c_mix)*p_coeff
      else
        mu_c = equi_mu + 2*(c - c_mix)/p_mix
      end if
    end associate
  end function get_df_chem_dc    

  pure subroutine calc_df_order_deta(order_parameters, df_order_deta)
    type(order_parameter_set), intent(in) :: order_parameters
    type(time_dependent_real), intent(inout) :: df_order_deta
    character(*), parameter :: proc_name = 'calc_df_order_deta'
    integer :: node

    call order_parameters%sqsum%check_time_is_equal(order_parameters%eta%real, proc_name)
    if(df_order_deta%time == order_parameters%sqsum%time) return
    do concurrent(node = 1:order_parameters%grid%n_nodes_total)
      df_order_deta%vals(node) = get_df_order_deta(order_parameters, node)
    end do
    df_order_deta%time = order_parameters%sqsum%time
  end subroutine calc_df_order_deta

  pure subroutine calc_df_chem_dc(c_field, phases, species_idx, df_chem_dc)
    type(phase_set), intent(in) :: phases
    type(fft_field), intent(in) :: c_field(:)
    integer, intent(in) :: species_idx
    type(time_dependent_real), intent(inout) :: df_chem_dc
    character(*), parameter :: proc_name = 'calc_df_chem_dc'
    integer :: node, i

    do i = 2,phases%n_species
      call c_field(i)%real%check_time_is_equal(c_field(i-1)%real, proc_name)
    end do
    do i = 2,phases%n_phases - 1
      call phases%phi(i)%check_time_is_equal(phases%phi(i-1), proc_name)
    end do
    call c_field(1)%real%check_time_is_equal(phases%phi(1), proc_name)
    if(df_chem_dc%time == c_field(1)%real%time) return

    do concurrent(node = 1:phases%grid%n_nodes_total)
      df_chem_dc%vals(node) = get_df_chem_dc(c_field, phases, species_idx, node)
    end do

    df_chem_dc%time = c_field(1)%real%time
  end subroutine calc_df_chem_dc

  pure subroutine calc_df_chem_deta_common(&
      c_field, order_parameters, species_properties, phases, phase_idx, prestep_c, dt, df_chem_deta_common)
    type(phase_set), intent(in) :: phases
    integer, intent(in) :: phase_idx
    type(order_parameter_set), intent(in) :: order_parameters
    type(species_properties_t), intent(in) :: species_properties(:)
    type(fft_field), intent(in) :: c_field(:)
    type(time_dependent_real), intent(inout) :: df_chem_deta_common
    logical, intent(in) :: prestep_c
    real(wp), intent(in) :: dt
    real(wp) :: mu_average
    integer :: i, species_idx
    character(*), parameter :: proc_name = 'calc_df_chem_deta_common'

    do i = 2,phases%n_species
      call c_field(i)%real%check_time_is_equal(c_field(i-1)%real, proc_name)
    end do
    do i = 2,phases%n_phases - 1
      call phases%phi(i)%check_time_is_equal(phases%phi(i-1), proc_name)
    end do
    if(prestep_c) then
      call c_field(1)%real%check_time_diff_within_tolerance(phases%phi(1), dt, proc_name)
    else
      call phases%phi(1)%check_time_is_equal(c_field(1)%real, proc_name)
    end if
    call phases%phi(1)%check_time_is_equal(order_parameters%eta%real, proc_name)
    call phases%phi(1)%check_time_is_equal(phases%c_mix%real, proc_name)
    call phases%phi(1)%check_time_is_equal(phases%p_mix%real, proc_name)

    if(df_chem_deta_common%time == order_parameters%eta%real%time) return

    df_chem_deta_common%vals = phases%phi(phase_idx)%f0 - phases%f0_mix%real%vals
    ! do concurrent fails here in combination with associate on ifort
    ! also described here: 
    ! https://community.intel.com/t5/Intel-Fortran-Compiler/DO-CONCURRENT-and-ASSOCIATE-triggers-Runtime-Errors/td-p/1197288
    do species_idx = 1, phases%n_species
      associate(c_mix => phases%c_mix(species_idx)%real%vals, &
          p_mix => phases%p_mix(species_idx)%real%vals, &
          p_coeff => phases%phi(phase_idx)%p_coeff(species_idx), &
          equal_p_coeff => phases%equal_p_coeff(species_idx), &
          equi_val => phases%phi(phase_idx)%equi_val(species_idx), &
          c => c_field(species_idx)%real%vals)

        if(species_properties(species_idx)%instant_diffusion) then
          mu_average = sum(c - c_mix)/sum(p_mix)
          if(equal_p_coeff) then
            df_chem_deta_common%vals = df_chem_deta_common%vals + (c_mix - equi_val)*mu_average
          else
            df_chem_deta_common%vals = df_chem_deta_common%vals &
              + mu_average/2*(c + c_mix - 2*equi_val - mu_average/(2*p_coeff))
          end if
        else
          if(equal_p_coeff) then
            df_chem_deta_common%vals = df_chem_deta_common%vals + &
              2*p_coeff*(c - c_mix)*(c_mix - equi_val)
          else
            df_chem_deta_common%vals = df_chem_deta_common%vals + &
              (c - c_mix)/p_mix*(c + c_mix - 2*equi_val - (c - c_mix)/p_mix/p_coeff)
          end if
        end if
      end associate
    end do
    df_chem_deta_common%vals = 2/order_parameters%sqsum%vals*df_chem_deta_common%vals
    df_chem_deta_common%time = order_parameters%eta%real%time
  end subroutine calc_df_chem_deta_common

  pure subroutine calc_df_chem_deta(order_parameters, df_chem_deta_common, df_chem_deta)
    type(order_parameter_set), intent(in) :: order_parameters
    type(time_dependent_real), intent(in) :: df_chem_deta_common
    type(time_dependent_real), intent(inout) :: df_chem_deta
    character(*), parameter :: proc_name = 'calc_df_chem_deta'

    call df_chem_deta_common%check_time_is_equal(order_parameters%eta%real, proc_name)
    if(df_chem_deta%time == df_chem_deta_common%time) return
    df_chem_deta%vals = order_parameters%eta%real%vals*df_chem_deta_common%vals
    df_chem_deta%time = df_chem_deta_common%time
  end subroutine calc_df_chem_deta 

  pure real(wp) function get_c_phase(phases, c_field, species_idx, phase_idx, node) result(c_phase)
    type(phase_set), intent(in) :: phases
    type(fft_field), intent(in) :: c_field
    integer, intent(in) :: species_idx, phase_idx, node
    character(*), parameter :: proc_name = 'get_c_phase'

    !call c_field%real%check_time_is_equal(phases%c_mix%real, proc_name)
    !call phases%c_mix%real%check_time_is_equal(phases%p_mix%real, proc_name)
    associate(c => c_field%real%vals(node), &
        c_mix => phases%c_mix(species_idx)%real%vals(node), &
        p_mix => phases%p_mix(species_idx)%real%vals(node), &
        p_coeff => phases%phi(phase_idx)%p_coeff(species_idx), &
        equi_val => phases%phi(phase_idx)%equi_val(species_idx))
      if(phases%equal_p_coeff(species_idx)) then
        c_phase = equi_val + c - c_mix
      else
        c_phase = equi_val + (c - c_mix)/(p_coeff*p_mix)
      end if
    end associate
  end function get_c_phase

  subroutine calc_numeric_energy_derivatives(order_parameters, phases, c_field, species_idx, variable_team, df_dc, df_deta)
    type(order_parameter_set), intent(in) :: order_parameters
    type(phase_set), intent(in) :: phases
    type(fft_field), intent(in) :: c_field(:)
    integer, intent(in) :: species_idx
    type(team_type), intent(in) :: variable_team
    type(phase_field_energies), intent(inout) :: df_dc(:), df_deta
    type(phase_field_energies) :: f_plus, f_minus
    type(order_parameter_set) :: order_parameters_delta
    type(phase_set) :: phases_delta
    type(fft_field) :: c_field_delta(size(c_field))
    integer :: i
    real(wp), parameter :: delta = 1E-6_wp

    order_parameters_delta = order_parameters
    phases_delta = phases
    do i = 1, size(c_field)
      call c_field_delta(i)%init(c_field(i)%real%grid, c_field(i)%real%var_name, c_field(i)%keep_prev_fft)
    end do
    c_field_delta%real%time = c_field%real%time
    phases_delta = phases
    order_parameters_delta = order_parameters
    do i = 1, size(c_field)
      f_plus = df_dc(i)
      f_minus = df_dc(i)
      call f_plus%reset()
      call f_minus%reset()
      c_field_delta = c_field
      call c_field_delta(i)%reset()
      c_field_delta(i)%real%vals = c_field(i)%real%vals + delta
      c_field_delta(i)%real%time = c_field(i)%real%time
      call f_plus%calc_total(order_parameters_delta, phases_delta, c_field_delta, species_idx, variable_team)
      call c_field_delta(i)%reset()
      c_field_delta(i)%real%vals = c_field(i)%real%vals - delta
      c_field_delta(i)%real%time = c_field(i)%real%time
      call f_minus%calc_total(order_parameters_delta, phases_delta, c_field_delta, species_idx, variable_team)
      df_dc(i) = (f_plus - f_minus)/(2*delta)
    end do
    c_field_delta = c_field
    do i = 1, order_parameters%n_eta
      call phases_delta%reset()
      change team(variable_team)
      if(team_number() == eta_team_num) then
        f_plus = df_deta
        f_minus = df_deta
        order_parameters_delta = order_parameters
        call order_parameters_delta%reset()
        order_parameters_delta%eta%real%time = order_parameters%eta%real%time
        if(order_parameters%eta%idx == i) then
          order_parameters_delta%eta%real%vals = order_parameters%eta%real%vals + delta
        end if
        call order_parameters_delta%calc_phi(phases_delta)
      end if
      end team
      call f_plus%reset()
      call f_minus%reset()
      call phases_delta%phi_broadcast()
      call f_plus%calc_total(order_parameters_delta, phases_delta, c_field_delta, species_idx, variable_team)
      call phases_delta%reset()
      change team(variable_team)
      if(team_number() == eta_team_num) then
        call order_parameters_delta%reset()
        order_parameters_delta%eta%real%time = order_parameters%eta%real%time
        if(order_parameters%eta%idx == i) then
          order_parameters_delta%eta%real%vals = order_parameters%eta%real%vals - delta
        end if
        call order_parameters_delta%calc_phi(phases_delta)
      end if
      end team
      call phases_delta%phi_broadcast()
      call f_minus%calc_total(order_parameters_delta, phases_delta, c_field_delta, species_idx, variable_team)
      change team(variable_team)
      if(team_number() == eta_team_num) then
        if(order_parameters%eta%idx == i) df_deta = (f_plus - f_minus)/(2*delta)
      end if
      end team
    end do
  end subroutine calc_numeric_energy_derivatives
end module mod_mu
