! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_phase_set
    use mod_realtype
    use mod_base_functions
    use mod_regular_grid
    use mod_time_dep_vars
    use mod_spectral
    use mod_phase_interface
    use mod_phase_properties
    use mod_parameters
    implicit none
    private
    
    character(*), parameter :: module_name = 'mod_phase_set'
    
    type, extends(time_dependent_real) :: phase_field
        type(phase_properties_t), pointer :: properties => null()
        real(wp), allocatable, dimension(:) :: omega
        real(wp), allocatable, dimension(:) :: p_coeff, k_coeff, equi_val
        real(wp), allocatable, dimension(:) :: x_mean
        real(wp) :: f0, amount
    end type phase_field

    type, public :: phase_set
        type(regular_grid) :: grid
        integer :: n_species, n_phases
        type(phase_field), allocatable, dimension(:) :: phi
        type(fft_field), allocatable, dimension(:) :: c_mix, p_mix, log_p_mix
        type(fft_field), allocatable :: f0_mix
        real(wp), allocatable, dimension(:) :: equi_mu
        logical, allocatable, dimension(:) :: equal_p_coeff
        character(2), dimension(:), allocatable :: species_name
        type(chemical_energy_model) :: chem_model
        integer :: porosity_img_idx, grain_img_idx
    contains
        procedure, public :: init
        procedure, public :: get_phi
        procedure, public :: phi_broadcast
        procedure, public :: phi_mix
        procedure, public :: calc_c_mix
        procedure, public :: calc_p_mix
        procedure, public :: calc_log_p_mix
        procedure, public :: calc_f0_mix
        procedure, public :: set_energy_params
        procedure, public :: get_p_max
        procedure, public :: get_brightness
        procedure, public :: output_composite_image
        procedure, public :: reset
    end type phase_set
    
    contains
    
    subroutine init(&
        this, grid, n_species, n_phases, species_name, p_settings, chem_model)
        class(phase_set), intent(inout) :: this
        type(phase_properties_t), intent(in), target :: p_settings(:)
        type(chemical_energy_model), intent(in) :: chem_model
        type(regular_grid), intent(in) :: grid
        integer, intent(in) :: n_species, n_phases
        character(*), dimension(:), intent(in) :: species_name
        integer :: i, species
    
        this%grid = grid
        this%n_species = n_species
        this%n_phases = n_phases
        allocate(this%phi(n_phases))
        allocate(this%c_mix(n_species))
        allocate(this%p_mix(n_species))
        allocate(this%log_p_mix(n_species))
        allocate(this%f0_mix)
        allocate(this%equal_p_coeff(n_species))
        allocate(this%species_name(n_species))
        allocate(this%equi_mu(n_species))
        this%species_name = species_name
        this%porosity_img_idx = n_species + 1
        this%grain_img_idx = n_species + n_phases
        this%chem_model = chem_model
    
        ! FIXME: initializing all of these is probably unnecessary
        call this%f0_mix%init(grid, f0_mix_name)
        do species = 1, n_species
            call this%c_mix(species)%init(grid, c_mix_name // '_' // this%species_name(species))
            call this%p_mix(species)%init(grid, p_mix_name // '_' // this%species_name(species))
            call this%log_p_mix(species)%init(grid, log_p_mix_name // '_' // this%species_name(species))
        end do
    
        do i = 1, n_phases
            associate(phase => this%phi(i))
                phase%properties => p_settings(i)
                if(i < n_phases) call phase%init(grid, phi_name // '_' // phase%properties%name)
                allocate(phase%omega(n_species))
                allocate(phase%p_coeff(n_species))
                allocate(phase%k_coeff(n_species))
                allocate(phase%equi_val(n_species))
                allocate(phase%x_mean(n_species))
                phase%omega = 0
                phase%p_coeff = 0
                phase%k_coeff = 0
                phase%equi_val = 0
                phase%properties%mol_per_formula_unit = 0
                phase%amount = 0
                phase%x_mean = 0
            end associate
        end do
    
        this%equal_p_coeff = .false.
        this%equi_mu = 0.0_wp
    end subroutine init
        
    pure elemental real(wp) function get_phi(this, phase_idx, node)
        class(phase_set), intent(in) :: this
        integer, intent(in) :: phase_idx, node
        integer :: i
    
        if(phase_idx < this%n_phases) then
            get_phi = this%phi(phase_idx)%vals(node)
        else if(phase_idx == this%n_phases) then
            get_phi = 1
            do concurrent(i = 1:this%n_phases - 1)
                get_phi = get_phi - this%phi(i)%vals(node)
            end do
        else
            error stop 'ERROR: get_phi: phase_idx is greater than amount of phases!'
        end if
    end function get_phi

    pure real(wp) function phi_mix(this, mix_vals, node)
        class(phase_set), intent(in) :: this
        real(wp), intent(in) :: mix_vals(:)
        integer, intent(in) :: node
        integer :: phase_idx
    
        phi_mix = mix_vals(this%n_phases)
        do concurrent(phase_idx = 1:this%n_phases - 1)
            phi_mix = phi_mix + this%phi(phase_idx)%vals(node)*(mix_vals(phase_idx) - mix_vals(this%n_phases))
        end do
    end function phi_mix
    
    subroutine phi_broadcast(this)
        class(phase_set), intent(inout) :: this
        integer :: phase_idx
    
        do phase_idx = 1, this%n_phases - 1
            call co_broadcast(this%phi(phase_idx)%vals, this%n_species + phase_idx)
            call co_broadcast(this%phi(phase_idx)%time, this%n_species + phase_idx)
        end do
    end subroutine phi_broadcast
  
    pure subroutine set_energy_params(this, f0, equi_mu, delta_mu, omega, expansion_val)
        class(phase_set), intent(inout) :: this
        real(wp), intent(in) :: f0(:), equi_mu(:), delta_mu(:,:), omega(:,:), expansion_val(:,:)
        integer :: phase_idx, species_idx
    
        this%equi_mu = equi_mu
        do phase_idx = 1, this%n_phases
            associate(phase => this%phi(phase_idx))
                phase%f0 = f0(phase_idx)
                phase%k_coeff = delta_mu(:, phase_idx)
                phase%omega = omega(:, phase_idx)
                if(phase%properties%V_m == 0.0_wp) then
                    phase%p_coeff = this%chem_model%p_coeff_max
                    phase%equi_val = expansion_val(:, phase_idx)
                else
                    phase%p_coeff = min(max(phase%omega/2*phase%properties%V_m,&
                    this%chem_model%p_coeff_min), this%chem_model%p_coeff_max)
                    phase%equi_val = expansion_val(:, phase_idx) - delta_mu(:, phase_idx)/(phase%omega*phase%properties%V_m)
                end if
            end associate
        end do
        this%equal_p_coeff = .true.
        !FIXME: this assumes that n_phases > 1
        do concurrent(phase_idx = 2:this%n_phases, species_idx = 1:this%n_species)
            this%equal_p_coeff(species_idx) = this%equal_p_coeff(species_idx)&
                .and. this%phi(phase_idx)%p_coeff(species_idx) == this%phi(phase_idx - 1)%p_coeff(species_idx)
        end do
    end subroutine set_energy_params
    
    pure elemental subroutine calc_c_mix(this, species_idx)  
        class(phase_set), intent(inout) :: this
        integer, intent(in) :: species_idx
        real(wp) :: c_mix_phase(this%n_species)
        character(*), parameter :: proc_name = 'phase_set_calc_c_mix'
        integer :: phase, node
        !dir$ vector aligned
    
        call this%phi(:this%n_phases - 1)%check_is_valid(proc_name)
        if(this%c_mix(species_idx)%real%time_is_equal(this%phi(1))) return
        do concurrent(phase = 1:this%n_phases)
            c_mix_phase(phase) = this%phi(phase)%equi_val(species_idx)
        end do
        do concurrent(node = 1:this%grid%n_nodes_total)
            this%c_mix(species_idx)%real%vals(node) = this%phi_mix(c_mix_phase, node)
        end do
        this%c_mix(species_idx)%real%time = this%phi(1)%time
    end subroutine calc_c_mix
    
    pure elemental subroutine calc_p_mix(this, species_idx)  
        class(phase_set), intent(inout) :: this
        integer, intent(in) :: species_idx
        real(wp) :: inv_p_coeff(this%n_species)
        character(*), parameter :: proc_name = 'phase_set_calc_p_mix'
        integer :: phase, node
        !dir$ vector aligned
    
        call this%phi(:this%n_phases - 1)%check_is_valid(proc_name)
        if(this%p_mix(species_idx)%real%time_is_equal(this%phi(1))) return
        do concurrent(phase = 1:this%n_phases)
            inv_p_coeff(phase) = 1/this%phi(phase)%p_coeff(species_idx)
        end do
        do concurrent(node = 1:this%grid%n_nodes_total)
            this%p_mix(species_idx)%real%vals(node) = this%phi_mix(inv_p_coeff, node)
        end do
        this%p_mix(species_idx)%real%time = this%phi(1)%time
    end subroutine calc_p_mix
    
    pure elemental subroutine calc_log_p_mix(this, species_idx)  
        class(phase_set), intent(inout) :: this
        integer, intent(in) :: species_idx
        character(*), parameter :: proc_name = 'phase_set_calc_log_p_mix'
        !dir$ vector aligned
    
        call this%calc_p_mix(species_idx)
        if(this%log_p_mix(species_idx)%real%time_is_equal(this%p_mix(species_idx)%real)) return
        this%log_p_mix(species_idx)%real%vals = log(this%p_mix(species_idx)%real%vals)
        this%log_p_mix(species_idx)%real%time = this%p_mix(species_idx)%real%time
    end subroutine calc_log_p_mix
    
    pure elemental subroutine calc_f0_mix(this)
        class(phase_set), intent(inout) :: this
        character(*), parameter :: proc_name = 'phase_set_calc_f0_mix'
        integer :: node
        !dir$ vector aligned
    
        call this%phi(:this%n_phases - 1)%check_is_valid(proc_name)
        if(this%f0_mix%real%time_is_equal(this%phi(1))) return
        do concurrent(node = 1:this%grid%n_nodes_total)
            this%f0_mix%real%vals(node) = this%phi_mix(this%phi%f0, node)
        end do
        this%f0_mix%real%time = this%phi(1)%time
    end subroutine calc_f0_mix
    
    pure real(wp) function get_p_max(this) result(p_max)
        class(phase_set), intent(in) :: this
        integer :: species, phase
    
        p_max = -huge(p_max)
        do concurrent(species = 1:this%n_species, phase = 1:this%n_phases)
            p_max = max(p_max, this%phi(phase)%p_coeff(species))
        end do
    end function get_p_max
    
    pure function get_brightness(this) result(brightness)
        class(phase_set), intent(in) :: this
        real(wp), dimension(this%n_phases) :: brightness
        integer :: phase
        
        do concurrent(phase = 1:this%n_phases)
            brightness(phase) = this%phi(phase)%properties%brightness
        end do
    end function get_brightness
    
    subroutine output_composite_image(this, sqsum)
        class(phase_set), intent(in) :: this
        type(time_dependent_real), intent(in) :: sqsum
        real(wp), allocatable, dimension(:) :: brightness, composite
        character(*), parameter :: composite_name = 'composite'
        character(*), parameter :: proc_name = 'output_composite_image'
        integer :: node, i
        
        do i = 2, this%n_phases - 1
            call this%phi(i)%check_time_is_equal(this%phi(i-1), proc_name)
        end do
        call this%phi(1)%check_time_is_equal(sqsum, proc_name)
        brightness = this%get_brightness()
        allocate(composite(this%grid%n_nodes_total))
        composite = 0
        do concurrent(node = 1:this%grid%n_nodes_total)
            composite(node) = sqsum%vals(node)*this%phi_mix(brightness, node)
        end do
        call output_image(composite, composite_name, this%grid%axes(1)%n_nodes, a_min_opt = 0.0_wp, a_max_opt = 1.0_wp)
    end subroutine output_composite_image
    
    pure subroutine reset(this)
        class(phase_set), intent(inout) :: this
        integer :: i
        
        if(allocated(this%f0_mix)) call this%f0_mix%reset()
        if(allocated(this%c_mix)) then
            do i = 1, size(this%c_mix)
                call this%c_mix(i)%reset()
            end do
        end if
        if(allocated(this%p_mix)) then
            do i = 1, size(this%p_mix)
                call this%p_mix(i)%reset()
            end do
        end if
        if(allocated(this%log_p_mix)) then
            do i = 1, size(this%log_p_mix)
                call this%log_p_mix(i)%reset()
            end do
        end if
        if(allocated(this%phi)) then
            do i = 1, size(this%phi)
                call this%phi(i)%invalidate()
            end do
        end if
    end subroutine reset
end module mod_phase_set
