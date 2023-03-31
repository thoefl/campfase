! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_phase_field_energies
    use mod_base_functions
    use mod_regular_grid
    use mod_time_dep_vars
    use mod_spectral
    use mod_phase_set
    use mod_order_parameter_set
    use mod_phase_properties
    use mod_parameters
    implicit none
    private
    
    type, public :: phase_field_energies
        type(regular_grid) :: grid
        type(time_dependent_real) :: chem, order, gradient, total
        real(wp) :: total_sum, total_sum_prev
    contains
        procedure, public :: init
        procedure, public :: calc_total
        procedure, public :: calc_chem
        procedure, public :: calc_order
        procedure, public :: calc_gradient
        procedure, public :: reset
    end type phase_field_energies
    
    public operator (-)
    public operator (/)
    
    interface operator (-)
        module procedure subtract
    end interface

    interface operator (/)
        module procedure scalar_div
    end interface
    
    contains
    
    subroutine init(this, name_prefix, name_suffix, grid)
        class(phase_field_energies), intent(inout) :: this
        type(regular_grid), intent(in) :: grid
        character(*), intent(in) :: name_prefix, name_suffix
    
        this%grid = grid
        call this%chem%init(grid, name_prefix // '_chem' // name_suffix)
        call this%gradient%init(grid, name_prefix // '_gradient' // name_suffix)
        call this%order%init(grid, name_prefix // '_order' // name_suffix)
        call this%total%init(grid, name_prefix // '_total' // name_suffix)
        this%total_sum = 0
        this%total_sum_prev = huge(0.0_wp)
    end subroutine init

    pure subroutine calc_chem(this, phases, c_field, species_idx)
        class(phase_field_energies), intent(inout) :: this
        type(phase_set), intent(inout) :: phases
        type(fft_field), intent(in) :: c_field
        integer, intent(in) :: species_idx
        character(*), parameter :: proc_name = 'calc_chem'
        integer :: i
    
        do i = 2, phases%n_phases - 1
            call phases%phi(i)%check_time_is_equal(phases%phi(i-1), proc_name)
        end do
        call phases%phi(1)%check_time_is_equal(c_field%real, proc_name)
    
        if(this%chem%time == c_field%real%time) return
        call phases%calc_f0_mix()
        call phases%calc_c_mix(species_idx)
        if(.not. phases%equal_p_coeff(species_idx)) call phases%calc_p_mix(species_idx)
    
        associate(f_chem => this%chem%vals, &
                    c => c_field%real%vals, &
                    equi_mu => phases%equi_mu(species_idx), &
                    p_coeff => phases%phi(1)%p_coeff(species_idx), &
                    c_mix => phases%c_mix(species_idx)%real%vals, &
                    p_mix => phases%p_mix(species_idx)%real, &
                    f0_mix => phases%f0_mix%real%vals)
                
            f_chem = f0_mix + equi_mu*c
            if(phases%equal_p_coeff(species_idx)) then
                f_chem = f_chem + (c - c_mix)**2*p_coeff
            else
                f_chem = f_chem + (c - c_mix)**2/p_mix%vals
            end if
        end associate
        this%chem%time = c_field%real%time
    end subroutine calc_chem

    subroutine calc_order(this, order_parameters)
        class(phase_field_energies), intent(inout) :: this
        type(order_parameter_set), intent(inout) :: order_parameters
        integer :: node
    
        call order_parameters%calc_sqsum()
        call order_parameters%calc_quadsum()
        if(this%order%time == order_parameters%eta%real%time) return
        do concurrent(node = 1:this%grid%n_nodes_total)
            this%order%vals(node) = f_order(order_parameters, node)
        end do
        this%order%time = order_parameters%eta%real%time
    end subroutine calc_order

    subroutine calc_gradient(this, order_parameters)
        class(phase_field_energies), intent(inout) :: this
        type(order_parameter_set), intent(inout) :: order_parameters
        character(*), parameter :: proc_name = 'calc_gradient'
        integer :: node
    
        call order_parameters%check_eta_is_valid(proc_name)
        if(this%gradient%time == order_parameters%eta%real%time) return
        call order_parameters%eta%calc_grad()
        do concurrent(node = 1:this%grid%n_nodes_total)
            this%gradient%vals(node) = order_parameters%eta%interface_data%get_kappa()/2*order_parameters%eta%nabla_square(node)
        end do
        this%gradient%time = order_parameters%eta%real%time
    end subroutine calc_gradient
    
    subroutine calc_total(this, order_parameters, phases, c_field, species_idx, variable_team)
        class(phase_field_energies), intent(inout) :: this
        type(order_parameter_set), intent(inout) :: order_parameters
        type(phase_set), intent(inout) :: phases
        type(fft_field), intent(in) :: c_field(:)
        integer, intent(in) :: species_idx
        type(team_type), intent(in) :: variable_team
    
        change team(variable_team)
            if(team_number() == eta_team_num) then
                call order_parameters%calc_phi(phases)
            end if
        end team
        call phases%phi_broadcast()
        change team(variable_team)
            if(team_number() == eta_team_num) then
                this%chem%vals = 0
                call this%calc_gradient(order_parameters)
                call this%calc_order(order_parameters)
            else if(species_idx > 0) then
                this%gradient%vals = 0
                call this%calc_chem(phases, c_field(species_idx), species_idx)
            end if
        end team
        call co_broadcast(this%order%vals, phases%porosity_img_idx)
        call co_sum(this%chem%vals)
        call co_sum(this%gradient%vals)
        call co_broadcast(this%order%time, phases%porosity_img_idx)
        call co_broadcast(this%chem%time, 1)
        call co_broadcast(this%gradient%time, phases%porosity_img_idx)
        if(.not. (this%chem%time == this%gradient%time .and. this%chem%time == this%order%time)) error stop &
            'ERROR: calc_total: Energies calculated at different times!'
        this%total%vals = this%chem%vals + this%gradient%vals + this%order%vals
        this%total%time = this%chem%time
        this%total_sum_prev = this%total_sum
        this%total_sum = sum(this%total%vals)
    end subroutine calc_total
    
    pure elemental type(phase_field_energies) function subtract(a, b) result(res)
        type(phase_field_energies), intent(in) :: a, b
    
        res = a
        res%chem%vals = a%chem%vals - b%chem%vals
        res%order%vals = a%order%vals - b%order%vals
        res%gradient%vals = a%gradient%vals - b%gradient%vals
        res%total%vals = a%total%vals - b%total%vals
        res%total_sum = a%total_sum - b%total_sum
        res%total_sum_prev = a%total_sum_prev - b%total_sum_prev
    end function subtract
    
    pure elemental type(phase_field_energies) function scalar_div(field, scalar) result(res)
        type(phase_field_energies), intent(in) :: field
        real(wp), intent(in) :: scalar
    
        res = field
        res%chem%vals = field%chem%vals/scalar
        res%order%vals = field%order%vals/scalar
        res%gradient%vals = field%gradient%vals/scalar
        res%total%vals = field%total%vals/scalar
        res%total_sum = field%total_sum/scalar
        res%total_sum_prev = field%total_sum_prev/scalar
    end function scalar_div
    
    pure elemental real(wp) function f_order(order_parameters, node)
        type(order_parameter_set), intent(in) :: order_parameters
        integer, intent(in) :: node
    
        associate(m_order => order_parameters%eta%interface_data%get_m_order(), &
                  eta => order_parameters%eta%real%vals(node), &
                  sqsum => order_parameters%sqsum%vals(node), &
                  quadsum => order_parameters%quadsum%vals(node), &
                  gamma => order_parameters%eta%interface_data%get_gamma())
            f_order = m_order/4*((1 - 2*gamma)*quadsum + 2*sqsum*(gamma*sqsum - 1) + 1)
        end associate
    end function f_order
    
    pure subroutine reset(this)
        class(phase_field_energies), intent(inout) :: this
        
        call this%chem%invalidate()
        call this%order%invalidate()
        call this%gradient%invalidate()
        call this%total%invalidate()
    end subroutine reset
end module mod_phase_field_energies
