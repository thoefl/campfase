! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_phase_properties
    use mod_base_functions
    use mod_crystal_indices
    implicit none
    private
    
    type, public :: phase_properties_t
        character(16) :: name
        real(wp) :: V_m, mobility_factor, mol_per_formula_unit, equi_fraction
        logical :: is_compset, calphad_use_phase_conc
        character(128) :: default_comp
        real(wp) :: brightness
        type(crystal_indices) :: idx
        integer :: n_order
    contains
        procedure, public :: get_atomic_volume
        procedure, public :: get_lattice_parameter
    end type phase_properties_t
    
    contains
    
    pure real(wp) elemental function get_atomic_volume(this, N_A) result(atomic_volume)
        class(phase_properties_t), intent(in) :: this
        real(wp), intent(in) :: N_A
        character(*), parameter :: proc_name = 'get_atomic_volume'
        
        if(this%V_m == 0) then
            atomic_volume = huge(0.0_wp)
            return
        end if
        atomic_volume = this%V_m/(this%mol_per_formula_unit*N_A)
    end function get_atomic_volume
    
    ! FIXME: this assumes a primitive cubic lattice and is inaccurate in many cases
    pure real(wp) elemental function get_lattice_parameter(this, N_A) result(lattice_parameter)
        class(phase_properties_t), intent(in) :: this
        real(wp), intent(in) :: N_A
        character(*), parameter :: proc_name = 'get_lattice_parameter'
        
        if(this%V_m == 0) then
            lattice_parameter = huge(0.0_wp)
            return
        end if
        lattice_parameter = this%get_atomic_volume(N_A)**(1/3.0_wp)
    end function get_lattice_parameter
end module mod_phase_properties