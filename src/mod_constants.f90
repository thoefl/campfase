! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_constants
    use mod_realtype
    use mod_unit_factors
    implicit none
    private
    
    type, public :: constants_t
        real(wp) :: R_gas, zero_celsius, N_A, k_B
    end type constants_t
    
    interface constants_t
        module procedure :: construct
    end interface constants_t
    
    contains
    
    pure type(constants_t) function construct(units) result(this)
        type(unit_factor_t), intent(in) :: units
        
        this%R_gas = 8.31446261815324_wp/units%R_gas
        this%zero_celsius = 273.15_wp
        this%N_A = 6.02214076E23_wp/units%N_A
        this%k_B = this%R_gas/this%N_A
    end function construct
end module mod_constants
