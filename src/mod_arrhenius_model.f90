! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_arrhenius_model
    use mod_base_functions, only: wp
    implicit none
    private
    
    type, public :: arrhenius_model
        private
        real(wp) :: A, E_a_R_gas
    contains
        procedure, public :: get_value
    end type arrhenius_model
    
    interface arrhenius_model
        module procedure :: construct
    end interface arrhenius_model
    
    contains
    
    pure elemental type(arrhenius_model) function construct(A, E_a, R_gas) result(this)
        real(wp), intent(in) :: A, E_a, R_gas
        
        this%A = A
        this%E_a_R_gas = E_a/R_gas
    end function construct
    
    pure elemental real(wp) function get_value(this, temperature) result(res)
        class(arrhenius_model), intent(in) :: this
        real(wp), intent(in) :: temperature
        
        res = this%A*exp(-this%E_a_R_gas/temperature)
    end function get_value
end module mod_arrhenius_model