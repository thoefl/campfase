! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_unit_factors
    use mod_base_functions, only: wp
    implicit none
    private
    
    type, public :: unit_factor_t
        real(wp) :: V_m, G_m, D, meter, kelvin
        real(wp) :: second, mole_FU, mole, joule, minute, hour, newton, kilogram, gram, pascal
        real(wp) :: t, if_energy, m_order, kappa, M_eta, df_chem_deta, mu_c, k_coeff, p_coeff, chem_pot, omega, c, R_gas, N_A
    end type unit_factor_t
    
    interface unit_factor_t
        module procedure construct
    end interface unit_factor_t
    
    contains
    
    pure type(unit_factor_t) function construct(V_m, G_m, D, meter, kelvin) result(this)
        real(wp), intent(in) :: V_m, G_m, D, meter, kelvin
        
        this%V_m = V_m
        this%G_m = G_m
        this%D = D
        this%meter = meter
        this%kelvin = kelvin
        
        this%second = this%meter**2/this%D
        this%mole_FU = this%meter**3/this%V_m
        this%mole = this%mole_FU
        this%joule = this%G_m*this%mole

        this%t = this%second                                    ! [s]
        this%if_energy = this%joule/this%meter**2               ! [J/m^2]
        this%m_order = this%joule/this%meter**3                 ! [J/m^3]
        this%kappa = this%joule/this%meter                      ! [J/m]
        this%M_eta = this%meter**3/(this%joule*this%second)     ! [m^3/(J*s)]
        this%df_chem_deta = this%joule/this%meter**3            ! [J/m^3]
        this%mu_c = this%mole/this%meter**5                     ! [mol/m^5]
        this%k_coeff = this%joule/this%mole                     ! [J/mol]
        this%p_coeff = this%joule*this%meter**3/this%mole**2    ! [J*m^3/mol^2]
        this%chem_pot = this%joule/this%mole                    ! [J/mol]
        this%omega = this%joule/this%mole                       ! [J/mol]
        this%c = this%mole/this%meter**3                        ! [mol/m^3]
        this%R_gas = this%joule/this%mole/this%kelvin           ! [J/mol/K]
        this%N_A = 1/this%mole                                  ! [1/mole]

        this%minute = this%second/60
        this%hour = this%minute/60
        this%newton = this%joule/this%meter
        this%kilogram = this%newton/this%meter*this%second**2
        this%gram = this%kilogram*1000
        this%pascal = this%newton/this%meter**2
    end function construct
end module mod_unit_factors