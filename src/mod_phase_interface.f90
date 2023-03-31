! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_phase_interface

    use mod_base_functions
    implicit none
    private
    
    type, public :: phase_interface
    private
        real(wp) :: width, energy, gamma
        real(wp) :: kappa, m_order
    contains
        procedure, public :: init
        procedure, public :: get_width
        procedure, public :: get_energy
        procedure, public :: get_gamma
        procedure, public :: get_kappa
        procedure, public :: get_m_order
    end type phase_interface
    
    contains
    
    pure subroutine init(this, width, energy, gamma)
        class(phase_interface), intent(inout) :: this
        real(wp), intent(in) :: width, energy, gamma
        
        this%width = width
        this%energy = energy
        this%gamma = gamma
        this%m_order = get_if_doublewell(energy, width, gamma)
        this%kappa = get_if_kappa(energy, width)
    end subroutine init
    
    pure elemental real(wp) function get_width(this) result(res)
        class(phase_interface), intent(in) :: this
        
        res = this%width
    end function get_width
    
    pure elemental real(wp) function get_energy(this) result(res)
        class(phase_interface), intent(in) :: this
        
        res = this%energy
    end function get_energy
    
    pure elemental real(wp) function get_gamma(this) result(res)
        class(phase_interface), intent(in) :: this
        
        res = this%gamma
    end function get_gamma
    
    pure elemental real(wp) function get_kappa(this) result(res)
        class(phase_interface), intent(in) :: this
        
        res = this%kappa
    end function get_kappa
    
    pure elemental real(wp) function get_m_order(this) result(res)
        class(phase_interface), intent(in) :: this
        
        res = this%m_order
    end function get_m_order
    
    ! Taken from Moelans et al., 2008: https://doi.org/10.1103/PhysRevB.78.024113, p.11
    pure elemental real(wp) function get_if_doublewell(energy, width, gamma) &
        result(if_doublewell)
        real(wp), intent(in) :: energy, width, gamma
    
        if_doublewell = 0.75_wp/get_f_saddle(gamma)*energy/width
    end function get_if_doublewell

    ! Taken from Moelans et al., 2008: https://doi.org/10.1103/PhysRevB.78.024113, p.4
    pure elemental real(wp) function get_f_saddle(gamma) result(f_saddle)
        real(wp), intent(in) :: gamma
    
        f_saddle = (2*gamma - 1)/(4*(2*gamma + 1))
    end function get_f_saddle
    
    ! Taken from Moelans et al., 2008: https://doi.org/10.1103/PhysRevB.78.024113, p.11
    pure elemental real(wp) function get_if_kappa(energy, width) result(kappa)
        real(wp), intent(in) :: energy, width
    
        kappa = 0.75_wp*width*energy
    end function get_if_kappa
end module mod_phase_interface