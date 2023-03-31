! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_species_properties
    use mod_base_functions, only: wp
    use mod_arrhenius_model
    implicit none
    private
    
    type, public :: species_properties_t
        character(2) :: name
        type(arrhenius_model) :: diffusion_coefficient
        real(wp) :: x_init, x_oc_init
        logical :: instant_diffusion
    end type species_properties_t
end module mod_species_properties