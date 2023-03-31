! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_initial_field_settings
    use mod_base_functions
    use mod_regular_grid
    use mod_solid_shape
    use mod_polycrystal
    use mod_phase_properties
    implicit none
    private

    type, public :: initial_field_settings
        integer(int8) :: eta_type
        type(solid_shape) :: solid_part
        real(wp) :: spinodal_noise
        real(wp), dimension(:), allocatable :: eta_shift
        type(polycrystal_t), allocatable :: polycrystal
        integer :: n_precipitates
        type(phase_properties_t), pointer :: matrix_phase, precipitate_phase, pore_phase
    end type initial_field_settings

end module mod_initial_field_settings