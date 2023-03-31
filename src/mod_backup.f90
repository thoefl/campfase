! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_backup
    use mod_base_functions
    use mod_stopwatch
    use mod_spectral
    use mod_phase_set
    use mod_order_parameter_set
    use mod_nucleation
    implicit none
    private
    
    type, public :: backup_t
        integer(int64) :: step
        integer, allocatable :: randseed(:)
        type(manual_stopwatch) :: sim_clock
        type(fft_field), allocatable :: c_field(:)
        type(phase_set) :: phases
        type(order_parameter_set) :: order_parameters
        type(nucleation_data) :: nucleation
    contains
        procedure, public :: store
        procedure, public :: restore
    end type backup_t
    
    contains
    
    subroutine store(this, step, sim_clock, c_field, phases, order_parameters, nucleation, randseed)
        class(backup_t), intent(inout) :: this
        integer(int64), intent(in) :: step
        type(manual_stopwatch), intent(in) :: sim_clock
        type(fft_field), intent(in) :: c_field(:)
        type(phase_set), intent(in) :: phases
        type(order_parameter_set), intent(in) :: order_parameters
        type(nucleation_data), intent(in) :: nucleation
        integer, intent(in) :: randseed(:)
        
        this%step = step
        this%sim_clock = sim_clock
        this%c_field = c_field
        this%phases = phases
        this%order_parameters = order_parameters
        this%nucleation = nucleation
        this%randseed = randseed
    end subroutine store
   
    subroutine restore(this, step, sim_clock, c_field, phases, order_parameters, nucleation, randseed)
        class(backup_t), intent(in) :: this
        integer(int64), intent(out) :: step
        type(manual_stopwatch), intent(out) :: sim_clock
        type(fft_field), intent(out) :: c_field(:)
        type(phase_set), intent(out) :: phases
        type(order_parameter_set), intent(out) :: order_parameters
        type(nucleation_data), intent(out) :: nucleation
        integer, intent(out) :: randseed(:)
        
        step = this%step
        sim_clock = this%sim_clock
        c_field = this%c_field
        phases = this%phases
        order_parameters = this%order_parameters
        nucleation = this%nucleation
        randseed = this%randseed
    end subroutine restore
end module mod_backup
