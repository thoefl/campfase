! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_timer_set
    use mod_base_functions
    use mod_timer
    use mod_step_timer
    implicit none
    private
    
    character(*), parameter :: module_name = 'mod_timer_set'
    
    type, public :: timer_set
        type(timer), allocatable :: wall_timer, sim_timer
        type(step_timer), allocatable :: step_timer
        logical :: timers_coupled = .false.
        logical :: force = .false.
        logical :: initialized = .false.
        logical :: trigger_on_negative_step = .false.
    contains
        procedure, public :: new_schedule
        procedure, public :: poll
        procedure, public :: sync
        procedure, public :: is_due
        procedure, public :: schedule
        procedure, public :: is_enabled
    end type timer_set
    
    contains
    
    subroutine new_schedule(this, step)
        class(timer_set), intent(inout) :: this
        integer(int64), intent(in) :: step
        
        if(allocated(this%step_timer)) call this%step_timer%set_step_trigger(step)
        if(allocated(this%wall_timer)) call this%wall_timer%set_t_trigger(0.0_wp)
        if(allocated(this%sim_timer)) call this%sim_timer%set_t_trigger(0.0_wp)
        this%initialized = .true.
    end subroutine new_schedule
    
    pure elemental logical function is_due(this, step)
        class(timer_set), intent(in) :: this
        integer(int64), intent(in) :: step
        character(*), parameter :: proc_name = 'is_due'
    
        if(.not. this%initialized) call error_msg(proc_name, 'Not initialized', module_name_opt = module_name)
        is_due = .false.
        if(step < 0 .and. .not. this%trigger_on_negative_step) return
        if(this%force) then
            is_due = .true.
            return
        end if
        if(allocated(this%wall_timer)) is_due = is_due .or. this%wall_timer%elapsed
        if(allocated(this%sim_timer)) is_due = is_due .or. this%sim_timer%elapsed
        if(allocated(this%step_timer)) is_due = is_due .or. this%step_timer%is_due(step)
    end function is_due
    
    subroutine schedule(this, step)
        class(timer_set), intent(inout) :: this
        integer(int64), intent(in) :: step
        character(*), parameter :: proc_name = 'schedule'
        
        if(.not. this%initialized) call error_msg(proc_name, 'Not initialized', module_name_opt = module_name)
        this%force = .false.
        if(allocated(this%wall_timer)) then
            if(this%timers_coupled .or. this%wall_timer%elapsed) call this%wall_timer%arm()
        end if
        if(allocated(this%sim_timer)) then
            if(this%timers_coupled .or. this%sim_timer%elapsed) call this%sim_timer%arm()
        end if
        if(allocated(this%step_timer)) then
            if(this%timers_coupled .or. this%step_timer%is_due(step)) call this%step_timer%arm(step)
        end if
    end subroutine schedule
    
    pure elemental logical function is_enabled(this)
        class(timer_set), intent(in) :: this
        
        is_enabled = allocated(this%wall_timer) .or. allocated(this%sim_timer) .or. allocated(this%step_timer)
    end function is_enabled
    
    subroutine poll(this)
        class(timer_set), intent(inout) :: this
        
        if(allocated(this%wall_timer)) call this%wall_timer%poll()
        if(allocated(this%sim_timer)) call this%sim_timer%poll()
    end subroutine poll
    
    subroutine sync(this)
        class(timer_set), intent(inout) :: this
        
        if(allocated(this%wall_timer)) call this%wall_timer%sync()
        if(allocated(this%sim_timer)) call this%sim_timer%sync()
    end subroutine sync
end module mod_timer_set