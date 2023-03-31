! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_step_timer
    use mod_base_functions
    implicit none
    private
    
    character(*), parameter :: module_name = 'mod_step_timer'
    
    type, public :: step_timer
        private
        integer(int64) :: step_trigger, dstep_trigger
        logical :: late_allowed = .false.
        logical :: initialized = .false.
    contains
        procedure, public :: init
        procedure, public :: arm
        procedure, public :: is_due
        procedure, public :: get_dstep_trigger
        procedure, public :: get_step_trigger
        procedure, public :: set_step_trigger
    end type step_timer
    
    contains
    
    pure subroutine init(this, dstep_trigger, late_allowed_opt)
        class(step_timer), intent(inout) :: this
        integer(int64), intent(in) :: dstep_trigger
        logical, intent(in), optional :: late_allowed_opt
        
        if(present(late_allowed_opt)) this%late_allowed = late_allowed_opt
        this%dstep_trigger = dstep_trigger
        this%step_trigger = huge(this%step_trigger)
        this%initialized = .true.
    end subroutine init
    
    pure subroutine arm(this, step_now)
        class(step_timer), intent(inout) :: this
        integer(int64), intent(in) :: step_now
        character(*), parameter :: proc_name = 'arm'
    
        if(.not. this%initialized) call error_msg(proc_name, 'Timer not initialized', module_name_opt = module_name)
        this%step_trigger = step_now + this%dstep_trigger
    end subroutine arm
    
    pure elemental logical function is_due(this, step_now)
        class(step_timer), intent(in) :: this
        integer(int64), intent(in) :: step_now
        character(*), parameter :: proc_name = 'is_due'
        
        if(.not. this%initialized) call error_msg(proc_name, 'Timer not initialized', module_name_opt = module_name)
        if(step_now == this%step_trigger) then
            is_due = .true.
            return
        else if(step_now > this%step_trigger) then
            if(.not. this%late_allowed) then
                call error_msg(proc_name, 'Timer is late: trigger = ' // convert_to_char(this%step_trigger) &
                    // ' < now = ' // convert_to_char(step_now) // ' (not allowed per setting "late_allowed")' &
                    , module_name_opt = module_name)
            else
                is_due = .true.
                return
            end if
        else
            is_due = .false.
        end if
    end function is_due
    
    pure elemental integer(int64) function get_dstep_trigger(this)
        class(step_timer), intent(in) :: this
        
        get_dstep_trigger = this%dstep_trigger
    end function get_dstep_trigger
    
    pure elemental integer(int64) function get_step_trigger(this)
        class(step_timer), intent(in) :: this
        
        get_step_trigger = this%step_trigger
    end function get_step_trigger
    
    pure elemental subroutine set_step_trigger(this, step_trigger)
        class(step_timer), intent(inout) :: this
        integer(int64), intent(in) :: step_trigger
        
        this%step_trigger = step_trigger
    end subroutine set_step_trigger
end module mod_step_timer