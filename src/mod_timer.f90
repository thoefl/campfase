! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_timer
    use mod_base_functions
    use mod_stopwatch
    implicit none
    private
    
    character(*), parameter :: module_name = 'mod_timer'
    
    type, public :: timer
        private
        class(stopwatch), pointer :: clock
        real(wp) :: t_trigger = huge(0.0_wp)
        real(wp) :: dt_trigger
        logical, public :: elapsed = .false.
    contains
        procedure, public :: init
        procedure, public :: arm
        procedure, public :: poll
        procedure, public :: sync
        procedure, public :: get_dt_until_trigger
        procedure, public :: get_dt_trigger
        procedure, public :: get_t_trigger
        procedure, public :: set_t_trigger
    end type timer
    
    contains
        
    subroutine init(this, clock, dt_trigger)
        class(timer), intent(inout) :: this
        class(stopwatch), intent(in), target :: clock
        real(wp), intent(in) :: dt_trigger
        
        this%clock => clock
        this%dt_trigger = dt_trigger
        this%elapsed = dt_trigger <= 0.0_wp
    end subroutine init
    
    subroutine arm(this)
        class(timer), intent(inout) :: this
        character(*), parameter :: proc_name = 'arm'
        
        if(.not. this%clock%is_running()) call error_msg(proc_name, 'Timer not initialized', module_name_opt = module_name)
        this%t_trigger = this%clock%get_time() + this%dt_trigger
    end subroutine arm
    
    subroutine poll(this)
        class(timer), intent(inout) :: this
        
        this%elapsed = this%t_trigger <= this%clock%get_time()
    end subroutine poll
    
    subroutine sync(this)
        class(timer), intent(inout) :: this
        
        call co_any(this%elapsed)
    end subroutine sync
    
    pure elemental real(wp) function get_dt_trigger(this)
        class(timer), intent(in) :: this
        
        get_dt_trigger = this%dt_trigger
    end function get_dt_trigger
    
    pure elemental real(wp) function get_t_trigger(this)
        class(timer), intent(in) :: this
        
        get_t_trigger = this%t_trigger
    end function get_t_trigger
    
    subroutine set_t_trigger(this, t_trigger)
        class(timer), intent(inout) :: this
        real(wp), intent(in) :: t_trigger
        
        this%t_trigger = t_trigger
        call this%poll()
        call this%sync()
    end subroutine set_t_trigger
    
    real(wp) function get_dt_until_trigger(this) result(dt)
        class(timer), intent(in) :: this
        
        dt = this%t_trigger - this%clock%get_time()
    end function get_dt_until_trigger
end module mod_timer