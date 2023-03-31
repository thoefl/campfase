! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_stopwatch
  use mod_base_functions
  implicit none
  private
  
  character(*), parameter :: module_name = 'mod_stopwatch'
  
  type, public :: stopwatch
    private
    logical :: initialized = .false.
    logical :: running = .false.
  contains
    procedure, public :: init
    procedure, public :: start
    procedure, public :: get_time
    procedure, public :: stop
    procedure, public :: is_running
    procedure, public :: tic
    procedure, public :: toc
  end type stopwatch
  
  type, extends(stopwatch), public :: system_stopwatch
    private
    integer(int64) :: rate, count_running, count_idle
  end type system_stopwatch
  
  type, extends(stopwatch), public :: manual_stopwatch
    private
    real(wp) :: time
  contains
    procedure, public :: set_time
    procedure, public :: increment_time
  end type manual_stopwatch
  
  type, extends(stopwatch), public :: cpu_stopwatch
    private
    real(wp) :: time_idle, time_running
  end type cpu_stopwatch
  
  contains
  
  subroutine init(this)
    class(stopwatch), intent(inout) :: this
    
    select type(this)
    type is(system_stopwatch)
      this%count_running = 0
      call system_clock(count_rate = this%rate)
    type is(manual_stopwatch)
      this%time = 0.0_wp
    type is(cpu_stopwatch)
      this%time_running = 0.0_wp
    end select
    this%initialized = .true.
  end subroutine init
  
  subroutine start(this)
    class(stopwatch), intent(inout) :: this
    integer(int64) :: count_now
    real(wp) :: t_cpu_now
    
    if(this%running) return
    select type(this)
    type is(system_stopwatch)
      call system_clock(count = count_now)
      this%count_idle = count_now - this%count_running
    type is(cpu_stopwatch)
      call cpu_time(t_cpu_now)
      this%time_idle = t_cpu_now - this%time_running
    end select
    this%running = .true.
  end subroutine start
  
  real(wp) function get_time(this) result(time)
    class(stopwatch), intent(in) :: this
    integer(int64) :: count_now, count_running
    real(wp) :: t_cpu_now
    character(*), parameter :: proc_name = 'poll'
   
    time = 0.0_wp
    if(.not. this%initialized) call error_msg(proc_name, 'Not initialized', module_name_opt = module_name)
    select type(this)
    type is(system_stopwatch)
      call system_clock(count = count_now)
      if(this%running) then
        count_running = count_now - this%count_idle
      else
        count_running = this%count_running
      end if
      time = count_running/real(this%rate, wp)
    type is(manual_stopwatch)
      time = this%time
    type is(cpu_stopwatch)
      call cpu_time(t_cpu_now)
      if(this%running) then
        time = t_cpu_now - this%time_idle
      else
        time = this%time_running
      end if
    class default
      call error_msg(proc_name, 'Stopwatch type not defined', module_name_opt = module_name)
    end select
  end function get_time
  
  subroutine stop(this)
    class(stopwatch), intent(inout) :: this
    character(*), parameter :: proc_name = 'stop'
    integer(int64) :: count_now
    real(wp) :: t_cpu_now
    
    if(.not. this%running) return
    select type(this)
    type is(system_stopwatch)
      call system_clock(count = count_now)
      this%count_running = count_now - this%count_idle
    type is(cpu_stopwatch)
      call cpu_time(t_cpu_now)
      this%time_running = t_cpu_now - this%time_idle
    end select
    this%running = .false.
  end subroutine stop
  
  pure logical function is_running(this)
    class(stopwatch), intent(in) :: this
    
    is_running = this%running
  end function is_running
  
  pure subroutine set_time(this, time)
    class(manual_stopwatch), intent(inout) :: this
    real(wp), intent(in) :: time
    
    this%time = time
  end subroutine set_time
  
  pure subroutine increment_time(this, dtime)
    class(manual_stopwatch), intent(inout) :: this
    real(wp), intent(in) :: dtime
    
    this%time = this%time + dtime
  end subroutine increment_time

  subroutine tic(this, counter)
    class(stopwatch), intent(in) :: this
    real(wp), intent(inout) :: counter

    counter = counter - this%get_time()
  end subroutine tic

  subroutine toc(this, counter)
    class(stopwatch), intent(in) :: this
    real(wp), intent(inout) :: counter

    counter = counter + this%get_time()
  end subroutine toc
end module mod_stopwatch
