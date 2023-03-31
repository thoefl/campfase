! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_time_step
  use mod_base_functions
  use mod_timer_set
  implicit none
  private

  type :: adaptive_settings
    logical :: enabled, use_bdfab
    real(wp) :: rel_err_tol, rel_err_max, max_factor, min_factor, max_factor_osc, bdfab_switch_rel_error_max
    real(wp) :: local_err = 0.0_wp
    real(wp) :: global_err_max
    type(timer_set) :: timer
    integer :: n_steps_richardson
    integer :: n_min_factor = 0
    integer :: substep = 0
    integer :: n_local_limiting = 0
  end type adaptive_settings

  type :: step_factor_t
    real(wp) :: val, max, bdfab_max
    real(wp) :: mean = 0.0_wp
    integer(int64) :: step = 0
  end type step_factor_t

  type, public :: time_step
    real(wp) :: val, init_val, saved_val
    type(adaptive_settings) :: adaptive
    logical :: prestep_c, bdfab_enabled, force_bdfab, use_bdfab
    type(step_factor_t) :: step_factor
    integer :: bdfab_stage = 0
  contains
    procedure, public :: set_step_factor
  end type time_step

contains

  subroutine set_step_factor(this, step_now, step_factor_now)
  class(time_step), intent(inout) :: this
    integer(int64), intent(in) :: step_now
    real(wp), intent(in) :: step_factor_now

    this%step_factor%max = max(this%step_factor%max, step_factor_now)
    if(this%use_bdfab) then
      this%step_factor%bdfab_max = max(this%step_factor%bdfab_max, step_factor_now)
    end if
    if(step_now == 0) then
      this%step_factor%mean = step_factor_now
    else
      this%step_factor%mean = &
        ((step_now - this%step_factor%step)*step_factor_now + this%step_factor%step*this%step_factor%mean)/step_now
    end if
    this%step_factor%step = step_now
    this%step_factor%val = step_factor_now
  end subroutine set_step_factor
end module mod_time_step
