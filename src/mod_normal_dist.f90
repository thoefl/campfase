! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

! Modified from https://fortran-lang.discourse.group/t/normal-random-number-generator/3724/15
! Released under MIT license (see https://fortran-lang.discourse.group/t/license-of-code-posted-to-the-fortran-discourse/188/29 )

module mod_normal_dist

use mod_base_functions

implicit none
private
public :: rnorm_box_muller, rnorm_box_muller_vec

contains
    
function rnorm_box_muller_vec(n, mu, sigma) result(variates)
    integer, intent(in) :: n     ! # of normal variates
    real(wp), intent(in), optional :: mu    ! target mean
    real(wp), intent(in), optional :: sigma ! target standard deviation
    real(wp) :: variates(n) ! normal variates
    integer :: i,j
    logical :: n_odd
    
    n_odd = mod(n,2) /= 0
    do i=1,n/2
        j = 2*i - 1
        variates(j:j+1) = rnorm_box_muller()
    end do
    if (n_odd) variates(n) = rnorm_box_muller_single_variate()
    if (present(sigma)) variates = sigma*variates
    if (present(mu)) variates = variates + mu
end function rnorm_box_muller_vec

function rnorm_box_muller() result(variates) ! coded formulas from https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
    ! return two uncorrelated standard normal variates
    real(wp) :: variates(2)
    real(wp) :: u(2), factor, arg
    do
        call random_number(u)
        if (u(1) > 0.0_wp) exit
    end do
    factor = sqrt(-2 * log(u(1)))
    arg = 2*pi*u(2)
    variates = factor * [cos(arg),sin(arg)]
end function rnorm_box_muller
!
function rnorm_box_muller_single_variate() result(variate)
    ! return a standard normal variate
    real(wp) :: variate
    real(wp) :: u(2), factor, arg
    do
        call random_number(u)
        if (u(1) > 0.0_wp) exit
    end do
    factor = sqrt(-2 * log(u(1)))
    arg = 2*pi*u(2)
    variate = factor * cos(arg)
end function rnorm_box_muller_single_variate
end module mod_normal_dist