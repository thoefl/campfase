! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_combinatorics
    use iso_fortran_env, only: int64
    use mod_base_functions, only: operator(//), error_msg
    use mod_rational_number
    implicit none
    private
    
    public :: binomial_coefficient
    
    contains
    
    pure elemental integer(int64) function binomial_coefficient(n, k) result(coeff)
        integer, intent(in) :: n, k
        type(rational_number) :: frac
        integer :: i
        character(*), parameter :: proc_name = 'binomial_coefficient'
    
        if(k > n) then
            coeff = 0
        else if(k == n .or. k == 0) then
            coeff = 1
        else
            frac = rational_number(n,k)
            do concurrent(i = 1:k-1)
                frac = frac*(n - i)/(k - i)
            end do
            if(frac%dividend > 1) call error_msg(proc_name, 'Integer division failed')
            if(frac%divisor > huge(coeff)) call error_msg(proc_name, 'Overflow error')
            coeff = frac%divisor
        end if
    end function binomial_coefficient
end module mod_combinatorics