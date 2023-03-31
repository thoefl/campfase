! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_boundary_condition
    
    use mod_base_functions
    implicit none
    private
    
    integer(int8), parameter, public :: bc_periodic = 0
    integer(int8), parameter, public :: bc_neumann = 1
    integer(int8), parameter, public :: bc_dirichlet = 2
    
    type, public :: boundary_condition
        integer(int8) :: type_min, type_max
        real(wp) :: val_min, val_max
        logical :: is_periodic
    end type boundary_condition
    
    interface boundary_condition
        module procedure :: boundary_condition_constructor
    end interface boundary_condition
    
    contains
    
    pure elemental type(boundary_condition) function &
        boundary_condition_constructor(type_min, type_max, val_min_opt, val_max_opt) result(bc)
        integer(int8), intent(in) :: type_min, type_max
        real(wp), intent(in), optional :: val_min_opt, val_max_opt
        
        if(type_min == bc_periodic .and. type_max == bc_periodic) then
            bc%is_periodic = .true.
        else if(count([type_min == bc_periodic, type_max == bc_periodic]) == 1) then
            error stop &
                'ERROR: boundary_condition_constructor : Boundary conditions have to be either all periodic or all non-periodic!'
        else
            bc%is_periodic = .false.
        end if
        bc%type_min = type_min
        bc%type_max = type_max
        if(present(val_min_opt)) then
            bc%val_min = val_min_opt
        end if
        if(present(val_max_opt)) then
            bc%val_max = val_max_opt
        end if
    end function boundary_condition_constructor
end module mod_boundary_condition