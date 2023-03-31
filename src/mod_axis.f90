! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_axis
    use mod_base_functions
    use mod_boundary_condition
    implicit none
    private
    
    type, public :: axis
        character(:), allocatable :: label
        integer :: n_nodes
        real(wp) :: spacing, length
        type(boundary_condition) :: bound_cond
    end type axis
    
    interface axis
        module procedure :: axis_constructor
    end interface
    
    contains
    
    pure elemental type(axis) function axis_constructor(n_nodes, spacing, bound_cond, label) result(this)
        integer, intent(in) :: n_nodes
        real(wp), intent(in) :: spacing
        type(boundary_condition), intent(in) :: bound_cond
        character(*), intent(in) :: label
        
        this%n_nodes = n_nodes
        this%spacing = spacing
        this%bound_cond = bound_cond
        this%label = label
        this%length = n_nodes*spacing
    end function axis_constructor
end module mod_axis