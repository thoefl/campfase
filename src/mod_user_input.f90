! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_user_input
    use iso_fortran_env, only: input_unit, output_unit
    implicit none
    private
    
    public :: binary_user_choice
    
    contains
    
    subroutine binary_user_choice(true_input, false_input, default_choice, bool, prompt_opt)
        character(*), intent(in) :: true_input, false_input
        logical, intent(in) :: default_choice
        logical, intent(out) :: bool
        character(*), intent(in), optional :: prompt_opt
        character(:), allocatable :: user_choice
    
        if(present(prompt_opt)) then
            write(output_unit, '(a)', advance = 'no') prompt_opt
            write(output_unit, '(a)', advance = 'no') ' ('
            if(default_choice) write(output_unit, '(a)', advance = 'no') '['
            write(output_unit, '(a)', advance = 'no') true_input
            if(default_choice) write(output_unit, '(a)', advance = 'no') ']'
            write(output_unit, '(a)', advance = 'no') '/'
            if(.not. default_choice) write(output_unit, '(a)', advance = 'no') '['
            write(output_unit, '(a)', advance = 'no') false_input
            if(.not. default_choice) write(output_unit, '(a)', advance = 'no') ']'
            write(output_unit, '(a)', advance = 'no') ')'
        end if
        allocate(character(max(len(true_input), len(false_input))) :: user_choice)
        do
            read (input_unit, '(a)') user_choice
            if(trim(user_choice) == trim(false_input)) then
                bool = .false.
                return
            else if(trim(user_choice) == trim(true_input)) then
                bool = .true.
                return
            else if(trim(user_choice) == '') then
                bool = default_choice
                return
            else
                write(output_unit, '(a)', advance = 'no') 'Input not understood, please enter "' 
                if(default_choice) write(output_unit, '(a)', advance = 'no') '['
                write(output_unit, '(a)', advance = 'no') true_input
                if(default_choice) write(output_unit, '(a)', advance = 'no') ']'
                write(output_unit, '(a)', advance = 'no') '" or "'
                if(.not. default_choice) write(output_unit, '(a)', advance = 'no') '['
                write(output_unit, '(a)', advance = 'no') false_input
                if(.not. default_choice) write(output_unit, '(a)', advance = 'no') ']'
                write(output_unit, '(a)', advance = 'yes') '".'
            end if
        end do
    end subroutine binary_user_choice
end module mod_user_input