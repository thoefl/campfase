! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_array_index
    implicit none
    private

    public :: count_true_false_transitions, mask_to_sections, sections_overlap, section_mask_overlap, sections_to_mask

    type, public :: array_index
        integer :: first, last
    end type array_index

    type, public, extends(array_index) :: strided_array_index
        integer :: stride
    end type strided_array_index

    contains

    pure integer function count_true_false_transitions(mask) result(count)
        logical, intent(in) :: mask(:)
        integer :: i

        count = 0
        if(all(mask) .or. .not. any(mask)) return
        do concurrent(i = 2:size(mask))
            if(mask(i) .neqv. mask(i - 1)) count = count + 1 
        end do
    end function count_true_false_transitions
    
    pure function mask_to_sections(mask) result(sections)
        logical, intent(in) :: mask(:)
        type(array_index), dimension(:), allocatable :: sections
        integer :: node, n_sections, n_trans, section_idx
        logical :: inside
        character(*), parameter :: proc_name = 'mask_to_sections'

        n_trans = count_true_false_transitions(mask)
        
        if(mod(n_trans, 2) /= 0) then
            n_sections = n_trans/2
        else
            n_sections = (n_trans + 1)/2
        end if

        if(mask(1)) n_sections = n_sections + 1

        allocate(sections(n_sections))
        if(n_sections == 0) return

        section_idx = 1
        inside = mask(1)
        if(inside) sections(1)%first = 1
        do node = 2, size(mask)
           if(inside .eqv. mask(node)) then
               cycle
           else if(inside) then
               sections(section_idx)%last = node - 1
               inside = .false.
               section_idx = section_idx + 1
           else
               sections(section_idx)%first = node
               inside = .true.
           end if
        end do
        if(inside) then
            sections(section_idx)%last = size(mask)
            section_idx = section_idx + 1
        end if
        if(section_idx - 1 /= n_sections) then
            error stop 'ERROR: ' // proc_name // ': Number of sections found inconsistent with expected amount!'
        end if
    end function mask_to_sections

    pure function sections_to_mask(sections, mask_size) result(mask)
        type(array_index), dimension(:), intent(in) :: sections
        integer, intent(in) :: mask_size
        logical, dimension(mask_size) :: mask
        integer :: idx

        mask = .false.
        do concurrent(idx = 1:size(sections))
            mask(sections(idx)%first:sections(idx)%last) = .true.
        end do
    end function sections_to_mask

    pure logical function sections_overlap(a, b)
        type(array_index), dimension(:), intent(in) :: a, b
        integer :: a_idx, b_idx

        sections_overlap = .false.
        do a_idx = 1, size(a)
            do b_idx = 1, size(b)
                if(a(a_idx)%last >= b(b_idx)%first .and. a(a_idx)%first <= b(b_idx)%last) then
                    sections_overlap = .true.
                    return
                end if
            end do
        end do
    end function sections_overlap

    pure logical function section_mask_overlap(sections, mask) result(overlap)
        type(array_index), dimension(:), intent(in) :: sections
        logical, dimension(:), intent(in) :: mask
        integer :: idx

        overlap = .false.
        do idx = 1, size(sections)
            if(any(mask(sections(idx)%first:sections(idx)%last))) then
                overlap = .true.
                return
            end if
        end do
    end function section_mask_overlap
end module mod_array_index
