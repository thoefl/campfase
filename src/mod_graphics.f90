! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_graphics
    
use mod_regular_grid
implicit none
private

public :: flood_fill_8, flood_fill_4, label_areas, label_nodes

contains

pure recursive subroutine flood_fill_8(grid, a, dim_idx, fill)
    type(regular_grid), intent(in) :: grid
    integer, intent(in) :: dim_idx(:)
    logical, contiguous, intent(in) :: a(:)
    logical, contiguous, intent(inout) :: fill(:)
    integer, dimension(grid%n_dims) :: offset, dim_idx_offset
    integer :: dim, linear_idx
    
    linear_idx = grid%dim_to_linear_idx(dim_idx)
    if(.not. a(linear_idx) .or. fill(linear_idx)) return
    fill(linear_idx) = .true.
    offset = -1
    do
        call grid%offset_dim_idx(dim_idx, offset, dim_idx_offset)
        call flood_fill_8(grid, a, dim_idx_offset, fill)
        if(all(offset == 1)) exit
        do dim = 1, grid%n_dims
            if(offset(dim) < 1) then
                offset(dim) = offset(dim) + 1
                exit
            else
                offset(dim) = -1
            end if
        end do
    end do
end subroutine flood_fill_8
    
pure recursive subroutine flood_fill_4(grid, a, dim_idx, fill)
    type(regular_grid), intent(in) :: grid
    integer, intent(in) :: dim_idx(:)
    logical, contiguous, intent(in) :: a(:)
    logical, contiguous, intent(inout) :: fill(:)
    integer, dimension(grid%n_dims) :: offset, dim_idx_offset
    integer :: dim, linear_idx, offset_val
    
    linear_idx = grid%dim_to_linear_idx(dim_idx)
    if(.not. a(linear_idx) .or. fill(linear_idx)) return
    fill(linear_idx) = .true.
    do dim = 1, grid%n_dims
        offset = 0
        do offset_val = -1, 1, 2
            offset(dim) = offset_val
            call grid%offset_dim_idx(dim_idx, offset, dim_idx_offset)
            call flood_fill_4(grid, a, dim_idx_offset, fill)
        end do
    end do
end subroutine flood_fill_4

subroutine label_areas(grid, a, labels, nlabels)
    type(regular_grid), intent(in) :: grid
    logical, intent(in) :: a(:)
    integer, intent(out) :: labels(:), nlabels
    integer :: curr_label, linear_idx, dim_idx(grid%n_dims)
    logical, allocatable :: fill_mask(:)
   
    allocate(fill_mask(size(a)))
    curr_label = 0
    labels = 0
    
    do linear_idx = 1, grid%n_nodes_total
        if(a(linear_idx) .and. labels(linear_idx) == 0) then
            fill_mask = .false.
            dim_idx = grid%linear_to_dim_idx(linear_idx)
            call flood_fill_4(grid, a, dim_idx, fill_mask)
            curr_label = curr_label + 1
            where(fill_mask) labels = curr_label
        end if
    end do
    nlabels = curr_label
end subroutine label_areas
    
subroutine label_nodes(grid, a, nodes, labels, nlabels, overlap_labels)
    type(regular_grid), intent(in) :: grid
    logical, intent(in) :: a(:)
    integer, intent(in) :: nodes(:)
    integer, intent(out) :: labels(:), overlap_labels(:), nlabels
    integer :: curr_label, linear_idx, dim_idx(grid%n_dims)
    logical, allocatable :: fill_mask(:)
    
    allocate(fill_mask(size(a)))
    curr_label = 0
    labels = 0
    overlap_labels = 0
    nlabels = size(nodes)
    
    do curr_label = 1, size(nodes)
        linear_idx = nodes(curr_label)
        if(labels(linear_idx) > 0) then
            overlap_labels(curr_label) = labels(linear_idx)
            nlabels = nlabels - 1
            cycle
        end if
        fill_mask = .false.
        dim_idx = grid%linear_to_dim_idx(linear_idx)
        call flood_fill_4(grid, a, dim_idx, fill_mask)
        where(fill_mask) labels = curr_label
    end do
end subroutine label_nodes
end module mod_graphics
