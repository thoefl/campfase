! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_remap_data
    use mod_realtype
    use mod_array_index
    implicit none
    private
    public :: broadcast_remap_data

    type :: remap_target
        logical :: eligible
        real(wp) :: badness
    end type remap_target

    type, public :: remap_grain
        integer :: n_targets
        type(remap_target), allocatable, dimension(:) :: targets
        integer :: node, target_idx, n_conflicts
        logical :: merge_risk
        type(array_index), allocatable :: sections(:)
        real(wp) :: separator
    contains
        procedure, public :: select_target
    end type remap_grain

    type, public :: remap_data
        integer :: n_grains
        type(remap_grain), allocatable, dimension(:) :: grains
    contains
        procedure, public :: get_n
        procedure, public :: get_grain_map
    end type remap_data
   

    interface remap_data
        module procedure :: remap_data_constructor
    end interface remap_data

    interface remap_grain
        module procedure :: remap_grain_constructor
    end interface remap_grain
    
    contains

    pure elemental type(remap_data) function remap_data_constructor(n_grains, n_targets) result(this)
        integer, intent(in) :: n_grains, n_targets
       
        this%n_grains = n_grains
        allocate(this%grains(n_grains))
        this%grains = remap_grain(n_targets)
    end function remap_data_constructor

    pure elemental type(remap_grain) function remap_grain_constructor(n_targets) result(this)
        integer, intent(in) :: n_targets

        this%n_targets = n_targets
        allocate(this%targets(n_targets))
        this%targets%eligible = .false.
        this%targets%badness = 0.0_wp
        this%node = 0
        this%separator = 0.0_wp
        this%merge_risk = .false.
        this%n_conflicts = 0
        this%target_idx = 0
    end function remap_grain_constructor

    pure elemental subroutine select_target(this)
        class(remap_grain), intent(inout) :: this
        
        this%target_idx = minloc(this%targets%badness, 1, this%targets%eligible)
    end subroutine select_target    
 
    pure integer elemental function get_n(this) result(n)
        class(remap_data), intent(in) :: this
        
        n = count(this%grains%node > 0)
    end function get_n

    pure function get_grain_map(this, n_nodes) result(grain_map)
        class(remap_data), intent(in) :: this
        integer, intent(in) :: n_nodes
        integer, dimension(n_nodes) :: grain_map       
        integer :: grain_idx, section_idx
        
        grain_map = 0
        do concurrent(grain_idx = 1:this%n_grains, this%grains(grain_idx)%node > 0)
            do concurrent(section_idx = 1:size(this%grains(grain_idx)%sections))
                grain_map(this%grains(grain_idx)%sections(section_idx)%first:this%grains(grain_idx)%sections(section_idx)%last) = &
                grain_map(this%grains(grain_idx)%sections(section_idx)%first:this%grains(grain_idx)%sections(section_idx)%last) &
                    + grain_idx
            end do
        end do
    end function get_grain_map

    subroutine broadcast_remap_data(this)
        type(remap_data), intent(inout) :: this(:)
        integer :: img, grain_idx, co_size
        do img = 1, size(this)
            if(this_image() == img) then
                co_size = size(this(img)%grains)
            else
                if(allocated(this(img)%grains)) deallocate(this(img)%grains)
            end if
            call co_broadcast(co_size, img)
            if(this_image() /= img) allocate(this(img)%grains(co_size))

            do grain_idx = 1, size(this(img)%grains)
                if(this_image() == img) then
                    co_size = size(this(img)%grains(grain_idx)%sections)
                else
                    if(allocated(this(img)%grains(grain_idx)%sections)) deallocate(this(img)%grains(grain_idx)%sections)
                end if
                call co_broadcast(co_size, img)
                if(this_image() /= img) then
                    allocate(this(img)%grains(grain_idx)%sections(co_size))
                    allocate(this(img)%grains(grain_idx)%targets(num_images()))
                end if
                
                call co_broadcast(this(img)%grains(grain_idx)%sections%first, img)
                call co_broadcast(this(img)%grains(grain_idx)%sections%last, img)
                call co_broadcast(this(img)%grains(grain_idx)%targets%eligible, img)
                call co_broadcast(this(img)%grains(grain_idx)%targets%badness, img)
                call co_broadcast(this(img)%grains(grain_idx)%node, img)
                call co_broadcast(this(img)%grains(grain_idx)%n_conflicts, img)
                call co_broadcast(this(img)%grains(grain_idx)%target_idx, img)
                call co_broadcast(this(img)%grains(grain_idx)%merge_risk, img)
            end do

        end do
    end subroutine broadcast_remap_data
end module mod_remap_data
