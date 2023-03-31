! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_nucleus

use mod_base_functions
use mod_particle
use mod_regular_grid
use mod_timer_set
implicit none
private

type, extends(spheric_particle), public :: nucleus
    type(timer_set) :: expiration
    type(array_index), allocatable :: sections(:)
contains
    procedure, public :: nucleus_init
    procedure, public :: apply_to_mask
    procedure, public :: remove_from_mask
    procedure, public :: fits_in_mask
end type nucleus

contains
    
subroutine nucleus_init(this, grid, id, radius, coordinates, expiration_template, formation_step)
    type(regular_grid), intent(in), target :: grid
    class(nucleus), intent(inout) :: this
    integer, intent(in) :: id, coordinates(:)
    integer(int64), intent(in) :: formation_step
    real(wp), intent(in) :: radius
    type(timer_set), intent(in) :: expiration_template
        
    call this%init(grid, id, radius, coordinates)
    this%sections = grid%get_circle_sections(coordinates, radius)
    this%expiration = expiration_template
    call this%expiration%new_schedule(formation_step)
    call this%expiration%schedule(formation_step)
end subroutine nucleus_init

pure logical function fits_in_mask(this, mask)
    class(nucleus), intent(in) :: this
    integer, intent(in) :: mask(:)
    integer :: i
    
    fits_in_mask = .true.
    associate(idx => this%sections)
        do i = 1, size(idx)
            if(any(mask(idx(i)%first:idx(i)%last) > 0)) then
                fits_in_mask = .false.
                return
            end if
        end do
    end associate
end function fits_in_mask

pure subroutine apply_to_mask(this, mask)
    class(nucleus), intent(in) :: this
    integer, intent(inout) :: mask(:)
    integer :: i
    character(*), parameter :: proc_name = 'apply_to_mask'
    
    associate(idx => this%sections)
        do concurrent(i = 1:size(idx))
            if(any(mask(idx(i)%first:idx(i)%last) > 0)) then
                call error_msg(proc_name, 'Mask already set in array section ' &
                    // convert_to_char(idx(i)%first) // ':' &
                    // convert_to_char(idx(i)%last))
            end if
            mask(idx(i)%first:idx(i)%last) = this%id
        end do
    end associate
end subroutine apply_to_mask

pure subroutine remove_from_mask(this, mask)
    class(nucleus), intent(in) :: this
    integer, intent(inout) :: mask(:)
    integer :: i
    
    associate(idx => this%sections)
        do concurrent(i = 1:size(idx))
            mask(idx(i)%first:idx(i)%last) = 0
        end do
    end associate
end subroutine remove_from_mask
    
!! Simple model for calculating nucleation probability in 2D
!! Source: Simmons et al., 2000: PFM of simultaneous nucleation and growth...
!!         https://doi.org/10.1016/S1359-6462(00)00517-0
!pure real(wp) elemental function &
!    nucl_probability(dC, dt_nucl, k_1, k_2, dx, dy) result(P)
!    real(wp), intent(in) :: dC, dt_nucl, k_1, k_2, dx, dy
!    
!    P = 1 - exp(-k_1*exp(-k_2/dC)*dx*dy*dt_nucl)
!end function nucl_probability
!
!subroutine nucl_add(node, phase_idx, nucl_node, nucl_phase_idx, nucl_sites, nucl_counter)
!    integer, intent(in) :: node, phase_idx
!    integer, intent(inout), allocatable, dimension(:) :: nucl_node, nucl_counter, nucl_phase_idx
!    integer, intent(inout) :: nucl_sites
!    
!    nucl_sites = nucl_sites + 1
!    nucl_node = [nucl_node, node]
!    nucl_phase_idx = [nucl_phase_idx, phase_idx]
!    nucl_counter = [nucl_counter, 0]
!end subroutine nucl_add
!
!subroutine nucl_del(node, nucl_node, nucl_phase_idx, nucl_sites, nucl_counter)
!    integer, intent(in) :: node
!    integer, intent(inout), allocatable, dimension(:) :: nucl_node, nucl_counter, nucl_phase_idx
!    integer, intent(inout) :: nucl_sites
!    
!    nucl_sites = nucl_sites - 1
!    nucl_node = [nucl_node(:node-1), nucl_node(node+1:)]
!    nucl_phase_idx = [nucl_phase_idx(:node-1), nucl_phase_idx(node+1:)]
!    nucl_counter = [nucl_counter(:node-1), nucl_counter(node+1:)]
!end subroutine nucl_del
!
!subroutine nucl_setmask(grid, nucl_node, nucl_phase_idx, phase_idx, nucl_size, nucl_mask, if_width)
!    type(regular_grid), intent(in) :: grid
!    integer, intent(in) :: nucl_node(:), nucl_phase_idx(:), phase_idx
!    real(wp), intent(in) :: nucl_size, if_width
!    real(wp), intent(out) :: nucl_mask(:)
!    integer, allocatable :: new_nucl_nodes(:)
!    integer :: i
!
!    nucl_mask = 0
!    do i = 1, size(nucl_node)
!        if(nucl_phase_idx(i) == phase_idx) then
!            if(verbosity_level >= 3) call print_message('Adding nucleus to mask, node = ' // convert_to_char(nucl_node(i)))
!            new_nucl_nodes = grid%get_circle_nodes()
!            nucl_mask = nucl_mask + &
!                nucl_circle(nucl_size, get_x(nucl_node(i)),&
!                get_y(nucl_node(i)), bc_a, if_width)
!        end if
!    end do
!    call output_image(nucl_mask, 'nucl_mask', nx, co_img_opt = .true.)
!end subroutine nucl_setmask
!
!pure logical function eligible_nucleus(c_field, phi_field, node, nucl_node, sim_mode, fc_nodes)
!    type(fft_field), intent(in) :: c_field(:)
!    type(fft_field), intent(in) :: phi_field(:)
!    integer, intent(in) :: nucl_node(:), node, fc_nodes
!    integer(int8), intent(in) :: sim_mode
!    real(wp) :: min_x, max_x, nucl_mindist, c_nearmax, &
!                nucl_nearestdist, existing_nucl_x(size(nucl_node)), &
!                existing_nucl_y(size(nucl_node)), node_x, node_y
!
!    eligible_nucleus = .true.
!    
!    !if(rho(node) < nucl_rho_lbound) then
!    !    eligible_nucleus = .false.
!    !    return
!    !end if
!    !
!    !if(c(node) > nucl_c_ubound) then
!    !    eligible_nucleus = .false.
!    !    return
!    !end if
!    !
!    !node_x = get_x(node)
!    !node_y = get_y(node)
!    !!if(bc_a%type_min == bc_periodic) then
!    !    min_x = x(1)
!    !    max_x = x(nx)
!    !!else
!    !!    min_x = x(1) + 2*nucl_size + if_width
!    !!    if(sim_mode == sim_fourier) then            
!    !!        max_x = x(nx - fc_nodes) - 2*(nucl_size + if_width)
!    !!    else
!    !!        max_x = x(nx) - 2*(nucl_size + if_width)
!    !!    end if
!    !!end if
!    !
!    !if(node_x < min_x .or. node_x > max_x) then
!    !    eligible_nucleus = .false.
!    !    return
!    !end if
!    !
!    nucl_mindist = nucl_size + if_nodes + sqrt(2.0_wp)
!    existing_nucl_x = get_x(nucl_node(:))
!    existing_nucl_y = get_y(nucl_node(:))
!    nucl_nearestdist = minval(point_distance(&
!        existing_nucl_x, existing_nucl_y, node_x, node_y, c_field(1)%bc))
!    
!    if(nucl_nearestdist <= nucl_mindist) then
!        eligible_nucleus = .false.
!        return
!    end if
!    !
!    !c_nearmax = maxval(pack(&
!    !    c, grid_distance(node_x, node_y, bc_a) <&
!    !        nucl_mindist))
!    !
!    !if(c_nearmax > 1.01_wp*c_matrix) then
!    !    eligible_nucleus = .false.
!    !    return
!    !end if
!end function eligible_nucleus
    
end module mod_nucleus
