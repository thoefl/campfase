! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_fcdata
    use mod_base_functions
    use mod_boundary_condition
    use mod_regular_grid
    
    implicit none
    
    type fourier_cont_data
        real(wp), dimension(:), allocatable :: cont_bc_l, cont_bc_r
        real(wp), dimension(:,:), allocatable :: cont_mat_l, cont_mat_r
        integer :: d_l, d_r
    end type fourier_cont_data
    
    interface fourier_cont_data
        module procedure :: fourier_cont_data_constructor
    end interface fourier_cont_data
    
    contains
    
    pure type(fourier_cont_data) function&
    fourier_cont_data_constructor(cont_bc_l, cont_bc_r, cont_mat_l, cont_mat_r) result(fc_data)
        real(wp), intent(in), dimension(:) :: cont_bc_l, cont_bc_r
        real(wp), intent(in), dimension(:,:) :: cont_mat_r, cont_mat_l
        
        if(.not. (size(cont_bc_l,1) == size(cont_bc_r,1)&
            .and. size(cont_bc_r,1) == size(cont_mat_l,1)&
            .and. size(cont_mat_l,1) == size(cont_bc_r,1))) then
                error stop 'ERROR: fourier_cont_data_constructor: Dimension 1 of all data encoord_tries needs to be equal!'
        end if
        fc_data%cont_bc_l = cont_bc_l
        fc_data%cont_bc_r = cont_bc_r
        fc_data%cont_mat_l = cont_mat_l
        fc_data%cont_mat_r = cont_mat_r
        fc_data%d_l = size(cont_mat_l,2) + 1
        fc_data%d_r = size(cont_mat_r,2) + 1
    end function fourier_cont_data_constructor
    
    ! All functions are designed for the right side interpolation function
    
    pure elemental real(wp) function smoothing_func(x, w)
        real(wp), intent(in) :: x, w
        
        if(x <= -1) then
            smoothing_func = 1
        else if(x >= 1) then
            smoothing_func = 0
        else
            smoothing_func = 0.5_wp*(1 - erf(w*x/sqrt(1 - x**2)))
        end if
    end function smoothing_func
    
    ! M_idx denotes the index along the M-axis, d_idx for the d-axis, both starting at 1
    pure real(wp) elemental function neumann_factor(M_idx, d_idx, d)
        integer, intent(in) :: M_idx, d_idx, d
        
        if(d_idx < d) then
            neumann_factor = 1 - (1.0_wp/M_idx + 1.0_wp/(d - d_idx))/harmonic(d - 1)
        else if(d_idx == d) then
            neumann_factor = 1/harmonic(d - 1)
        else
            error stop 'ERROR: neumann_factor: index d_idx has to be <= d!'
        end if
    end function neumann_factor
    
    pure real(wp) elemental function cont_mat_builder(bc_type, M, M_idx, d_idx, d, w) result(cont)
        integer(int8), intent(in) :: bc_type
        integer, intent(in) :: M, M_idx, d_idx, d
        real(wp), intent(in) :: w
        integer(int64) :: base, factor
        integer :: i
        real(wp) :: smooth_x
        
        base = M_idx
        factor = (-1)**(d + 1)
        smooth_x = 2*real(M_idx - 1, wp)/(M - 1) - 1
        do i = 1, d - 1
            base = (base*(M_idx + i))//i
            if(base < 0) error stop 'ERROR: cont_mat_builder: overflow occured: base = ' // convert_to_char(base)
        end do
        do i = 1, d_idx - 1
            if(sign(int(1, int64), factor*i) /= sign(int(1, int64), factor)) then
                error stop 'ERROR: cont_mat_builder: overflow occured: factor = ' // convert_to_char(factor)
            end if
            factor = (factor*(i - d))//i
        end do
        if(bc_type == bc_dirichlet) then
            cont = (base*factor//(M_idx + d - d_idx))*smoothing_func(smooth_x, w)
        else if(bc_type == bc_neumann) then
            cont = (base*factor//(M_idx + d - d_idx))*smoothing_func(smooth_x, w)*neumann_factor(M_idx, d_idx, d)
        else
            error stop 'ERROR: cont_mat_builder: unrecognized boundary condition type!'
        end if
    end function cont_mat_builder
    
    subroutine set_fcdata(grid, fc_nodes, d_l, d_r, w, fc_data)
        integer, intent(in) :: d_l, d_r, fc_nodes
        type(regular_grid), intent(in) :: grid
        real(wp), intent(in) :: w
        type(fourier_cont_data), intent(out) :: fc_data
        real(wp) :: cont_mat_l(fc_nodes, d_l-1), cont_mat_r(fc_nodes, d_r-1), cont_bc_l(fc_nodes), cont_bc_r(fc_nodes)
        integer :: m_idx, d_idx
        
        do concurrent(m_idx = 1:fc_nodes, d_idx = 1:d_r-1)
            cont_mat_r(m_idx, d_idx) = cont_mat_builder(grid%axes(1)%bound_cond%type_max, fc_nodes, M_idx, d_idx, d_r, w)
        end do
        do concurrent(m_idx = 1:fc_nodes)
            cont_bc_r(m_idx) = cont_mat_builder(grid%axes(1)%bound_cond%type_max, fc_nodes, M_idx, d_r, d_r, w)
        end do
        if(grid%axes(1)%bound_cond%type_max == bc_neumann) cont_bc_r = cont_bc_r/(grid%axes(1)%n_nodes - fc_nodes - 1)
        do concurrent(m_idx = 1:fc_nodes, d_idx = 1:d_l-1)
            cont_mat_l(m_idx, d_idx) = cont_mat_builder(grid%axes(1)%bound_cond%type_min, fc_nodes, M_idx, d_idx, d_l, w)
        end do
        cont_mat_l = cont_mat_l(fc_nodes:1:-1,d_l-1:1:-1)
        do concurrent(m_idx = 1:fc_nodes)
            cont_bc_l(m_idx) = cont_mat_builder(grid%axes(1)%bound_cond%type_min, fc_nodes, M_idx, d_l, d_l, w)
        end do
        cont_bc_l = cont_bc_l(fc_nodes:1:-1)
        if(grid%axes(1)%bound_cond%type_min == bc_neumann) cont_bc_l = -cont_bc_l/(grid%axes(1)%n_nodes - fc_nodes - 1)
        fc_data = fourier_cont_data(cont_bc_l, cont_bc_r, cont_mat_l, cont_mat_r)
    end subroutine set_fcdata
end module mod_fcdata