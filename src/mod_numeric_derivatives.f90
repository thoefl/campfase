! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_numeric_derivatives
    
    use mod_mu
    use mod_pfm
    implicit none
    contains
    
    !function num_df_dc(c_field, phi_field, eta, sqsum, m_order, kappa, species_idx)
    !    type(concentration_field), intent(in) :: c_field(:)
    !    type(fft_field), intent(in) :: phi_field(:), eta
    !    real(wp), intent(in) :: m_order, kappa, sqsum(:)
    !    integer, intent(in) :: species_idx
    !    type(concentration_field), dimension(size(c_field)) :: c_plus, c_minus
    !    real(wp), dimension(size(c_field(1)%vals)) :: f_plus, f_minus, num_df_dc
    !    real(wp), parameter :: delta_num = 1E-8_wp
    !    
    !    c_plus = c_field
    !    c_minus = c_field
    !    c_plus(species_idx)%vals = c_plus(species_idx)%vals + delta_num
    !    c_minus(species_idx)%vals = c_minus(species_idx)%vals - delta_num
    !    f_plus = f_total(c_plus, phi_field, eta, sqsum, m_order, kappa, fft)
    !    f_minus = f_total(c_minus, phi_field, eta, sqsum, m_order, kappa, fft)
    !    num_df_dc = (f_plus - f_minus)/(2*delta_num)
    !end function num_df_dc
    !
    !function num_mu_c(c_field, phi_field, eta, sqsum, m_order, kappa, species_idx)
    !    type(concentration_field), intent(in) :: c_field(:)
    !    type(fft_field), intent(in) :: phi_field(:), eta
    !    real(wp), intent(in) :: m_order, kappa, sqsum(:)
    !    integer, intent(in) :: species_idx
    !    real(wp), dimension(size(c_field(1)%vals)) :: df_dc, num_mu_c_x, num_mu_c_y&
    !        , num_mu_c, c_phase, c_phase_x, c_phase_y, num_mu_c_xx, num_mu_c_yy
    !    real(wp), parameter :: delta_num = 1E-8_wp
    !    integer :: i
    !    
    !    num_mu_c_x = 0
    !    num_mu_c_y = 0
    !    df_dc = num_df_dc(c_field, phi_field, eta, sqsum, m_order, kappa, fft, species_idx)
    !    do i = 1, n_phases
    !        c_phase = df_dc/(2*c_field(species_idx)%p_coeff(i))&
    !            + c_field(species_idx)%equi_val(i) - c_field(species_idx)%k_coeff(i)/(2*c_field(species_idx)%p_coeff(i))
    !        !c_phase = get_c_phase(phi_field, c_field(species_idx), i)
    !        call fft%grad_x(get_phi(phi_field, i), c_phase_x)
    !        call fft%grad_y(get_phi(phi_field, i), c_phase_y)
    !        num_mu_c_x = num_mu_c_x - c_phase*c_phase_x
    !        num_mu_c_y = num_mu_c_y - c_phase*c_phase_y
    !    end do
    !    call fft%grad_x(num_mu_c_x, num_mu_c_xx)
    !    call fft%grad_y(num_mu_c_y, num_mu_c_yy)
    !    num_mu_c = num_mu_c_xx + num_mu_c_yy
    !end function num_mu_c
    !   
    !function num_df_chem_deta(c_field, eta_start, phase_idx, diff_idx)
    !    type(concentration_field), intent(in) :: c_field(:)
    !    type(fft_field) :: phi_field(n_phases - 1)
    !    real(wp), dimension(:), intent(in) :: eta_start
    !    integer, intent(in) :: phase_idx, diff_idx
    !    real(wp), dimension(size(c_field(1)%vals)) :: f_plus, f_minus, num_df_chem_deta, sqsum, eta_square, eta
    !    real(wp), parameter :: d_eta = 1E-8_wp
    !    
    !    if(phase_idx == 0) then
    !        eta = 0.0_wp
    !    else if(this_image() == diff_idx) then
    !        eta = eta_start + d_eta
    !    else
    !        eta = eta_start
    !    end if
    !    call phi_calc(eta, eta_square, sqsum, phi_field)
    !    call field_broadcast(phi_field)
    !    
    !    f_plus = f_chem(c_field, phi_field)
    !    if(this_image() == diff_idx) then
    !        eta = eta_start - d_eta
    !    end if
    !    call phi_calc(eta, eta_square, sqsum, phi_field)
    !    call field_broadcast(phi_field)
    !    f_minus = f_chem(c_field, phi_field)
    !    num_df_chem_deta = (f_plus - f_minus)/(2*d_eta)
    !end function num_df_chem_deta
    !
    !function num_df_order_deta(c_field, eta_start, m_order, phase_idx, diff_idx)
    !    type(concentration_field), intent(in) :: c_field(:)
    !    type(fft_field) :: phi_field(n_phases - 1)
    !    real(wp), dimension(:), intent(in) :: eta_start
    !    real(wp), intent(in) :: m_order
    !    integer, intent(in) :: phase_idx, diff_idx
    !    real(wp), dimension(size(c_field(1)%vals)) :: f_plus, f_minus, num_df_order_deta, sqsum, quadsum, eta_square, eta
    !    real(wp), parameter :: d_eta = 1E-8_wp
    !    
    !    if(phase_idx == 0) then
    !        eta = 0.0_wp
    !    else if(this_image() == diff_idx) then
    !        eta = eta_start + d_eta
    !    else
    !        eta = eta_start
    !    end if
    !    call phi_calc(eta, eta_square, sqsum, phi_field)
    !    call field_broadcast(phi_field)
    !    quadsum = eta**4
    !    call co_sum(quadsum)
    !    
    !    f_plus = f_order(sqsum, quadsum, m_order)
    !    if(this_image() == diff_idx) then
    !        eta = eta_start - d_eta
    !    end if
    !    call phi_calc(eta, eta_square, sqsum, phi_field)
    !    call field_broadcast(phi_field)
    !    quadsum = eta**4
    !    call co_sum(quadsum)
    !    f_minus = f_order(sqsum, quadsum, m_order)
    !    num_df_order_deta = (f_plus - f_minus)/(2*d_eta)
    !end function num_df_order_deta
    !
    !function num_df_deta_nograd(c_field, eta_start, m_order, phase_idx, diff_idx)
    !    type(concentration_field), intent(in) :: c_field(:)
    !    real(wp), dimension(:), intent(in) :: eta_start
    !    real(wp), intent(in) :: m_order
    !    integer, intent(in) :: phase_idx, diff_idx
    !    real(wp), dimension(size(c_field(1)%vals)) :: num_df_deta_nograd
    !    
    !    num_df_deta_nograd = num_df_chem_deta(c_field, eta_start, phase_idx, diff_idx)&
    !        + num_df_order_deta(c_field, eta_start, m_order, phase_idx, diff_idx)
    !end function num_df_deta_nograd
end module mod_numeric_derivatives