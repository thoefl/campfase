! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

!TODO: include T/c dependence of molar volume, thermal expansion, ...

module mod_parameters

  use, intrinsic :: iso_fortran_env, only: int8
  use mod_realtype
  use mod_timer_set
  use mod_piecewise_linear
  implicit none

  type grain_reindex_performance_setting
    logical :: enabled
    integer :: n_eta, n_eta_min_active
    real(wp) :: fill_ratio
  end type

  type grain_separators
    real(wp) :: intra, inter, fill, rel_tol
  end type grain_separators

  type grain_reindex_settings
    type(timer_set) :: timer
    type(grain_reindex_performance_setting) :: for_performance
    type(grain_separators) :: separators
  end type

  type chemical_energy_model
    real(wp) :: p_coeff_min, p_coeff_max
  end type chemical_energy_model

  type restore_file_settings
    logical :: enabled, c_only
    integer :: idx
    character(:), allocatable :: hdf_filename
  end type restore_file_settings

  type benchmark_settings
    logical :: enabled
    integer :: count
  end type benchmark_settings

  type check_settings
    logical :: equilibrium, energy_model
    integer :: debug_level
  end type check_settings

  type fourier_continuation_settings
    integer :: n_nodes
    real(wp) :: w
    logical :: lazy_boundary
  end type fourier_continuation_settings

  type opencalphad_settings
    character(:), allocatable :: tdb_filename
  end type opencalphad_settings

  !! f0 = DG*FU/V_m
  !! p = omega/2*V_m/FU
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Scaling factors
  !real(wp), parameter :: V_m_factor = 1E-5_wp                     ! [m^3/mol_FU]
  !real(wp), parameter :: G_m_factor = 1E3_wp                      ! [J/mol_at]
  !real(wp), parameter :: D_factor = 1.0E-14_wp                    ! [m^2/s]
  !
  !real(wp), parameter :: meter = 1E-7_wp
  !real(wp), parameter :: second = meter**2/D_factor
  !real(wp), parameter :: mole_FU = meter**3/V_m_factor
  !real(wp), parameter :: mole = mole_FU
  !real(wp), parameter :: joule = G_m_factor*mole
  !real(wp), parameter :: kelvin = 1.0_wp
  !
  !real(wp), parameter :: units%t = second                          ! [s]
  !real(wp), parameter :: if_energy_unit = joule/meter**2          ! [J/m^2]
  !real(wp), parameter :: m_order_unit = joule/meter**3            ! [J/m^3]
  !real(wp), parameter :: kappa_unit = joule/meter                 ! [J/m]
  !real(wp), parameter :: M_eta_unit = meter**3/(joule*second)     ! [m^3/(J*s)]
  !real(wp), parameter :: units%df_chem_deta = joule/meter**3       ! [J/m^3]
  !real(wp), parameter :: units%mu_c = mole/meter**5                ! [mol/m^5]
  !real(wp), parameter :: units%k_coeff = joule/mole                ! [J/mol]
  !real(wp), parameter :: units%p_coeff = joule*meter**3/mole**2    ! [J*m^3/mol^2]
  !real(wp), parameter :: chem_pot_unit = joule/mole               ! [J/mol]
  !real(wp), parameter :: units%omega = joule/mole                  ! [J/mol]
  !real(wp), parameter :: c_unit = mole/meter**3                   ! [mol/m^3]
  !real(wp), parameter :: R_gas_unit = Joule/mole/Kelvin           ! [J/mol/K]
  !real(wp), parameter :: N_A_unit = 1/mole                        ! [1/mole]
  !
  !real(wp), parameter :: minute = 60*second
  !real(wp), parameter :: hour = 60*minute
  !real(wp), parameter :: newton = joule/meter
  !real(wp), parameter :: kilogram = newton/meter*second**2
  !real(wp), parameter :: gram = kilogram/1000
  !real(wp), parameter :: pascal = newton/meter**2

  real(wp), parameter :: stab_coeff_bdfab = 1.0_wp
  real(wp), parameter :: stab_coeff_be = 0.0_wp

  real(wp), parameter :: eta_min = -1.1_wp
  real(wp), parameter :: eta_max = 1.1_wp
  real(wp), parameter :: c_min = -1.0_wp
  real(wp), parameter :: c_max = 100.0_wp
  real(wp), parameter :: phi_max = 1.0_wp
  real(wp), parameter :: phi_min = 0.0_wp

  integer, parameter :: hdf_max_saves = 1E6

  !real(wp), parameter :: oscillation_lim_factor = 1/dt_incmax
  !real(wp), parameter :: dt_incmax_osc = (1/oscillation_lim_factor)**(2.0_wp*real(ndt,wp)/oscillation_lim_adt)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Other physical properties
  !!FIXME
  !integer(int8), parameter :: phase_mobility_model = manual_mobility_model
  !!FIXME
  !real(wp), dimension(n_phases), parameter :: M_factor = [1.0E5_wp, 1.0E5_wp, 1.0E5_wp]*D_factor/M_eta_unit
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Gradient energy coefficient parameters
  !! if_width:     Target interface width, based upon which kappa is calculated
  !!               dependent on temperature
  !! delta_num_GT ~ sigma*V_m_factor/(2*p)
  !!sigma*V_m/(r_nbc*df*(c_nbc - c_sol))
  !
  !real(wp), parameter :: if_nodes = 5.0_wp
  !real(wp), parameter :: if_energy_setval = 1.0_wp/if_energy_unit
  !real(wp), parameter :: if_energy_settemp = max(temp_hold, temp_start)      ![K]
  !real(wp), parameter :: gamma = 1.5_wp
  !real(wp), parameter :: p_coeff_min = 1.0E-3_wp/G_m_factor
  !real(wp), parameter :: p_coeff_max = 2.5E4_wp/G_m_factor
  !real(wp), parameter :: k_coeff_max = 1.0E16_wp/chem_pot_unit
  !real(wp), parameter :: k_coeff_min = -1.0E16_wp/chem_pot_unit
  !real(wp), parameter :: p_coeff_max_ratio_species = 1.0E10_wp
  !real(wp), parameter :: p_coeff_max_ratio_total = 0.0_wp
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Nucleation parameters
  !! Nucleation is implemented by imposing a penalty to the chemical potential,
  !! thus causing a gradient and encouraging diffusion into the nucleus.
  !! Nucleus stability is mainly influenced by supersaturation dC, nucleus size
  !! and the gradient energy coefficient (lower values cause lower critical size).
  !! The free energy penalty is shaped by a tanh function, approximating
  !! equlibrium.
  !!
  !! Nucleation probability is calculated according to equation
  !! P = 1 - exp(-k_1*exp(-k_2/dC)*dx*dy*dt)
  !!
  !! k_1: Preexponential factor, defines maximum nucleation rate
  !! k_2: Exponential factor, defines rate dependence on oversaturation dC
  !! nucl_size:      Initial size for the nucleus.
  !! nucl_penalty:   Maximum of mu penalty enforced causing nucleation
  !! nucl_target:    Used as a criterion when to lift enforced mu penalty
  !! nucl_initseeds: Number of nucleation nodes in random positions seeded at
  !!                 the beginning of the simulation. Seeding should be
  !!                 deterministic for given number of grid%n_nodes_total nodes.
  !real(wp), parameter :: k_1 = -5.0E12_wp
  !real(wp), parameter :: k_2 = 0.2_wp
  !real(wp), parameter :: nucl_size = 4.0_wp*dx
  !real(wp), parameter :: equilib_size = 2.57E-6_wp
  !real(wp), parameter :: nucl_target = max(dt_errtol, 1E-6_wp)
  !integer, parameter :: nucl_initseeds = 16
  !logical, parameter :: nucl_init_random = .true.
  !integer, parameter :: nucl_max_steps = 1
  !integer, parameter :: voronoi_attempts = 5
  !real(wp), parameter :: voronoi_eta_min_dist = 1.2_wp*eta_dist_lim
  !real(wp), parameter :: voronoi_rel_density = 0.95_wp
  !real(wp), parameter :: voronoi_initial_fill_ratio = 0.1_wp
  !real(wp), parameter :: voronoi_grain_r_sigma = 0.3_wp
  !real(wp), parameter :: voronoi_grain_r_median = 2.5E-6_wp/meter
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Derived parameters
  !integer, parameter :: eta_team_num = n_species + 1
  !integer, parameter :: c_team_num = n_phases + 1
  !
  !integer, parameter :: grid%n_nodes_total = nx*ny
  !integer(int16), parameter :: x(nx) = [(i, i=0, nx-1)]
  !integer(int16), parameter :: y(ny) = [(i, i=0, ny-1)]
  !integer, parameter :: nx_fft = nx/2 + 1
  !integer, parameter :: ny_fft = ny
  !integer, parameter :: grid%n_nodes_total_fft = nx_fft*ny_fft
  !real(wp), parameter :: fft_dx = 2*pi/nx
  !real(wp), parameter :: fft_dy = 2*pi/ny
  !real(wp), parameter :: fft_kx(nx_fft) = [(i*fft_dx, i=0, nx/2)]
  !real(wp), parameter :: fft_ky(ny_fft) = cshift([(i*fft_dy, i=-max(ny/2, 1) + 1, ny/2)], ny/2 - 1)
  !real(wp), parameter :: fft_kx_zero(nx_fft) = [(i*fft_dx, i=0, nx/2 - 1), 0.0_wp]
  !real(wp), parameter :: fft_ky_zero(ny_fft) = cshift([(i*fft_dy, i=-max(ny/2, 1) + 1, ny/2 - 1), 0.0_wp], ny/2 - 1)
  !real(wp), parameter :: fft_kx_minus(nx_fft) = [(i*fft_dx, i=0, nx/2 - 1), -nx/2*fft_dx]
  !real(wp), parameter :: fft_ky_minus(ny_fft) = cshift([(i*fft_dy, i=-max(ny/2, 1) + 1, ny/2 - 1), -ny/2*fft_dy], ny/2 - 1)
  !integer, parameter :: fft_maxfreqnode_xy = nx_fft*(1 + ny_fft/2)
  !integer, parameter :: fft_maxfreqnode_x = nx_fft
  !integer, parameter :: fft_maxfreqnode_y = nx_fft*ny_fft/2 + 1
  !integer, parameter :: n_osc_nodes = 3
  !integer, parameter :: fft_osc_check_nodes(n_osc_nodes) = [fft_maxfreqnode_x, fft_maxfreqnode_y, fft_maxfreqnode_xy]
  !character(2), parameter :: fft_osc_check_names(n_osc_nodes) = [character(2) :: 'x', 'y', 'xy']
  !real(wp), parameter :: fft_osc_check_errtol(n_osc_nodes) = [fft_osc_errtol_x, fft_osc_errtol_y, fft_osc_errtol_xy]
  !! Could be initialized here, but very slow for large nx*ny
  !!real(wp), parameter :: fft_k2(grid%n_nodes_total_fft) = [((i + j, i=1, nx_fft), j=1, ny_fft)]
  !!real(wp), parameter :: fft_k4(grid%n_nodes_total_fft) = fft_k2**2
  !

  integer(int8), parameter :: init_monocrystal = 0
  integer(int8), parameter :: init_polycrystal = 1
  integer(int8), parameter :: init_single_boundary_pin = 2
  integer(int8), parameter :: init_boundary_check = 3
  integer(int8), parameter :: init_spheres = 4
  integer(int8), parameter :: init_planar_x = 5
  integer(int8), parameter :: init_planar_y = 6
  integer(int8), parameter :: init_planar_xy = 7
  integer(int8), parameter :: init_quad_boundary_pin = 8
  integer(int8), parameter :: init_spinodal = 9
  integer(int8), parameter :: init_spinodal_dense = 10
  integer(int8), parameter :: init_sine = 11

  integer(int8), parameter :: init_full_equilib = -1
  integer(int8), parameter :: init_base_equilib = -2
  integer(int8), parameter :: init_base_flat = -3
  integer(int8), parameter :: init_equilib_flat = -4
  integer(int8), parameter :: init_full_equilib_free_eta = -5
  integer(int8), parameter :: init_dissolve_test = -6
  integer(int8), parameter :: init_defined_conc = -7

  integer(int8), parameter :: manual_energy = 0
  integer(int8), parameter :: crit_fit_energy = 1
  integer(int8), parameter :: parabolic_fit_energy = 2
  integer(int8), parameter :: centered_barrier_energy = 3
  integer(int8), parameter :: parabolic_fit_solonly_energy = 4

  integer(int8), parameter :: manual_mobility_model = 0
  integer(int8), parameter :: diff_limit_mobility_model = 1

  character(*), parameter :: div_line = '================================================================================'

  character(*), parameter :: mu_name = 'mu'
  character(*), parameter :: c_name = 'c'
  character(*), parameter :: eta_name = 'eta'
  character(*), parameter :: phi_name = 'phi'
  character(*), parameter :: c_equi_name = 'c_equi'
  character(*), parameter :: f0_name = 'f_0'
  character(*), parameter :: p_coeff_name = 'p_coeff'
  character(*), parameter :: k_coeff_name = 'k_coeff'
  character(*), parameter :: omega_name = 'omega'
  character(*), parameter :: x_name = 'x'
  character(*), parameter :: y_name = 'y'
  character(*), parameter :: t_now_name = 'time'
  character(*), parameter :: dt_name = 'dt'
  character(*), parameter :: step_name = 'step'
  character(*), parameter :: backup_step_name = 'step_backup'
  character(*), parameter :: n_t_error_name = 'n_t_error'
  character(*), parameter :: c_error_name = 'c_error'
  character(*), parameter :: eta_error_name = 'eta_error'
  character(*), parameter :: fft_maxfreq_x_name = 'fft_maxfreq_x'
  character(*), parameter :: fft_maxfreq_y_name = 'fft_maxfreq_y'
  character(*), parameter :: fft_maxfreq_xy_name = 'fft_maxfreq_xy'
  character(*), parameter :: temp_now_name = 'temperature'
  character(*), parameter :: M_name = 'mobility'
  character(*), parameter :: D_name = 'diff_coeff'
  character(*), parameter :: t_cpu_name = 'cpu time'
  character(*), parameter :: t_wall_name = 'wall time'
  character(*), parameter :: kappa_name = 'kappa'
  character(*), parameter :: f_total_name = 'f_total'
  character(*), parameter :: f_total_sum_name = 'f_total_sum'
  character(*), parameter :: nx_name = 'nx'
  character(*), parameter :: ny_name = 'ny'
  character(*), parameter :: dx_name = 'dx'
  character(*), parameter :: dy_name = 'dy'
  character(*), parameter :: V_m_name = 'molar_volume'
  character(*), parameter :: G_m_factor_name = 'G_m_factor'
  character(*), parameter :: D_factor_name = 'diff_coeff_factor'
  character(*), parameter :: nucl_k1_name = 'nucl_k1'
  character(*), parameter :: nucl_k2_name = 'nucl_k2'
  character(*), parameter :: nucl_size_name = 'nucl_size'
  character(*), parameter :: if_width_name = 'if_width'
  character(*), parameter :: if_energy_name = 'if_sigma'
  character(*), parameter :: if_nodes_name = 'if_nodes'
  character(*), parameter :: nucl_penalty_name = 'nucl_penalty'
  character(*), parameter :: nucl_target_name = 'nucl_target'
  character(*), parameter :: nucl_node_name = 'nucl_node'
  character(*), parameter :: randseed_name = 'random seed'
  character(*), parameter :: dt_errtol_name = 'dt_errtol'
  character(*), parameter :: nucl_dt_name = 'dt_nucleation'
  character(*), parameter :: m_order_name = 'm_order'
  character(*), parameter :: eta_idx_name = 'eta_idx'
  character(*), parameter :: mu_eta_name = 'mu_eta'
  character(*), parameter :: mu_eta_common_name = 'mu_eta_common'
  character(*), parameter :: M_eta_name = 'M_eta'
  character(*), parameter :: nucl_attempts_name = 'nucl_attempts'
  character(*), parameter :: p_nucl_name = 'p_nucl'
  character(*), parameter :: sqsum_name = 'eta_square_sum'
  character(*), parameter :: eta_error_field_name = 'eta_error_field'
  character(*), parameter :: c_error_field_name = 'c_error_field'
  character(*), parameter :: c_ni_name = 'c_ni'
  character(*), parameter :: n_eta_limit_name = 'n_eta_limit'
  character(*), parameter :: n_c_limit_name = 'n_c_limit'
  character(*), parameter :: n_c_ni_limit_name = 'n_c_ni_limit'
  character(*), parameter :: n_osc_limit_name = 'n_osc_limit'
  character(*), parameter :: osc_error_name = 'osc_error'
  character(*), parameter :: osc_maxfreq_name = 'osc_maxfreq'
  character(*), parameter :: eta_sol_name = 'eta_sol'
  character(*), parameter :: eta_nbc_name = 'eta_nbc'
  character(*), parameter :: eta_pore_name = 'eta_pore'
  character(*), parameter :: use_bdfab_name = 'use_bdfab'
  character(*), parameter :: prestep_c_name = 'prestep_c'
  character(*), parameter :: stab_coeff_bdfab_name = 'stabilisation_coeff_bdfab'
  character(*), parameter :: stab_coeff_be_name = 'stabilisation_coeff_be'
  character(*), parameter :: surf_x_node_name = 'surface_x_node'
  character(*), parameter :: error_name = 'error'
  character(*), parameter :: error_field_name = 'error_field'
  character(*), parameter :: n_limit_name = 'n_limit'
  character(*), parameter :: common_name = 'common'
  character(*), parameter :: f_grad_name = 'f_gradient'
  character(*), parameter :: c_mix_name = 'c_mix'
  character(*), parameter :: p_mix_name = 'p_mix'
  character(*), parameter :: log_p_mix_name = 'log_p_mix'
  character(*), parameter :: f0_mix_name = 'f0_mix'
  character(*), parameter :: df_chem_deta_common_name = 'df_chem_deta_common'
  character(*), parameter :: df_order_deta_name = 'df_order_deta'

  character(*), parameter :: valid_filename_chars = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'

  character(*), parameter :: json_extension = '.json'

  integer, parameter :: eta_team_num = huge(0) - 1
  integer, parameter :: c_team_num = eta_team_num - 1
  integer, parameter :: primary_team_num = eta_team_num - 2
  integer, parameter :: secondary_team_num = eta_team_num - 3
end module mod_parameters
