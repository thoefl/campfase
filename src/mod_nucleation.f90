! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_nucleation
  use mod_base_functions
  use mod_regular_grid
  use mod_nucleus
  use mod_spectral
  use mod_mu
  use mod_order_parameter_set
  use mod_timer_set
  use mod_phase_properties
  use mod_constants
  use mod_phase_set
  use mod_species_properties
  use mod_random
  use mod_arrhenius_model
  implicit none
  private

  integer, parameter :: vacant_id = -1
  integer, parameter :: unset = -huge(1)
  character(*), parameter :: nucl_mask_name_suffix = '_nucleation_mask'
  character(*), parameter :: module_name = 'mod_nucleation'

  type, public :: nucleation_data
    type(timer_set) :: timer, expiration_template
    type(regular_grid) :: grid
    type(nucleus), allocatable :: nuclei(:)
    integer, allocatable :: eta_labels(:)
    type(order_parameter_t) :: eta_mask
    integer :: n_nuclei_max = unset
    integer :: queue_start_idx = 1
    integer :: queue_end_idx = 1
    real(wp) :: r_nodes, last_t
    real(wp), allocatable :: manual_prefactor, manual_activation_energy
    logical :: eta_mask_current = .false.
  contains
    procedure, public :: init
    procedure, public :: queue_is_empty
    procedure, public :: enqueue_nucleus
    procedure, public :: dequeue_nucleus
    procedure, public :: remove_expired_nuclei
    procedure, public :: calc_eta_mask
    procedure, public :: get_n_nuclei
    procedure, public :: get_n_vacant
    procedure, public :: get_nucl_candidates
  end type nucleation_data

contains

  subroutine init(this, grid, eta)
  class(nucleation_data), intent(inout) :: this
    type(order_parameter_t), intent(in) :: eta
    type(regular_grid), intent(in) :: grid
    character(*), parameter :: proc_name = 'init'

    if(this%n_nuclei_max == unset) call error_msg(proc_name, 'n_nuclei_max not set', module_name_opt = module_name)
    allocate(this%nuclei(this%n_nuclei_max))
    allocate(this%eta_labels(grid%n_nodes_total))

    this%eta_mask = eta
    this%eta_mask%real%time = 0.0_wp
    call this%eta_mask%fft%invalidate()
    this%nuclei%id = vacant_id
    this%grid = grid
    this%eta_labels = 0
  end subroutine init

  subroutine enqueue_nucleus(this, new_nucleus)
  class(nucleation_data), intent(inout) :: this
    type(nucleus), intent(in) :: new_nucleus
    character(*), parameter :: proc_name = 'enqueue_nucleus'

    if(this%queue_end_idx == mod_idx(this%queue_start_idx - 1, this%n_nuclei_max)) then
      call error_msg(proc_name, 'Nuclei queue overflow')
    end if

    associate(vacant_nucleus => this%nuclei(this%queue_end_idx))
      if(vacant_nucleus%id /= vacant_id) then
        call error_msg(proc_name, 'Nucleus at queue index ' // convert_to_char(this%queue_end_idx) // ' is occupied')
      end if
      vacant_nucleus = new_nucleus
    end associate
    call new_nucleus%apply_to_mask(this%eta_labels)
    this%queue_end_idx = mod_idx(this%queue_end_idx + 1, this%n_nuclei_max)
    this%eta_mask_current = .false.
  end subroutine enqueue_nucleus

  subroutine dequeue_nucleus(this)
  class(nucleation_data), intent(inout) :: this
    character(*), parameter :: proc_name = 'dequeue_nucleus'

    if(this%queue_is_empty()) then
      call error_msg(proc_name, 'Nuclei queue is empty')
    end if
    associate(obsolete_nucleus => this%nuclei(this%queue_start_idx))
      obsolete_nucleus%id = vacant_id
      call obsolete_nucleus%remove_from_mask(this%eta_labels)
    end associate
    this%queue_start_idx = mod_idx(this%queue_start_idx + 1, this%n_nuclei_max)
    this%eta_mask_current = .false.
  end subroutine dequeue_nucleus

  pure elemental logical function queue_is_empty(this)
  class(nucleation_data), intent(in) :: this

    queue_is_empty = this%queue_end_idx == this%queue_start_idx
  end function queue_is_empty

  subroutine remove_expired_nuclei(this, step)
  class(nucleation_data), intent(inout) :: this
    integer(int64), intent(in) :: step
    integer :: i

    if(this%queue_is_empty()) return
    do i = this%queue_start_idx, this%queue_end_idx - 1
      if(this%nuclei(i)%expiration%is_due(step)) call this%dequeue_nucleus()
    end do
  end subroutine remove_expired_nuclei

  subroutine calc_eta_mask(this, phase_idx)
  class(nucleation_data), intent(inout) :: this
    integer, intent(in) :: phase_idx

    if(this%eta_mask_current) return
    where(this%eta_labels == phase_idx)
      this%eta_mask%real%vals = 1.0_wp
    elsewhere(this%eta_labels == 0)
      this%eta_mask%real%vals = 0.0_wp
    elsewhere
      this%eta_mask%real%vals = -1.0_wp
    end where
    call this%eta_mask%fft%invalidate()
    call this%eta_mask%sharp_to_diffuse_interfaces()
    this%eta_mask_current = .true.
  end subroutine calc_eta_mask

  function get_nucl_candidates(this, phases, species_properties, nucl_phase_idx, temp_now, t_now, if_energy, delta_G_trans, &
      constants) result(nucl_candidates)

  class(nucleation_data), intent(in) :: this
    type(phase_set), intent(in) :: phases
    type(species_properties_t), intent(in) :: species_properties(:)
    integer, intent(in) :: nucl_phase_idx
    real(wp), intent(in) :: temp_now, t_now, if_energy, delta_G_trans(:)
    type(constants_t), intent(in) :: constants
    integer, allocatable, dimension(:) :: nucl_candidates
    real(wp) :: diff_limit, n_atoms_mix, atomic_volume_mix, lattice_param_mix, &
      nucl_rate, nucl_probability, prefactor, activation_energy
    real(wp), dimension(phases%n_phases) :: n_atoms, lattice_param, atomic_volume
    integer :: phase_idx, n_candidates, node
    logical :: override_prefactor, override_activation_energy
    type(arrhenius_model) :: atomic_rate
    character(*), parameter :: proc_name = 'get_nucl_candidates'

    allocate(nucl_candidates(this%get_n_vacant()))
    n_candidates = 0

    diff_limit = minval(species_properties%diffusion_coefficient%get_value(temp_now))
    do phase_idx = 1, phases%n_phases
      atomic_volume(phase_idx) = phases%phi(phase_idx)%properties%get_atomic_volume(constants%N_A)
      lattice_param(phase_idx) = phases%phi(phase_idx)%properties%get_lattice_parameter(constants%N_A)
    end do

    where(atomic_volume == huge(0.0_wp))
      n_atoms = 0.0_wp
      atomic_volume = 1E3*maxval(atomic_volume, atomic_volume < huge(0.0_wp))
    elsewhere
      ! This assumes that non-3d simulations have a "depth" corresponding to the minimum grid distance defined
      n_atoms = product(this%grid%axes%spacing)*minval(this%grid%axes%spacing)**(3 - this%grid%n_dims)/atomic_volume
      ! Alternatively, one could assume a single atomic layer, but nucleation probability is drastically reduced
      !n_atoms = product(this%grid%axes%spacing)/atomic_volume**(this%grid%n_dims/3.0_wp)
    end where

    where(lattice_param == huge(0.0_wp))
      lattice_param = 10*maxval(lattice_param, lattice_param < huge(0.0_wp))
    end where

    override_prefactor = allocated(this%manual_prefactor)
    override_activation_energy = allocated(this%manual_activation_energy)
    if(override_prefactor) prefactor = this%manual_prefactor

    do node = 1, this%grid%n_nodes_total
      if(delta_G_trans(node) >= 0) cycle
      n_atoms_mix = phases%phi_mix(n_atoms, node)
      atomic_volume_mix = phases%phi_mix(atomic_volume, node)
      lattice_param_mix = phases%phi_mix(lattice_param, node)
      if(.not. override_prefactor) then
        prefactor = zeldovich_factor(this%grid%n_dims, temp_now, if_energy, atomic_volume_mix, lattice_param_mix, &
          delta_G_trans(node), constants%k_B) &
          *frequency_factor(this%grid%n_dims, diff_limit, if_energy, lattice_param_mix, delta_G_trans(node))
      end if
      activation_energy = critical_energy(this%grid%n_dims, if_energy, delta_G_trans(node), lattice_param_mix)
      if(override_activation_energy) then
        activation_energy = activation_energy*this%manual_activation_energy
      end if
      atomic_rate = arrhenius_model(prefactor, activation_energy, constants%k_B)
      nucl_rate = (1 - phases%get_phi(nucl_phase_idx, node))*n_atoms_mix*atomic_rate%get_value(temp_now)
      nucl_probability = 1 - exp(-nucl_rate*(t_now - this%last_t))
      if(random_event(nucl_probability)) then
        print *, this_image(), node, 'YES'
        n_candidates = n_candidates + 1
        if(n_candidates > size(nucl_candidates)) then
          call error_msg(proc_name, 'Number of nucleation candidates larger than maximum number of nuclei')
        end if
        nucl_candidates(n_candidates) = node
      end if
    end do
    call resize_array(nucl_candidates, n_candidates)
  end function get_nucl_candidates

  pure integer function get_n_nuclei(this) result(n_nuclei)
  class(nucleation_data), intent(in) :: this

    if(this%queue_is_empty()) then
      n_nuclei = 0
    else
      n_nuclei = this%queue_end_idx - this%queue_start_idx
    end if
  end function get_n_nuclei

  pure integer function get_n_vacant(this) result(n_vacant)
  class(nucleation_data), intent(in) :: this

    n_vacant = this%n_nuclei_max - this%get_n_nuclei()
  end function get_n_vacant

  pure elemental real(wp) function critical_energy(n_dims, if_energy, delta_G_trans, lattice_param)
    integer, intent(in) :: n_dims
    real(wp), intent(in) :: if_energy, delta_G_trans, lattice_param
    character(*), parameter :: proc_name = 'critical_energy'

    if(delta_G_trans >= 0) then
      critical_energy = huge(0.0_wp)
      return
    end if
    select case(n_dims)
    case(1)
      critical_energy = 2*if_energy*lattice_param**2
    case(2)
      critical_energy = -pi*if_energy**2/delta_G_trans*lattice_param
    case(3)
      critical_energy = 16/3.0_wp*pi*if_energy**3/delta_G_trans**2
    case default
      critical_energy = 0.0_wp
      call error_msg(proc_name, 'Not defined for given amount of spatial dimensions: ' // convert_to_char(n_dims))
    end select
  end function critical_energy

  pure elemental real(wp) function critical_radius(n_dims, if_energy, delta_G_trans)
    integer, intent(in) :: n_dims
    real(wp), intent(in) :: if_energy, delta_G_trans
    character(*), parameter :: proc_name = 'critical_radius'

    if(delta_G_trans >= 0) then
      critical_radius = huge(0.0_wp)
      return
    end if
    select case(n_dims)
    case(1)
      critical_radius = 0.0_wp
    case(2)
      critical_radius = -if_energy/delta_G_trans
    case(3)
      critical_radius = -2*if_energy/delta_G_trans
    case default
      critical_radius = 0.0_wp
      call error_msg(proc_name, 'Not defined for given amount of spatial dimensions: ' // convert_to_char(n_dims))
    end select
  end function critical_radius

  pure elemental real(wp) function critical_n_atoms(n_dims, if_energy, delta_G_trans, atomic_volume, lattice_param)
    integer, intent(in) :: n_dims
    real(wp), intent(in) :: if_energy, delta_G_trans, atomic_volume, lattice_param
    character(*), parameter :: proc_name = 'critical_n_atoms'

    if(delta_G_trans >= 0) then
      critical_n_atoms = huge(0.0_wp)
      return
    end if
    select case(n_dims)
    case(1)
      critical_n_atoms = 2*critical_radius(n_dims, if_energy, delta_G_trans)*pi/atomic_volume*lattice_param**2
    case(2)
      critical_n_atoms = critical_radius(n_dims, if_energy, delta_G_trans)**2*pi/atomic_volume*lattice_param
    case(3)
      critical_n_atoms = 3/4.0_wp*critical_radius(n_dims, if_energy, delta_G_trans)**3*pi/atomic_volume
    case default
      critical_n_atoms = 0.0_wp
      call error_msg(proc_name, 'Not defined for given amount of spatial dimensions: ' // convert_to_char(n_dims))
    end select
  end function critical_n_atoms

  ! Formulae for zeldovich factor and frequency factor quoted from (modified):
  ! https://doi.org/10.1016/j.actamat.2015.08.041

  pure elemental real(wp) function &
      zeldovich_factor(n_dims, temp_now, if_energy, atomic_volume, lattice_param, delta_G_trans, k_B)

    integer, intent(in) :: n_dims
    real(wp), intent(in) :: temp_now, if_energy, delta_G_trans, atomic_volume, lattice_param, k_B
    character(*), parameter :: proc_name = 'zeldovich_factor'

    if(delta_G_trans >= 0) then
      zeldovich_factor = 0.0_wp
      return
    end if
    zeldovich_factor = sqrt(critical_energy(n_dims, if_energy, delta_G_trans, lattice_param) &
      /(4*pi*k_B*temp_now*critical_n_atoms(n_dims, if_energy, delta_G_trans, atomic_volume, lattice_param)**2))
  end function zeldovich_factor

  pure elemental real(wp) function frequency_factor(n_dims, diff_limit, if_energy, lattice_param, delta_G_trans)
    integer, intent(in) :: n_dims
    real(wp), intent(in) :: if_energy, delta_G_trans, diff_limit, lattice_param
    character(*), parameter :: proc_name = 'frequency_factor'

    if(delta_G_trans >= 0) then
      frequency_factor = 0.0_wp
      return
    end if
    select case(n_dims)
    case(1)
      frequency_factor = 2*pi/lattice_param**2*diff_limit
    case(2)
      frequency_factor = 2*critical_radius(n_dims, if_energy, delta_G_trans)*pi/lattice_param**3*diff_limit
    case(3)
      frequency_factor = 4*critical_radius(n_dims, if_energy, delta_G_trans)**2*pi/lattice_param**4*diff_limit
    case default
      frequency_factor = 0.0_wp
      call error_msg(proc_name, 'Not defined for given amount of spatial dimensions: ' // convert_to_char(n_dims))
    end select
  end function frequency_factor
end module mod_nucleation
