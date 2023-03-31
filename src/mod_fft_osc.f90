! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_fft_osc
    
    use mod_base_functions
    use mod_regular_grid
    use mod_globals, only: verbosity_level
    
    implicit none
    private
    
    type, public :: fft_oscillation_data
        character(:), allocatable :: name
        integer :: node, n_osc, n_osc_max, error_idx, n_limit, error_nsave
        real(wp) :: errtol, prev_maxfreq
        real(wp), allocatable, dimension(:) :: error
    contains
        procedure :: init
        procedure :: update
        procedure :: get_error
        procedure :: is_limiting
    end type fft_oscillation_data
    
    type, public :: fft_grid_max_freq_oscillations
        type(regular_grid) :: real_grid, fft_grid
        integer :: n_freqs
        type(fft_oscillation_data), allocatable, dimension(:) :: max_freqs
        character(:), allocatable :: var_name
    contains
        procedure :: init => max_freq_init
        procedure :: update_all
    end type fft_grid_max_freq_oscillations
    
    contains
    
    pure subroutine init(this, name, node, errtol, n_osc_max, error_nsave)
        class(fft_oscillation_data), intent(inout) :: this
        character(*), intent(in) :: name
        integer, intent(in) :: node, error_nsave, n_osc_max
        real(wp), intent(in) :: errtol
        
        allocate(this%error(error_nsave))
        this%name = name
        this%node = node
        this%error = 0.0_wp
        this%error_idx = error_nsave
        this%error_nsave = error_nsave
        this%n_osc = 0
        this%errtol = errtol
        this%n_osc_max = n_osc_max
        this%prev_maxfreq = 0.0_wp
        this%n_limit = 0
    end subroutine init
    
    pure subroutine max_freq_init(this, var_name, real_grid, fft_grid, errtol, n_osc_max, error_nsave)
        class(fft_grid_max_freq_oscillations), intent(inout) :: this
        character(*), intent(in) :: var_name
        type(regular_grid), intent(in) :: real_grid, fft_grid
        real(wp), intent(in) :: errtol
        integer, intent(in) :: error_nsave, n_osc_max
        type(int_vector), allocatable, dimension(:) :: max_freq_coords, cartesian_coord_base
        type(string_vector), allocatable, dimension(:) :: max_freq_labels, cartesian_label_base
        type(string) :: name
        integer :: dim, freq, node
        character(*), parameter :: axes_label_separator = '-'
        
        allocate(cartesian_coord_base(real_grid%n_dims))
        allocate(cartesian_label_base(real_grid%n_dims))
        this%var_name = var_name
        do dim = 1, real_grid%n_dims
            allocate(cartesian_coord_base(dim)%v(2))
            allocate(cartesian_label_base(dim)%s(2))
            cartesian_coord_base(dim)%v(1) = 1
            cartesian_coord_base(dim)%v(2) = max_freq_node(real_grid%axes(dim)%n_nodes)
            cartesian_label_base(dim)%s(1)%c = ''
            cartesian_label_base(dim)%s(2)%c = real_grid%axes(dim)%label
        end do
        max_freq_coords = cartesian_product(cartesian_coord_base)
        max_freq_labels = cartesian_product(cartesian_label_base)
        ! Exclude node at coord = 1 for all axes (therefore also freq + 1 below)
        this%n_freqs = size(max_freq_coords) - 1
        allocate(this%max_freqs(this%n_freqs))
        do freq = 1, this%n_freqs
            node = fft_grid%dim_to_linear_idx(max_freq_coords(freq + 1)%v)
            name = concatenate_string_vector(max_freq_labels(freq + 1), axes_label_separator, discard_empty_opt = .true.)
            call this%max_freqs(freq)%init(name%c, node, errtol, n_osc_max, error_nsave)
        end do
    end subroutine max_freq_init
    
    pure subroutine update_all(this, a_fft)
        class(fft_grid_max_freq_oscillations), intent(inout) :: this
        complex(wp), intent(in) :: a_fft(:)
        integer :: i
        
        do i = 1, size(this%max_freqs)
            call this%max_freqs(i)%update(a_fft)
        end do
    end subroutine update_all
    
    pure subroutine update(this, a_fft)
        class(fft_oscillation_data), intent(inout) :: this
        complex(wp), intent(in) :: a_fft(:)
        real(wp) :: maxfreq
        
        !FIXME: Imaginary part relevant?
        maxfreq = a_fft(this%node)%re
        this%error_idx = modulo(this%error_idx, this%error_nsave) + 1
        this%error(this%error_idx) = maxfreq - this%prev_maxfreq
        if(this%error(this%error_idx) > this%errtol .and. this%n_osc <= 0) then
            this%n_osc = -this%n_osc + 1
        else if(this%error(this%error_idx) < -this%errtol .and. this%n_osc >= 0) then
            this%n_osc = -this%n_osc - 1
        else
            this%n_osc = 0
        end if
        this%prev_maxfreq = maxfreq
    end subroutine update
    
    pure function get_error(this) result(error)
        class(fft_oscillation_data), intent(in) :: this
        real(wp) :: error(this%error_nsave)
        
        error = cshift(this%error, this%error_idx - this%error_nsave)
    end function get_error
    
    pure integer elemental function max_freq_node(axis_n_nodes)
        integer, intent(in) :: axis_n_nodes
        
        max_freq_node = (axis_n_nodes//2) + 1
    end function max_freq_node
    
    impure elemental logical function is_limiting(this)
        class(fft_oscillation_data), intent(inout) :: this

        is_limiting = abs(this%n_osc) > this%n_osc_max
        if(is_limiting) then
            if(verbosity_level >= 3) call print_message('Oscillations detected in ' // trim(this%name) // ' FFT')
            this%n_limit = this%n_limit + 1
        end if
    end function is_limiting
end module mod_fft_osc
