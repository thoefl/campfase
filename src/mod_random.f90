! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_random
  use iso_fortran_env, only: int64
  use mod_base_functions, only: wp
  implicit none
  private

  public :: get_random_seed, restore_random_seed, scramble, random_int, random_reals, init_random_seed, random_event, random_real

  interface scramble
    module procedure :: scramble_int
    module procedure :: scramble_real
  end interface scramble

contains

  function get_random_seed() result(seed)
    integer, dimension(:), allocatable :: seed
    integer :: seed_size

    call random_seed(size = seed_size)
    allocate(seed(seed_size))
    call random_seed(get = seed)
  end function get_random_seed

  subroutine restore_random_seed(seed)
    integer, intent(in), dimension(:) :: seed
    integer :: seed_size

    call random_seed(size = seed_size)
    if(size(seed) /= seed_size) error stop &
      'ERROR: restore_random_seed: Size of seed provided different from current seed size!'
    call random_seed(put = seed)
  end subroutine restore_random_seed

  integer function random_int(min_val, max_val)
    integer, intent(in) :: min_val, max_val
    real(wp) :: rand

    call RANDOM_NUMBER(rand)
    random_int = floor(min_val + rand*(max_val - min_val + 1))
  end function random_int

  real(wp) function random_real(min_val, max_val)
    real(wp), intent(in) :: min_val, max_val

    call random_number(random_real)
    random_real = min_val + random_real*(max_val - min_val)
  end function random_real

  subroutine random_reals(rand, min_val, max_val)
    real(wp), intent(in) :: min_val, max_val
    real(wp), intent(inout) :: rand(:)

    call RANDOM_NUMBER(rand)
    rand = min_val + rand*(max_val - min_val)
  end subroutine random_reals

  subroutine scramble_int(a, n_rounds)
    integer, intent(inout) :: a(:)
    integer, intent(in) :: n_rounds
    integer :: round, idx_1, idx_2, store

    do round = 1, n_rounds
      do idx_1 = 1, size(a)
        idx_2 = random_int(1, size(a))
        store = a(idx_2)
        a(idx_2) = a(idx_1)
        a(idx_1) = store
      end do
    end do
  end subroutine scramble_int

  subroutine scramble_real(a, n_rounds)
    real(wp), intent(inout) :: a(:)
    integer, intent(in) :: n_rounds
    integer :: round, idx_1, idx_2
    real(wp) :: store

    do round = 1, n_rounds
      do idx_1 = 1, size(a)
        idx_2 = random_int(1, size(a))
        store = a(idx_2)
        a(idx_2) = a(idx_1)
        a(idx_1) = store
      end do
    end do
  end subroutine scramble_real

  ! Modified from: https://gcc.gnu.org/onlinedocs/gcc-4.9.1/gfortran/RANDOM_005fSEED.html
  subroutine init_random_seed(repeatable, image_distinct)
    logical, intent(in) :: repeatable, image_distinct
    integer(int64) :: scalar_seed
    integer, allocatable :: seed(:)
    integer :: i, seed_size, rand_discard
    integer(int64), parameter :: repeatable_seed = 104729
    integer, parameter :: n_discard = 1000

    call random_seed(size = seed_size)
    allocate(seed(seed_size))

    if(repeatable) then
      scalar_seed = repeatable_seed
    else
      call system_clock(count=scalar_seed)
    end if
    if(image_distinct) then
      do i = 1, (this_image() - 1)*seed_size
        call random_lcg(scalar_seed, rand_discard)
      end do
    else if(num_images() > 1 .and. .not. repeatable) then
      call co_broadcast(scalar_seed, 1)
    end if
    do i = 1, n_discard
      call random_lcg(scalar_seed, rand_discard)
    end do
    do i = 1, seed_size
      call random_lcg(scalar_seed, seed(i))
    end do
    call random_seed(put = seed)
  end subroutine init_random_seed

  ! Modified from: https://gcc.gnu.org/onlinedocs/gcc-4.9.1/gfortran/RANDOM_005fSEED.html
  subroutine random_lcg(seed, rand)
    integer(int64), intent(inout) :: seed
    integer, intent(out) :: rand
    if (seed == 0) then
      seed = 104729
    else
      seed = mod(seed, 4294967296_int64)
    end if
    seed = mod(seed * 279470273_int64, 4294967291_int64)
    rand = int(mod(seed, int(huge(0), int64)), kind(0))
  end subroutine random_lcg

  impure elemental logical function random_event(probability)
    real(wp), intent(in) :: probability
    real(wp) :: rand

    call random_number(rand)
    random_event = rand < probability
  end function
end module mod_random
