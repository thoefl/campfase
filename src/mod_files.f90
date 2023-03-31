! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_files
  use iso_fortran_env, only: int64
  implicit none
  private

  public :: copy_file
contains

  subroutine copy_file(in_filename, out_filename)
    character(*), intent(in) :: in_filename, out_filename
    integer(int64), parameter :: buffer_max_bytes = 4096
    character(:), allocatable :: buffer
    integer :: in_unit, out_unit
    integer(int64) :: byte_pos, total_bytes, n_full_buffers, remainder_bytes, i

    byte_pos = 1
    open(status = 'old', newunit = in_unit, file = in_filename, form = 'unformatted', access = 'stream', action = 'read')
    open(status = 'new', newunit = out_unit, file = out_filename, form = 'unformatted', access = 'stream', action = 'write')
    inquire(in_unit, size = total_bytes)
    n_full_buffers = total_bytes/buffer_max_bytes
    remainder_bytes = mod(total_bytes, buffer_max_bytes)
    if(n_full_buffers > 0) then
      allocate(character(buffer_max_bytes) :: buffer)
      do i = 1, n_full_buffers
        read(in_unit, pos = byte_pos) buffer
        write(out_unit) buffer
        byte_pos = byte_pos + buffer_max_bytes
      end do
      deallocate(buffer)
    end if
    if(remainder_bytes > 0) then
      allocate(character(remainder_bytes) :: buffer)
      read(in_unit, pos = byte_pos) buffer
      write(out_unit) buffer
    end if
    close(in_unit)
    close(out_unit)
  end subroutine copy_file
end module mod_files
