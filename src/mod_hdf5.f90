! Copyright 2020-2023 Thomas Hoefler
! This file is part of CAMPFASE.
! CAMPFASE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! CAMPFASE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with CAMPFASE.
! If not, see <https://www.gnu.org/licenses/>.

module mod_hdf5
    use hdf5
    use mod_base_functions
    use mod_axis
    use mod_timer_set
    use mod_globals, only: verbosity_level
    implicit none
    private
    
    character(*), parameter :: module_name = 'mod_hdf5'
    
    type, public :: hdf_file
        integer :: n_dims, compression_level
        logical :: save_on_dt_min_factor, extended
        type(timer_set) :: output, eta_output
        character(:), allocatable :: filename
        integer(hid_t) :: realtype, file_handle
        integer(hsize_t), allocatable, dimension(:) :: maxdims, startdims, arraychunk, scalarchunk
        integer(hsize_t), allocatable, dimension(:) :: input_dims, input_maxdims, input_offset
    contains
        procedure, public :: touch => hdf_file_touch
        procedure, public :: open => hdf_file_open
        procedure, public :: close => hdf_file_close
    end type hdf_file
        
    type, public :: hdf_set
        type(hdf_file), pointer :: file => null()
        integer(hid_t) :: datatype
        integer(hsize_t), allocatable :: chunksize(:)
        character(len=:), allocatable :: var_name
        integer(hid_t) :: dset_handle
        integer(int32) :: n_saved = 0
    contains
        procedure, private :: init_array => hdf_set_init_array
        procedure, private :: init_scalar => hdf_set_init_scalar
        procedure, private :: append_realarray => hdf_set_append_realarray
        procedure, private :: append_real => hdf_set_append_real
        procedure, private :: append_integer => hdf_set_append_integer
        procedure, private :: append_intarray => hdf_set_append_intarray
        procedure, private :: append_long => hdf_set_append_long
        procedure, private :: append_logical => hdf_set_append_logical
        generic, public :: init => init_array, init_scalar
        generic, public :: append => append_realarray, append_real, append_integer, append_intarray, append_long, append_logical
    end type hdf_set
    
    contains
    
    subroutine hdf_file_touch(this)
        class(hdf_file), intent(inout) :: this
        integer(int32) :: stat
        
        ! Create a HDF5 file
        call h5fcreate_f(this%filename, H5F_ACC_EXCL_F, this%file_handle, stat)
        call error_handler(stat, 'ERROR: Could not create HDF5 file.')
        call h5fclose_f(this%file_handle, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 file.')
    end subroutine hdf_file_touch
    
    subroutine hdf_file_open(this)
        class(hdf_file), intent(inout) :: this
        integer(int32) :: stat
        
        call h5fopen_f(this%filename, H5F_ACC_RDWR_F, this%file_handle, stat)
        call error_handler(stat, 'ERROR: Could not open HDF5 file.')
    end subroutine hdf_file_open
    
    subroutine hdf_file_close(this)
        class(hdf_file), intent(inout) :: this
        integer(int32) :: stat
        
        call h5fclose_f(this%file_handle, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 file.')
    end subroutine hdf_file_close

    subroutine hdf_set_init_scalar(this, file, var, var_name, datatype_opt)
        class(hdf_set), intent(inout) :: this
        type(hdf_file), intent(in), target :: file
        class(*), intent(in) :: var
        character(*), intent(in) :: var_name
        integer(hid_t), intent(in), optional :: datatype_opt
        integer(hsize_t) :: chunksize(file%n_dims)
        integer(hid_t) :: datatype
        integer(hid_t) :: dataspace
        integer(hid_t) :: props
        integer(int32) :: stat
        
        if(verbosity_level >= 4) print '(a)', 'Initializing HDF5 set ' // var_name
        call h5screate_simple_f(file%n_dims, file%startdims, dataspace, stat, file%maxdims)
        call error_handler(stat, 'ERROR: Could not create HDF5 dataspace.')
        
        if(present(datatype_opt)) then
            datatype = datatype_opt
        else
            select type(var)
            type is(logical)
                datatype = H5T_STD_I8LE
            type is(integer(int8))
                datatype = H5T_STD_I8LE
            type is(integer(int16))
                datatype = H5T_STD_I16LE
            type is(integer(int32))
                datatype = H5T_STD_I32LE
            type is(integer(int64))
                datatype = H5T_STD_I64LE
            type is(real(real32))
                datatype = H5T_IEEE_F32LE
            type is(real(real64))
                datatype = H5T_IEEE_F64LE
            class default
                error stop 'ERROR: hdf5 createset: Undefined variable type!'
            end select
        end if
    
        chunksize = 1
        
        ! Set hdf5 property list to change creation properties
        call h5pcreate_f(H5P_DATASET_CREATE_F, props, stat)
        call error_handler(stat,&
            'ERROR: Could not create HDF5 property list.')
        
        ! Enable compression of the data
        call h5pset_deflate_f(props, file%compression_level, stat)
        call error_handler(stat, 'ERROR: Could not enable HDF5 compression.')
        
        ! Make the data chunked
        call h5pset_chunk_f(props, file%n_dims, chunksize, stat) 
        call error_handler(stat, 'ERROR: Could not make HDF5 file chunked')
        
        ! Create the dataset, using the property set above
        call h5dcreate_f(file%file_handle, trim(var_name), datatype, &
            dataspace, this%dset_handle, stat, props)
        call error_handler(stat, 'ERROR: Could not create HDF5 dataset.')
        
        this%file => file
        this%datatype = datatype
        this%var_name = trim(var_name)
        this%chunksize = chunksize
        
        call h5pclose_f(props, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 property list.')
        call h5dclose_f(this%dset_handle, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataset.')
        call h5sclose_f(dataspace, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataspace.')
    end subroutine hdf_set_init_scalar
        
    subroutine hdf_set_init_array(this, file, spatial_dims, var, var_name, datatype_opt)
        class(hdf_set), intent(inout) :: this
        type(hdf_file), intent(in), target :: file
        integer, intent(in) :: spatial_dims(:)
        class(*), intent(in), dimension(:) :: var
        character(*), intent(in) :: var_name
        integer :: n_spatial_dims
        integer(hid_t), intent(in), optional :: datatype_opt
        integer(hid_t) :: datatype
        integer(hid_t) :: dataspace
        integer(hid_t) :: props
        integer(int32) :: stat
        
        if(verbosity_level >= 4) print '(a)', 'Initializing HDF5 set ' // var_name
        call h5screate_simple_f(file%n_dims, file%startdims, dataspace, stat, file%maxdims)
        call error_handler(stat, 'ERROR: Could not create HDF5 dataspace.')
        
        if(present(datatype_opt)) then
            datatype = datatype_opt
        else
            select type(var)
            type is(logical)
                datatype = H5T_STD_I8LE
            type is(integer(int8))
                datatype = H5T_STD_I8LE
            type is(integer(int16))
                datatype = H5T_STD_I16LE
            type is(integer(int32))
                datatype = H5T_STD_I32LE
            type is(integer(int64))
                datatype = H5T_STD_I64LE
            type is(real(real32))
                datatype = H5T_IEEE_F32LE
            type is(real(real64))
                datatype = H5T_IEEE_F64LE
            class default
                error stop 'ERROR: hdf5 createset: Undefined variable type!'
            end select
        end if
        
        this%file => file
        this%datatype = datatype
        this%var_name = var_name
        n_spatial_dims = size(spatial_dims)
        allocate(this%chunksize(file%n_dims))
        this%chunksize = 1
        this%chunksize(1:n_spatial_dims) = spatial_dims
        
        ! Set hdf5 property list to change creation properties
        call h5pcreate_f(H5P_DATASET_CREATE_F, props, stat)
        call error_handler(stat,&
            'ERROR: Could not create HDF5 property list.')
        
        ! Enable compression of the data
        call h5pset_deflate_f(props, file%compression_level, stat)
        call error_handler(stat, 'ERROR: Could not enable HDF5 compression.')
        
        ! Make the data chunked
        call h5pset_chunk_f(props, file%n_dims, this%chunksize, stat) 
        call error_handler(stat, 'ERROR: Could not make HDF5 file chunked')
        
        ! Create the dataset, using the property set above
        call h5dcreate_f(file%file_handle, var_name, datatype, &
            dataspace, this%dset_handle, stat, props)
        call error_handler(stat, 'ERROR: Could not create HDF5 dataset.')
        
        call h5pclose_f(props, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 property list.')
        call h5dclose_f(this%dset_handle, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataset.')
        call h5sclose_f(dataspace, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataspace.')
    end subroutine hdf_set_init_array
       
    subroutine hdf_set_append_realarray(this, a)
        class(hdf_set), intent(inout) :: this
        real(wp), intent(in) :: a(:)
        integer(hsize_t), dimension(this%file%n_dims) :: offset, write_chunk, expand_chunk
        integer(int32) :: stat
        integer(hid_t) :: dataspace
        integer(hid_t) :: memspace
        
        if(verbosity_level >= 4) print '(a)', 'Saving variable ' // this%var_name
        
        call h5dopen_f(this%file%file_handle, this%var_name, this%dset_handle, stat)
        call error_handler(stat, 'ERROR: Could not open HDF5 dataset.')
        
        this%n_saved = this%n_saved + 1
        offset = 0
        offset(this%file%n_dims) = this%n_saved - 1
        write_chunk = this%chunksize
        expand_chunk = this%chunksize
        expand_chunk(this%file%n_dims) = this%n_saved
        call h5dset_extent_f(this%dset_handle, expand_chunk, stat)
        call error_handler(stat, 'ERROR: Could not get HDF5 file extent.')
        call h5dget_space_f(this%dset_handle, dataspace, stat)
        call error_handler(stat, 'ERROR: Could not get HDF5 dataspace handle.')
        call h5screate_simple_f &
            (this%file%n_dims, write_chunk, memspace, stat)
        call error_handler(stat, 'ERROR: Could not create HDF5 memory space.')
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F,&
            offset, write_chunk, stat)
        call error_handler(stat,&
            'ERROR: Could not create HDF5 hyperslab selection.')
        call h5dwrite_f(this%dset_handle, this%datatype, a,&
            write_chunk, stat, memspace, dataspace)
        call error_handler(stat, 'ERROR: Could not write to HDF5 file.')
        
        call h5dclose_f(this%dset_handle, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataset.')
        call h5sclose_f(dataspace, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataspace.')
    end subroutine hdf_set_append_realarray
    
    subroutine hdf_set_append_real(this, a)
        class(hdf_set), intent(inout) :: this
        real(wp), intent(in) :: a
        integer(hsize_t), dimension(this%file%n_dims) :: offset, expand_chunk
        integer(int32) :: stat
        integer(hid_t) :: dataspace
        integer(hid_t) :: memspace
        
        if(verbosity_level >= 4) print '(a)', 'Saving variable ' // this%var_name
        
        call h5dopen_f(this%file%file_handle, this%var_name, this%dset_handle, stat)
        call error_handler(stat, 'ERROR: Could not open HDF5 dataset.')
        
        this%n_saved = this%n_saved + 1
        offset = 0
        offset(this%file%n_dims) = this%n_saved - 1
        expand_chunk = this%chunksize
        expand_chunk(this%file%n_dims) = this%n_saved
        call h5dset_extent_f(this%dset_handle, expand_chunk, stat)
        call error_handler(stat, 'ERROR: Could not get HDF5 file extent.')
        call h5dget_space_f(this%dset_handle, dataspace, stat)
        call error_handler(stat, 'ERROR: Could not get HDF5 dataspace handle.')
        call h5screate_simple_f &
            (this%file%n_dims, this%chunksize, memspace, stat)
        call error_handler(stat, 'ERROR: Could not create HDF5 memory space.')
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F,&
            offset, this%chunksize, stat)
        call error_handler(stat,&
            'ERROR: Could not create HDF5 hyperslab selection.')
        call h5dwrite_f(this%dset_handle, this%datatype, a,&
            this%chunksize, stat, memspace, dataspace)
        call error_handler(stat, 'ERROR: Could not write to HDF5 file.')
        call h5dclose_f(this%dset_handle, stat)
        
        call error_handler(stat, 'ERROR: Could not close HDF5 dataset.')
        call h5sclose_f(dataspace, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataspace.')
    end subroutine hdf_set_append_real
    
    subroutine hdf_set_append_integer(this, a)
        class(hdf_set), intent(inout) :: this
        integer(int32), intent(in) :: a
        integer(hsize_t) :: offset(this%file%n_dims)
        integer(int32) :: stat
        integer(hid_t) :: dataspace
        integer(hid_t) :: memspace
        integer(hsize_t) :: expand_chunk(this%file%n_dims)
        
        if(verbosity_level >= 4) print '(a)', 'Saving variable ' // this%var_name
        
        call h5dopen_f(this%file%file_handle, this%var_name, this%dset_handle, stat)
        call error_handler(stat, 'ERROR: Could not open HDF5 dataset.')
        
        this%n_saved = this%n_saved + 1
        offset = 0
        offset(this%file%n_dims) = this%n_saved - 1
        expand_chunk = this%chunksize
        expand_chunk(this%file%n_dims) = this%n_saved
        call h5dset_extent_f(this%dset_handle, expand_chunk, stat)
        call error_handler(stat, 'ERROR: Could not get HDF5 file extent.')
        call h5dget_space_f(this%dset_handle, dataspace, stat)
        call error_handler(stat, 'ERROR: Could not get HDF5 dataspace handle.')
        call h5screate_simple_f &
            (this%file%n_dims, this%chunksize, memspace, stat)
        call error_handler(stat, 'ERROR: Could not create HDF5 memory space.')
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F,&
            offset, this%chunksize, stat)
        call error_handler(stat,&
            'ERROR: Could not create HDF5 hyperslab selection.')
        call h5dwrite_f(this%dset_handle, this%datatype, a,&
            this%chunksize, stat, memspace, dataspace)
        call error_handler(stat, 'ERROR: Could not write to HDF5 file.')
        call h5dclose_f(this%dset_handle, stat)
        
        call error_handler(stat, 'ERROR: Could not close HDF5 dataset.')
        call h5sclose_f(dataspace, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataspace.')
    end subroutine hdf_set_append_integer
    
    subroutine hdf_set_append_intarray(this, a)
        class(hdf_set), intent(inout) :: this
        integer(int32), intent(in) :: a(:)
        integer(hsize_t) :: offset(this%file%n_dims)
        integer(int32) :: stat
        integer(hid_t) :: dataspace
        integer(hid_t) :: memspace
        integer(hsize_t) :: expand_chunk(this%file%n_dims)
        
        if(verbosity_level >= 4) print '(a)', 'Saving variable ' // this%var_name
        
        call h5dopen_f(this%file%file_handle, this%var_name, this%dset_handle, stat)
        call error_handler(stat, 'ERROR: Could not open HDF5 dataset.')
        
        this%n_saved = this%n_saved + 1
        offset = 0
        offset(this%file%n_dims) = this%n_saved - 1
        expand_chunk = this%chunksize
        expand_chunk(this%file%n_dims) = this%n_saved
        call h5dset_extent_f(this%dset_handle, expand_chunk, stat)
        call error_handler(stat, 'ERROR: Could not get HDF5 file extent.')
        call h5dget_space_f(this%dset_handle, dataspace, stat)
        call error_handler(stat, 'ERROR: Could not get HDF5 dataspace handle.')
        call h5screate_simple_f &
            (this%file%n_dims, this%chunksize, memspace, stat)
        call error_handler(stat, 'ERROR: Could not create HDF5 memory space.')
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F,&
            offset, this%chunksize, stat)
        call error_handler(stat,&
            'ERROR: Could not create HDF5 hyperslab selection.')
        call h5dwrite_f(this%dset_handle, this%datatype, a,&
            this%chunksize, stat, memspace, dataspace)
        call error_handler(stat, 'ERROR: Could not write to HDF5 file.')
        
        call h5dclose_f(this%dset_handle, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataset.')
        call h5sclose_f(dataspace, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataspace.')
    end subroutine hdf_set_append_intarray
    
    subroutine hdf_set_append_long(this, a)
        class(hdf_set), intent(inout) :: this
        integer(int64), intent(in) :: a
        integer(hsize_t) :: offset(this%file%n_dims)
        integer(int32) :: stat
        integer(hid_t) :: dataspace
        integer(hid_t) :: memspace
        integer(hsize_t) :: expand_chunk(this%file%n_dims)
        
        if(verbosity_level >= 4) print '(a)', 'Saving variable ' // this%var_name
        
        call h5dopen_f(this%file%file_handle, this%var_name, this%dset_handle, stat)
        call error_handler(stat, 'ERROR: Could not open HDF5 dataset.')
        
        this%n_saved = this%n_saved + 1
        offset = 0
        offset(this%file%n_dims) = this%n_saved - 1
        expand_chunk = this%chunksize
        expand_chunk(this%file%n_dims) = this%n_saved
        call h5dset_extent_f(this%dset_handle, expand_chunk, stat)
        call error_handler(stat, 'ERROR: Could not get HDF5 file extent.')
        call h5dget_space_f(this%dset_handle, dataspace, stat)
        call error_handler(stat, 'ERROR: Could not get HDF5 dataspace handle.')
        call h5screate_simple_f &
            (this%file%n_dims, this%chunksize, memspace, stat)
        call error_handler(stat, 'ERROR: Could not create HDF5 memory space.')
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F,&
            offset, this%chunksize, stat)
        call error_handler(stat,&
            'ERROR: Could not create HDF5 hyperslab selection.')
        call h5dwrite_f(this%dset_handle, this%datatype, a,&
            this%chunksize, stat, memspace, dataspace)
        call error_handler(stat, 'ERROR: Could not write to HDF5 file.')
        
        call h5dclose_f(this%dset_handle, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataset.')
        call h5sclose_f(dataspace, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataspace.')
    end subroutine hdf_set_append_long
    
    subroutine hdf_set_append_logical(this, a)
        class(hdf_set), intent(inout) :: this
        logical, intent(in) :: a
        integer(hsize_t) :: offset(this%file%n_dims)
        integer(int32) :: stat, a_int
        integer(hid_t) :: dataspace
        integer(hid_t) :: memspace
        integer(hsize_t) :: expand_chunk(this%file%n_dims)
        
        if(verbosity_level >= 4) print '(a)', 'Saving variable ' // this%var_name
        
        call h5dopen_f(this%file%file_handle, this%var_name, this%dset_handle, stat)
        call error_handler(stat, 'ERROR: Could not open HDF5 dataset.')
        
        if(a) then
            a_int = 1
        else
            a_int = 0
        end if
        this%n_saved = this%n_saved + 1
        offset = 0
        offset(this%file%n_dims) = this%n_saved - 1
        expand_chunk = this%chunksize
        expand_chunk(this%file%n_dims) = this%n_saved
        call h5dset_extent_f(this%dset_handle, expand_chunk, stat)
        call error_handler(stat, 'ERROR: Could not get HDF5 file extent.')
        call h5dget_space_f(this%dset_handle, dataspace, stat)
        call error_handler(stat, 'ERROR: Could not get HDF5 dataspace handle.')
        call h5screate_simple_f &
            (this%file%n_dims, this%chunksize, memspace, stat)
        call error_handler(stat, 'ERROR: Could not create HDF5 memory space.')
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F,&
            offset, this%chunksize, stat)
        call error_handler(stat,&
            'ERROR: Could not create HDF5 hyperslab selection.')
        call h5dwrite_f(this%dset_handle, this%datatype, a_int,&
            this%chunksize, stat, memspace, dataspace)
        call error_handler(stat, 'ERROR: Could not write to HDF5 file.')
        
        call h5dclose_f(this%dset_handle, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataset.')
        call h5sclose_f(dataspace, stat)
        call error_handler(stat, 'ERROR: Could not close HDF5 dataspace.')
    end subroutine hdf_set_append_logical
    
    subroutine error_handler(stat, errmsg)
        integer(int32), intent(in) :: stat
        character(*), intent(in) :: errmsg
        
        if(stat /= 0) error stop errmsg
    end subroutine error_handler
end module mod_hdf5