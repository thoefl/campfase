project('campfase', 'fortran', version: '0.1', default_options: ['buildtype=custom', 'fortran_std=f2018'], license: 'GPL-3.0', meson_version: '>=0.63')
#project('campfase', 'fortran', version: '0.1', license: 'GPL-3.0', meson_version: '>=0.55')

jsonfortran_dep = dependency('jsonfortran', required: false)
if not jsonfortran_dep.found()
  jsonfortran_proj = subproject('jsonfortran')
  jsonfortran_dep = jsonfortran_proj.get_variable('jsonfortran_dep')
endif

ocasi_dep = dependency('OCASI', required: false)
if not ocasi_dep.found()
  ocasi_proj = subproject('opencalphad')
  ocasi_dep = ocasi_proj.get_variable('ocasi_dep')
endif

fftw_dep = dependency('fftw3')
fftw_omp_dep = dependency('fftw3_omp')
zlib_dep = dependency('zlib')
hdf_dep = dependency('hdf5', language : 'fortran')
openmp_dep = dependency('openmp')

add_project_arguments(['-fcoarray=lib'], language: 'fortran')

campfase_sources = files()
subdir('src')

buildtype = get_option('buildtype')
if buildtype == 'custom'
  add_project_arguments(['-Ofast', '-march=native', '-fopenmp'], language : 'fortran')
endif

campfase_dependencies = [zlib_dep, fftw_omp_dep, hdf_dep, ocasi_dep, jsonfortran_dep, openmp_dep]

fortran_comp = meson.get_compiler('fortran').get_id()
if fortran_comp == 'intel'
  message('Using intel compiler, adding flags...')
  add_project_arguments(['-parallel', '-xHost'], language : 'fortran')
elif fortran_comp == 'gcc'
  thread_dep = dependency('threads')
  coarray_dep = dependency('coarray')
  campfase_dependencies += [thread_dep, coarray_dep]
endif

executable('campfase', campfase_sources, dependencies : campfase_dependencies)

