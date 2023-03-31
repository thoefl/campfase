project('jsonfortran', 'fortran', license: 'JSON-Fortran', default_options: ['fortran_std=f2008'], meson_version: '>=0.51')

subdir('src')

jsonfortran_lib = static_library('jsonfortran', source_files)
jsonfortran_dep = declare_dependency(link_with: jsonfortran_lib)
