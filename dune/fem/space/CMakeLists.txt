dune_install(combineddiscretefunctionspace.hh combinedspace.hh discontinuousgalerkin.hh
             finitevolume.hh fourier.hh lagrange.hh
             padaptivespace.hh rannacherturek.hh)

dune_add_subdirs(basisfunctionset combinedspace common
                 discontinuousgalerkin dofmapper finitevolume
                 fourier lagrange mapper
                 padaptivespace rannacherturek shapefunctionset
                 test)

dune_add_library(space OBJECT 
    shapefunctionset/legendrepolynomials.cc
    shapefunctionset/orthonormal/orthonormalbase_1d.cc
    shapefunctionset/orthonormal/orthonormalbase_2d.cc
    shapefunctionset/orthonormal/orthonormalbase_3d.cc
    lagrange/shapefunctionset.cc )
