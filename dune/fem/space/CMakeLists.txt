dune_install(
    brezzidouglasmarini.hh
    combinedspace.hh
    discontinuousgalerkin.hh
    finitevolume.hh
    fourier.hh
    lagrange.hh
    padaptivespace.hh
    rannacherturek.hh
    raviartthomas.hh
    p1bubble.hh)

dune_add_subdirs(
    basisfunctionset
    combinedspace
    common
    discontinuousgalerkin
    finitevolume
    fourier
    hpdg
    lagrange
    localfiniteelement
    mapper
    padaptivespace
    rannacherturek
    raviartthomas
    shapefunctionset
    test)

target_sources(dunefem PRIVATE
    shapefunctionset/legendrepolynomials.cc
    shapefunctionset/orthonormal/orthonormalbase_1d.cc
    shapefunctionset/orthonormal/orthonormalbase_2d.cc
    shapefunctionset/orthonormal/orthonormalbase_3d.cc
    lagrange/shapefunctionset.cc )
