dune_install(basisfunctionset.hh codegen.hh default.hh evaluatecaller.hh evaluatecallerdeclaration.hh evaluatecallerdefaultimpl.hh
             functor.hh piolatransformation.hh proxy.hh simple.hh transformation.hh
             transformed.hh tuple.hh vectorial.hh)

exclude_from_headercheck( proxy.hh )

dune_add_subdirs(hpdg test)
