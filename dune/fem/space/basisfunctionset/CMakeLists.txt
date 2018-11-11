dune_install(basisfunctionset.hh codegen.hh default.hh default_codegen.hh evaluatecaller.hh evaluatecaller_spec.hh
             functor.hh piolatransformation.hh proxy.hh simple.hh transformation.hh
             transformed.hh tuple.hh vectorial.hh)

exclude_from_headercheck( default_codegen.hh )

dune_add_subdirs(hpdg test)
