dune_install(basicfilterwrapper.hh
             domainfilter.hh
             filter.hh
             inversefilter.hh
             radialfilter.hh
             simple.hh
             threadfilter.hh
             voidfilter.hh)

exclude_from_headercheck(
  threadfilter.hh
  )
