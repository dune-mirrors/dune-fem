dune_install(adaptiveleafgridpart.hh adaptiveleafindexset.hh codimindexset.hh
             filteredgridpart.hh filter.hh geogridpart.hh geometrygridpart.hh
             idgridpart.hh leafgridpart.hh levelgridpart.hh)

dune_add_subdirs(common filter filteredgridpart geogridpart idgridpart test geometrygridpart)
