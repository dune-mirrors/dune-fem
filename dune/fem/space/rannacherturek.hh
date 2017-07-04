#ifndef DUNE_FEM_SPACE_RANNACHERTUREK_HH
#define DUNE_FEM_SPACE_RANNACHERTUREK_HH

#ifdef USE_OLD_RANNACHERTUREK_SPACE

#include <dune/fem/space/rannacherturek/space.hh>

#else

#include <dune/fem/space/localfiniteelement/space.hh>
#include <dune/fem/space/rannacherturek/localfemap.hh>

namespace Dune
{
  namespace Fem
  {

    template< class FunctionSpace, class GridPart, template< class > class Storage = CachingStorage >
    using RannacherTurekDiscreteFunctionSpace
    = LocalFiniteElementDiscreteFunctionSpace<
    RannacherTurekLocalFiniteElementMap< GridPart, FunctionSpace >, FunctionSpace, Storage >;

  }
}

#endif
#endif // #ifndef DUNE_FEM_SPACE_RANNACHERTUREK_HH
