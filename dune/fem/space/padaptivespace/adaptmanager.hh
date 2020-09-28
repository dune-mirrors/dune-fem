#ifndef DUNE_FEM_SPACE_PADAPTIVE_ADAPTMANAGER_HH
#define DUNE_FEM_SPACE_PADAPTIVE_ADAPTMANAGER_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>

#include "declaration.hh"
#include "restrictprolong.hh"


namespace Dune
{

  namespace Fem
  {

    // DefaultLocalRestrictProlong
    // ---------------------------

    template< class FS, class GP, int ord, class S >
    class DefaultLocalRestrictProlong< Fem::PAdaptiveLagrangeSpace< FS, GP, ord, S > >
    : public PLagrangeLocalRestrictProlong< typename GP::GridType, Fem::PAdaptiveLagrangeSpace< FS, GP, ord, S > >
    {
    public:
      DefaultLocalRestrictProlong ( const Fem::PAdaptiveLagrangeSpace< FS, GP, ord, S > &space )
      : PLagrangeLocalRestrictProlong< typename GP::GridType, Fem::PAdaptiveLagrangeSpace< FS, GP, ord, S > >( space )
      {}
    };


    template< class FunctionSpaceImp, class GridPartImp, int polOrd, class StorageImp >
    class DefaultLocalRestrictProlong< Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > >
    : public DiscontinuousGalerkinLocalRestrictProlong< Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp >, false > // invert mass matrix or not
    {
    public:
      typedef DiscontinuousGalerkinLocalRestrictProlong< Fem::PAdaptiveDGSpace<
                                                         FunctionSpaceImp,
                                                         GridPartImp,
                                                         polOrd, StorageImp >, false > BaseType ;
      DefaultLocalRestrictProlong ( const Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > & space )
        : BaseType( space )
      {}
    };

    template< class FunctionSpaceImp, class GridPartImp, class StorageImp >
    class DefaultLocalRestrictProlong< Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
    : public ConstantLocalRestrictProlong< Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
    {
    public:
      DefaultLocalRestrictProlong ( const Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > & )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_ADAPTMANAGER_HH
