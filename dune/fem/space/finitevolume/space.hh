#ifndef DUNE_FEM_SPACE_FINITEVOLUME_SPACE_HH
#define DUNE_FEM_SPACE_FINITEVOLUME_SPACE_HH

#include <dune/common/deprecated.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/function/localfunction/average.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/discontinuousgalerkin/generic.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>

#include "basisfunctionsets.hh"
#include "declaration.hh"

namespace Dune
{

  namespace Fem
  {

    // FiniteVolumeSpaceTraits
    // -----------------------

    template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
    struct FiniteVolumeSpaceTraits
    {
      typedef FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = codim;

      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;
      typedef FiniteVolumeBasisFunctionSets< EntityType, typename FunctionSpaceType::RangeType > BasisFunctionSetsType;
      typedef typename BasisFunctionSetsType::BasisFunctionSetType BasisFunctionSetType;

      typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;
      static const int localBlockSize = FunctionSpaceType::dimRange;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };



    // FiniteVolumeSpace
    // -----------------

    template< class FunctionSpace, class GridPart, int codim = 0, template< class > class Storage = SimpleStorage >
    class FiniteVolumeSpace
    : public GenericDiscontinuousGalerkinSpace< FiniteVolumeSpaceTraits< FunctionSpace, GridPart, codim, Storage > >
    {
      typedef FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > ThisType;
      typedef GenericDiscontinuousGalerkinSpace< FiniteVolumeSpaceTraits< FunctionSpace, GridPart, codim, Storage > > BaseType;

    public:
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::GridPartType */
      typedef typename BaseType::GridPartType GridPartType;
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::EntityType */
      typedef typename BaseType::EntityType EntityType;

      /** \brief basis function sets type */
      typedef typename BaseType::BasisFunctionSetsType BasisFunctionSetsType;

      explicit FiniteVolumeSpace ( GridPartType &gridPart,
                                   const InterfaceType commInterface = InteriorBorder_All_Interface,
                                   const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, BasisFunctionSetsType(), commInterface, commDirection )
      {}

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      static DFSpaceIdentifier type () { return FiniteVolumeSpace_id; }

      /** \copydoc Dune::Fem::GenericDiscontinuousGalerkinSpace::interpolate */
      template< class LocalFunction, class LocalDofVector >
      void interpolate ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        typename BaseType::RangeType value;
        LocalAverage< LocalFunction, GridPartType >::apply( localFunction, value );
        for( int i = 0; i < BaseType::FunctionSpaceType::dimRange; ++i )
          localDofVector[ i ] = value[ i ];
      }
    };



    // DefaultLocalRestrictProlong for FiniteVolumeSpace
    // -------------------------------------------------

    template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
    struct DefaultLocalRestrictProlong< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
    : public ConstantLocalRestrictProlong< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
    {
      DefaultLocalRestrictProlong ( const FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > & )
      {}
    };



    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct hasFixedPolynomialOrder< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct hasStaticPolynomialOrder< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
        static const int order = 0;
      };

      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct isContinuous< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct isLocalized< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false; // there is no method 'shapeFunctionSet( const EntityType & )'
      };

      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct isParallel< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::isParallel< GridPart >::v;
      };

      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct isAdaptive< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct threadSafe< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int codim, template< class > class Storage >
      struct viewThreadSafe< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FINITEVOLUME_SPACE_HH
