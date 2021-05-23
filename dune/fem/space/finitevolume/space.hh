#ifndef DUNE_FEM_SPACE_FINITEVOLUME_SPACE_HH
#define DUNE_FEM_SPACE_FINITEVOLUME_SPACE_HH

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/common/hybrid.hh>
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
#include "interpolation.hh"

namespace Dune
{

  namespace Fem
  {

    // FiniteVolumeSpaceTraits
    // -----------------------

    template< class FunctionSpace, class GridPart, int codim, class Storage >
    struct FiniteVolumeSpaceTraits
    {
      typedef FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > DiscreteFunctionSpaceType;

      typedef GridPart GridPartType;
      typedef GridFunctionSpace< GridPartType, FunctionSpace > FunctionSpaceType;

      static const int codimension = codim;

      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;
      typedef FiniteVolumeBasisFunctionSets< EntityType, typename FunctionSpaceType::RangeType > BasisFunctionSetsType;
      typedef typename BasisFunctionSetsType::BasisFunctionSetType BasisFunctionSetType;

      typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;
      typedef Hybrid::IndexRange< int, FunctionSpaceType::dimRange > LocalBlockIndices;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };



    // FiniteVolumeSpace
    // -----------------

    template< class FunctionSpace, class GridPart, int codim = 0, class Storage = SimpleStorage >
    class FiniteVolumeSpace
    : public GenericDiscontinuousGalerkinSpace< FiniteVolumeSpaceTraits< FunctionSpace, GridPart, codim, Storage > >
    {
      typedef FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > ThisType;
      typedef GenericDiscontinuousGalerkinSpace< FiniteVolumeSpaceTraits< FunctionSpace, GridPart, codim, Storage > > BaseType;

    public:
      /** \brief maximum polynomial order of the space, here 0 since basis functions are constant */
      static const int polynomialOrder = 0;

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::GridPartType */
      typedef typename BaseType::GridPartType GridPartType;
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::EntityType */
      typedef typename BaseType::EntityType EntityType;

      /** \brief basis function sets type */
      typedef typename BaseType::BasisFunctionSetsType BasisFunctionSetsType;
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::BasisFunctionSetType */
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      /** \brief local interpolation type */
      typedef FiniteVolumeLocalInterpolation< GridPart, typename BasisFunctionSetType::RangeType > InterpolationType;

      explicit FiniteVolumeSpace ( GridPartType &gridPart,
                                   const InterfaceType commInterface = InteriorBorder_All_Interface,
                                   const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, BasisFunctionSetsType(), commInterface, commDirection )
      {}

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      static DFSpaceIdentifier type () { return FiniteVolumeSpace_id; }

      /** \brief return local interpolation */
      InterpolationType interpolation () const
      {
        return InterpolationType();
      }

      /** \brief return local interpolation */
      static InterpolationType interpolation ( const EntityType &entity )
      {
        return InterpolationType();
      }

      /** \brief extend size of space beyond what the index set is delivering */
      void extendSize( const size_t extension ) { this->blockMapper().extendSize( extension ); }
    };


   /** \brief Local Mass Matrix for FV space */
    template <class FunctionSpaceImp, class GridPartImp, int polOrd,
              class BaseFunctionStorageImp,
              class VolumeQuadratureImp>
    class LocalMassMatrix<
      FiniteVolumeSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >,
      VolumeQuadratureImp >
      : public LocalMassMatrixImplementationDgOrthoNormal<
          FiniteVolumeSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >, VolumeQuadratureImp >
    {
      typedef FiniteVolumeSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp > DiscreteFunctionSpaceImp;
      typedef LocalMassMatrixImplementationDgOrthoNormal< DiscreteFunctionSpaceImp, VolumeQuadratureImp > BaseType;
    public:
      using BaseType :: BaseType;
    };





    // DefaultLocalRestrictProlong for FiniteVolumeSpace
    // -------------------------------------------------

    template< class FunctionSpace, class GridPart, int codim, class Storage >
    class DefaultLocalRestrictProlong< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
    : public ConstantLocalRestrictProlong< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
    {
    public:
      DefaultLocalRestrictProlong ( const FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > & )
      {}
    };



    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct hasFixedPolynomialOrder< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct hasStaticPolynomialOrder< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
        static const int order = 0;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct isContinuous< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct isHierarchic< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct isLocalized< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false; // there is no method 'shapeFunctionSet( const EntityType & )'
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct isAdaptive< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct threadSafe< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct viewThreadSafe< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FINITEVOLUME_SPACE_HH
