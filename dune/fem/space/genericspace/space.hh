#ifndef DUNE_FEM_SPACE_GENERICSPACE_SPACE_HH
#define DUNE_FEM_SPACE_GENERICSPACE_SPACE_HH

//- dune-fem includes
#include <dune/fem/space/basisfunctionset/basisfunctionset.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/dofmapper/indexsetdofmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/localfunctions.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class GridPart, class LocalFiniteElementFactory >
    class GenericDiscreteFunctionSpace;



    // GenericDiscreteFunctionSpaceTraits
    // ----------------------------------

    template< class GridPart, class LocalFiniteElementFactory >
    struct GenericDiscreteFunctionSpaceTraits
    {
      typedef GenericDiscreteFunctionSpace< GridPart, LocalFiniteElementFactory > DiscreteFunctionSpaceType;

      typedef GridPart GridPartType;

      static const int codimension = 0;
      static const int polynomialOrder = 111;

      static const int localBlockSize = 1;

      typedef IndexSetDofMapper< GridPartType > BlockMapperType;
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };

    private:
      friend class GenericDiscreteFunctionSpace< GridPart, LocalFiniteElementFactory >;

      typedef LocalFiniteElementFactory LocalFiniteElementFactoryType;

      typedef typename LocalFiniteElementFactory::LocalFiniteElementType LocalFiniteElementType;

      typedef typename LocalFiniteElementType::Traits::LocalBasis LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolation LocalInterpolationType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficients LocalCoefficientsType;

      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      typedef LocalFunctionsShapeFunctionSet< LocalBasisType > ShapeFunctionSetType;
      typedef DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetImplType;

    public:
      typedef BasisFunctionSet< EntityType, typename BasisFunctionSetImplType::RangeType, BasisFunctionSetImplType > BasisFunctionSetType;
    };



    // GenericDiscreteFunctionSpace
    // ----------------------------

    template< class GridPart, class LocalFiniteElementFactory >
    class GenericDiscreteFunctionSpace
    : public DiscreteFunctionSpaceDefault< GenericDiscreteFunctionSpaceTraits< GridPart, LocalFiniteElementFactory > >
    {
      typedef GenericDiscreteFunctionSpace< GridPart, LocalFiniteElementFactory > ThisType;
      typedef DiscreteFunctionSpaceDefault< GenericDiscreteFunctionSpaceTraits< GridPart, LocalFiniteElementFactory > > BaseType;

      typedef GenericDiscreteFunctionSpaceTraits< GridPart, LocalFiniteElementFactory > TraitsType;

      typedef typename TraitsType::BasisFunctionSetImplType BasisFunctionSetImplType;
      typedef typename TraitsType::ShapeFunctionSetType ShapeFunctionSetType;

    public:
      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename BaseType::BasisFunctionSpaceType BasisFunctionSetType;
      typedef typename BaseType::BlockMapperType BlockMapperType;
      typedef typename BaseType::MapperType MapperType;

      using BaseType::gridPart;

      explicit GenericDiscreteFunctionSpace ( GridPart &gridPart, InterfaceType commInterface, CommunicationDirection commDirection )
      : BaseType( gridPart, commInterface, commDirection ),
        mapper_( blockMapper_ )
      {}

      bool continuous () const { return false; }
      bool continuous ( const IntersectionType &intersection ) { return false; }

      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( BasisFunctionSetImplType( entity, shapeFunctionSet( entity ) ) );
      }

      BlockMapperType &blockMapper () const { return blockMapper_; }
      MapperType &mapper () const { return mapper_; }

    private:
      const ShapeFunctionSetType &shapeFunctionSet ( const EntityType &entity ) const
      {
      }

      // ...  shapeFunctionSets_;
      BlockMapperType blockMapper_;
      MapperType mapper_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_GENERICSPACE_SPACE_HH
