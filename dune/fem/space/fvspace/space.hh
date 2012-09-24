#ifndef DUNE_FEM_SPACE_FVSPACE_SPACE_HH
#define DUNE_FEM_SPACE_FVSPACE_SPACE_HH

// dune-common includes
#include <dune/common/typetraits.hh>

// dune-fem includes
#include <dune/fem/space/common/discretefunctionspace.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

// local includes
#include "shapefunctionset.hh"

/*
  \file
  \brief Finite Volume space
  \author Christoph Gersbacher
*/


namespace
{
  template< template< class > class Storage >
  struct ShowWarning;

  template<>
  struct ShowWarning< Dune::Fem::CachingStorage >
  {
    static const bool v = true;
  };

  template<>
  struct ShowWarning< Dune::Fem::SimpleStorage >
  {
    static const bool v = false;
  };

} // namespace

namespace Dune
{

  namespace Fem
  {

    // Forward declaration
    // -------------------

    template< class FunctionSpace, class GridPart, int polOrder,
              template< class > class Storage >
    class FiniteVolumeSpace;


    // FiniteVolumeSpaceTraits
    // -----------------------

    template< class FunctionSpace, class GridPart, int polOrder,
              template< class > class Storage >
    struct FiniteVolumeSpaceTraits
    {
      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      typedef typename GridPartType::IndexSetType IndexSetType;

      static const int codimension = 0;
      typedef typename GridPartType::template Codim< codimension >::IteratorType IteratorType;
      typedef typename IteratorType::Entity EntityType;
      
      static const int dimRange = FunctionSpaceType::dimRange;
      static const int dimLocal = GridPartType::dimension; 

      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;

      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

      typedef FiniteVolumeSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef FiniteVolumeShapeFunctionSet< typename ToLocalFunctionSpace< FunctionSpace, dimLocal >::Type, polOrder > ShapeFunctionSetType;
      typedef DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

      static const int localBlockSize = dimRange;

      typedef CodimensionMapper< GridPartType, 0 > BlockMapperType;
      typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };



    // FiniteVolumeSpace
    // -----------------
 
    template< class FunctionSpace, class GridPart, int polOrder,
              template< class > class Storage = SimpleStorage >
    struct FiniteVolumeSpace
    : public DiscreteFunctionSpaceDefault< FiniteVolumeSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef FiniteVolumeSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > Traits;

    private:
      typedef FiniteVolumeSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

      static const InterfaceType defaultInterface = InteriorBorder_All_Interface;
      static const CommunicationDirection defaultDirection =  ForwardCommunication;

    public:
      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::GridType GridType;
      typedef typename BaseType::IndexSetType IndexSetType;
      typedef typename BaseType::IteratorType IteratorType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::MapperType MapperType;
      typedef typename BaseType::BlockMapperType BlockMapperType;

      explicit FiniteVolumeSpace( GridPartType &gridPart,
                                  const InterfaceType commInterface = defaultInterface,
                                  const CommunicationDirection commDirection = defaultDirection )
      {
        deprecationWarning( Dune::integral_constant< bool, ShowWarning< Storage >::v >() );
      }

      DFSpaceIdentifier type () const
      {
        return FiniteVolumeSpace_id;
      }

      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity );
      }

      bool contains ( const int codim ) const
      {
        return blockMapper().contains( codim );
      }

      bool continuous () const
      {
        return ( polOrder != 0 );
      }

      int order () const
      {
        return polOrder;
      }

      MapperType &mapper () const
      {
        return mapper_;
      }

      BlockMapperType &blockMapper () const
      {
        return blockMapper_;
      }

    private:
      void DUNE_DEPRECATED_MSG( "Caching disabled for FiniteVolumeSpace." )
      deprecationWarning ( Dune::integral_constant< bool, true > ) {}

      void
      deprecationWarning ( Dune::integral_constant< bool, false > ) {}        

      BlockMapperType& blockMapper_;
      mutable MapperType mapper_; 
    };



    // DefaultLocalRestrictProlong for FiniteVolumeSpace
    // -------------------------------------------------

    template< class FunctionSpace, class GridPart, template< class > class Storage >
    struct DefaultLocalRestrictProlong< FiniteVolumeSpace< FunctionSpace, GridPart, 0, Storage > >
    : public ConstantLocalRestrictProlong< FiniteVolumeSpace< FunctionSpace, GridPart, 0, Storage > >
    {
      DefaultLocalRestrictProlong ( const FiniteVolumeSpace< FunctionSpace, GridPart, 0, Storage > & )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FVSPACE_SPACE_HH
