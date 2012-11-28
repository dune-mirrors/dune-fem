#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SPACE_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SPACE_HH

// C++ includes
#include <algorithm>

// dune-fem includes
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>

// local includes
#include "dofmapper.hh"
#include "shapefunctionsetstorage.hh"

/**
  @file
  @author Christoph Gersbacher
  @brief Please doc me
*/


namespace AnisotropicDG
{

  // Forward declaration
  // -------------------

  template< class FunctionSpace, class GridPart, int maxOrder, template< class > class Storage >
  class DiscreteFunctionSpace;



  // DiscreteFunctionSpaceTraits
  //----------------------------

  template< class FunctionSpace, class GridPart, int maxOrder, template< class > class Storage >
  struct DiscreteFunctionSpaceTraits
  {
    typedef AnisotropicDG::DiscreteFunctionSpace< FunctionSpace, GridPart, maxOrder, Storage > DiscreteFunctionSpaceType;

    typedef FunctionSpace FunctionSpaceType;
    typedef GridPart GridPartType;

    static const int dimRange = FunctionSpaceType::dimRange;
    static const int dimLocal = GridPartType::dimension;

    static const int codimension = 0;

    typedef typename MultiIndexSet< dimLocal, maxOrder >::MultiIndexType MultiIndexType;

    typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
    typedef typename Dune::Fem::ToLocalFunctionSpace< ScalarFunctionSpaceType, dimLocal >::Type ShapeFunctionSpaceType;

    typedef ShapeFunctionSetStorage< ShapeFunctionSpaceType, maxOrder, Storage > ScalarShapeFunctionSetStorageType;
    typedef typename ScalarShapeFunctionSetStorageType::ShapeFunctionSetType ScalarShapeFunctionSetType;

    typedef Dune::Fem::ShapeFunctionSetProxy< ScalarShapeFunctionSetType > ScalarShapeFunctionSetProxyType;
    typedef Dune::Fem::VectorialShapeFunctionSet< ScalarShapeFunctionSetProxyType, typename FunctionSpaceType::RangeType > ShapeFunctionSetType;

    typedef typename GridPart::template Codim< 0 >::EntityType EntityType;
    typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

    typedef DofMapper< GridPartType, maxOrder > BlockMapperType;
    static const int localBlockSize = dimRange;
    typedef Dune::Fem::NonBlockMapper< BlockMapperType, localBlockSize > MapperType;

    template< class DiscreteFunction, class Operation = Dune::Fem::DFCommunicationOperation::Copy >
    struct CommDataHandle
    {
      typedef Dune::Fem::DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      typedef Operation OperationType;
    };
  };



  // DiscreteFunctionSpace
  // ---------------------

  template< class FunctionSpace, class GridPart, int maxOrder, template< class > class Storage >
  struct DiscreteFunctionSpace
  : public Dune::Fem::DiscreteFunctionSpaceDefault< AnisotropicDG::DiscreteFunctionSpaceTraits< FunctionSpace, GridPart, maxOrder, Storage > >
  {
    typedef AnisotropicDG::DiscreteFunctionSpaceTraits< FunctionSpace, GridPart, maxOrder, Storage > Traits; 

  private:
    typedef AnisotropicDG::DiscreteFunctionSpace< FunctionSpace, GridPart, maxOrder, Storage > ThisType;
    typedef Dune::Fem::DiscreteFunctionSpaceDefault< Traits > BaseType;

  public:
    typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

    typedef typename BaseType::GridPartType GridPartType;
    typedef typename BaseType::EntityType EntityType;

    typedef typename Traits::ShapeFunctionSetType ShapeFunctionSetType;
    typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

    typedef typename BaseType::MapperType MapperType;
    typedef typename BaseType::BlockMapperType BlockMapperType;

    typedef typename Traits::MultiIndexType MultiIndexType;

  private:
    typedef typename Traits::ScalarShapeFunctionSetStorageType ScalarShapeFunctionSetStorageType;
    typedef typename Traits::ScalarShapeFunctionSetType ScalarShapeFunctionSetType;

  public:
    using BaseType::gridPart;

    DiscreteFunctionSpace ( GridPartType &gridPart,
                            const MultiIndexType &multiIndex,
                            const Dune::InterfaceType commInterface,
                            const Dune::CommunicationDirection commDirection )
    : BaseType( gridPart, commInterface, commDirection ),
      blockMapper_( gridPart, multiIndex ),
      mapper_( blockMapper_ ),
      multiIndex_( multiIndex )
    {
      scalarShapeFunctionSets_.insert( multiIndex_ );
    }

    Dune::Fem::DFSpaceIdentifier type () const
    {
      return Dune::Fem::GenericSpace_id;
    }

    BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
    {
      return BasisFunctionSetType( entity, shapeFunctionSet( entity ) );
    }

    ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity ) const
    {
      // const MultiIndexType &multiIndex = mapper().order( entity );
      return shapeFunctionSet( multiIndex_ );
    }
    
    ShapeFunctionSetType shapeFunctionSet ( const MultiIndexType &multiIndex ) const
    {
      return ShapeFunctionSetType( &scalarShapeFunctionSets_[ multiIndex ] );
    }
    
    bool contains ( const int codim ) const
    {
      return blockMapper().contains( codim );
    }

    bool continuous () const
    {
      return false;
    }

    int order () const
    {
      return maxOrder; // mapper().maxOrder();
    }

    int order ( const EntityType &entity ) const
    {
      // const MultiIndexType &multiIndex = mapper().order( entity );
      return *(std::max_element( multiIndex_.begin(), multiIndex_.end() ) );
    }

    MapperType &mapper () const
    {
      return mapper_;
    }

    BlockMapperType &blockMapper () const
    {
      return blockMapper_;
    }

    bool multipleGeometryTypes () const
    {
      return false;
    }

    bool multipleBasisFunctionSets () const
    {
      return true;
    }

  private:
    DiscreteFunctionSpace ( const ThisType & );
    ThisType &operator= ( const ThisType & );

    mutable BlockMapperType blockMapper_;
    mutable MapperType mapper_;
    MultiIndexType multiIndex_;
    ScalarShapeFunctionSetStorageType scalarShapeFunctionSets_;
  };

}  // namespace AnisotropicDG



namespace Dune
{

  namespace Fem
  {

    // AnisotropicDGSpace
    // ------------------

    template< class FunctionSpace, class GridPart, int maxOrder, template< class > class Storage = CachingStorage >
    class AnisotropicDGSpace
    : public AnisotropicDG::DiscreteFunctionSpace< FunctionSpace, GridPart, maxOrder, Storage >
    {
      typedef AnisotropicDGSpace< FunctionSpace, GridPart, maxOrder, Storage > ThisType;
      typedef AnisotropicDG::DiscreteFunctionSpace< FunctionSpace, GridPart, maxOrder, Storage > BaseType;

      static const Dune::InterfaceType defaultInterface = Dune::InteriorBorder_InteriorBorder_Interface;
      static const Dune::CommunicationDirection defaultDirection = Dune::ForwardCommunication;

    public:
      static const int polynomialOrder = maxOrder;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::MultiIndexType MultiIndexType;

      AnisotropicDGSpace ( GridPartType &gridPart,
                           const MultiIndexType &multiIndex = MultiIndexType( maxOrder ),
                           const Dune::InterfaceType commInterface = defaultInterface,
                           const Dune::CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, multiIndex, commInterface, commDirection )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SPACE_HH
