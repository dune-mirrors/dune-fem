#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SPACE_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SPACE_HH

// C++ includes
#include <algorithm>

// dune-common includes
#include <dune/common/nullptr.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
#include "dofmapper.hh"
#include "shapefunctionset.hh"

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

    static const int dimDomain = FunctionSpaceType::dimDomain;
    static const int dimRange = FunctionSpaceType::dimRange;
    static const int dimLocal = GridPartType::dimension;

    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;

    typedef typename GridPart::template Codim< 0 >::EntityType EntityType;

    typedef Dune::Fem::FunctionSpace< DomainFieldType, RangeFieldType, dimLocal, dimLocal > ShapeFunctionSpaceType;
    typedef ShapeFunctionSet< ShapeFunctionSpaceType, maxOrder > ShapeFunctionSetType;
    typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

    typedef typename MultiIndexSet< dimLocal, maxOrder >::MultiIndexType MultiIndexType;

    static const int localBlockSize = 1;
    typedef DofMapper< GridPartType, maxOrder > BlockMapperType;
    typedef BlockMapperType MapperType;

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
    typedef typename BaseType::GridType GridType;
    typedef typename BaseType::IndexSetType IndexSetType;
    typedef typename BaseType::IteratorType IteratorType;
    typedef typename BaseType::EntityType EntityType;
    typedef typename BaseType::IntersectionType IntersectionType;

    typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

    typedef typename BaseType::MapperType MapperType;
    typedef typename BaseType::BlockMapperType BlockMapperType;

    typedef typename Traits::MultiIndexType MultiIndexType;

  private:
    typedef typename Traits::ShapeFunctionSet ShapeFunctionSetType;

  public:
    using BaseType::gridPart;

    DiscreteFunctionSpace ( GridPartType &gridPart,
                            const MultiIndexType &multiIndex,
                            const Dune::InterfaceType commInterface,
                            const Dune::CommunicationDirection commDirection )
    : BaseType( gridPart, commInterface, commDirection ),
      mapper_( gridPart, multiIndex )
    { }

    Dune::Fem::DFSpaceIdentifier type () const
    {
      return Dune::Fem::GenericSpace_id;
    }

    BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
    {
      const MultiIndexType &multiIndex = mapper().multiIndex( entity );
      return BasisFunctionSetType( entity, ShapeFunctionSetType( multiIndex ) );
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
      return maxOrder;
    }

    int order ( const EntityType &entity ) const
    {
      const MultiIndexType &multiIndex = mapper().multiIndex( entity );
      return std::max_element( multiIndex.begin(), multiIndex.end() );
    }

    MapperType &mapper () const
    {
      return mapper_;
    }

    BlockMapperType &blockMapper () const
    {
      return mapper();
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

    MapperType mapper_;
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
      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::MultiIndexType MultiIndexType;

      AnisotropicDGSpace ( GridPartType &gridPart,
                           const MultiIndexType &multiIndex = MultiIndexType( maxOrder ),
                           const Dune::InterfaceType commInterface = defaultInterface,
                           const Dune::CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, multiIndex, commInterface, commDirection )
      { }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SPACE_HH
