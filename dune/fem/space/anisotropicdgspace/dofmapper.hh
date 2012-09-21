#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_DOFMAPPER_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_DOFMAPPER_HH

// C++ includes
#include <cstddef>

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-fem includes
#include <dune/fem/space/mapper/dofmapper.hh>

// local includes
#include "multiindexset.hh"
#include "utility.hh"

/**
  @file
  @author Christoph Gersbacher
  @brief Please doc me
*/


namespace AnisotropicDG
{

  // Internal forward declaration 
  // ----------------------------

  template< class GridPart, int maxPolOrder >
  class DofMapper;



  // DofMapperTraits
  // ---------------

  template< class GridPart, int maxPolOrder >
  struct DofMapperTraits
  {
    //! \brief type of DoF mapper
    typedef AnisotropicDG::DofMapper< GridPart, maxPolOrder > DofMapperType;
    //! \brief element type
    typedef typename GridPart::template Codim< 0 >::EntityType ElementType;
  };



  // DofMapper
  // ---------

  template< class GridPart, int maxPolOrder >
  class DofMapper
  : public Dune::Fem::DofMapper< AnisotropicDG::DofMapperTraits< GridPart, maxPolOrder > >
  {
    typedef DofMapper< GridPart, maxPolOrder > ThisType;
    typedef Dune::Fem::DofMapper< AnisotropicDG::DofMapperTraits< GridPart, maxPolOrder > > BaseType;

  public:
    typedef typename BaseType::Traits Traits;
    typedef typename BaseType::ElementType ElementType;

    typedef GridPart GridPartType;
    static const int dimension = GridPartType::dimension;

    typedef typename MultiIndexSet< dimension, maxPolOrder >::MultiIndexType MultiIndexType;

    explicit DofMapper ( const GridPartType &gridPart, const MultiIndexType multiIndex )
    : gridPart_( gridPart ),
      multiIndex_( multiIndex )
    {
      // compute total number of DoFs
      size_ = computeSize( gridPart_ );
    }


    // DofMapper interface methods
    // ---------------------------

    std::size_t size () const
    {
      return computeSize();
    }

    bool contains ( const int codim ) const
    {
      return ( codim == 0 );
    }

    bool fixedDataSize ( int codim ) const
    {
      return ( contains( codim ) ? false : true );
    }

    template< class Functor >
    void mapEach ( const ElementType &element, Functor f ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method mapEach() not implemented yet" );
    }
    
    template< class Entity, class Functor >
    void mapEachEntityDof ( const Entity &entity, Functor f ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method mapEachEntityDof() not implemented yet" );
    }

    int maxNumDofs () const
    {
      return size_;
    }

    template< class Entity >
    int numDofs ( const Entity &entity ) const
    {
      assert( Entity::codimension == 0 );
      return NumShapeFunctions< dimension, maxOrder >::count( multiIndex_ );
    }

    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      return ( Entity::codimension == 0 ? numDofs( entity ) : 0 );
    }


    // Non-interface methods
    // ---------------------

    // return (multi index valued) polynomial order for element
    MultiIndexType &order ( const ElementType &element ) const
    {
      return multiIndex_;
    }

    // return true, if only one shape function set is used 
    bool fixedOrder () const
    {
      return true;
    }
   
    // return maximum polyonmial order for given entity
    typename MultiIndexType ::value_type maxOrder () const
    {
      return maxPolOrder;
    }

  protected:
    const GridPartType &gridPart () const
    {
      return gridPart_;
    }

  private:
    // compute total number of DoFs
    std::size_t computeSize () const
    {
      std::size_t size = 0;
      typedef typename GridPartType::template Codim< 0 >::Iterator IteratorType;
      const IteratorType end = gridPart().template end< 0 >();
      for( IteratorType it = gridPart().template begin< 0 >(); it != end; ++it )
        size += numDofs( *it );
      return size;
    }

    const GridPartType &gridPart_;
    MultiIndexType multiIndex_;
    size_t size_;
  };

} // namespace AnisotropicDG

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_DOFMAPPER_HH
