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

/**
  @file
  @author Christoph Gersbacher
  @brief Please doc me
*/


namespace AnisotropicDG
{

  // Internal forward declaration 
  // ----------------------------

  template< class GridPart, int maxOrder >
  class DofMapper;



  // DofMapperTraits
  // ---------------

  template< class GridPart, int maxOrder >
  struct DofMapperTraits
  {
    //! \brief type of DoF mapper
    typedef AnisotropicDG::DofMapper< GridPart, maxOrder > DofMapperType;
    //! \brief element type
    typedef typename GridPart::template Codim< 0 >::EntityType ElementType;
  };



  // DofMapper
  // ---------

  template< class GridPart, int maxOrder >
  class DofMapper
  : public Dune::Fem::AdaptiveDofMapper< AnisotropicDG::DofMapperTraits< GridPart, maxOrder > >
  {
    typedef DofMapper< GridPart, maxOrder > ThisType;
    typedef Dune::Fem::AdaptiveDofMapper< AnisotropicDG::DofMapperTraits< GridPart, maxOrder > > BaseType;

  public:
    typedef typename BaseType::Traits Traits;
    typedef typename BaseType::ElementType ElementType;

    typedef GridPart GridPartType;
    static const int dimension = GridPartType::dimension;

    typedef typename MultiIndexSet< dimension, maxOrder >::MultiIndexType MultiIndexType;

    explicit DofMapper ( const GridPartType &gridPart )
    : gridPart_( gridPart ),
      size_( computeSize( gridPart_ ) )
    {}


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
      DUNE_THROW( Dune::NotImplemented, "Method maxNumDofs() not implemented yet" );
    }

    int numDofs ( const ElementType &element ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method numDofs() not implemented yet" );
    }

    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      return ( Entity::codimension == 0 ? numDofs() : 0 );
    }


    // AdaptvieDofMapper interface methods
    // -----------------------------------

    int numberOfHoles ( const int block ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method numberOfHoles() not implemented yet" );
    }

    int oldIndex ( const int hole, const int block ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method oldIndex() not implemented yet" );
    }

    int newIndex ( const int hole, const int block ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method newIndex() not implemented yet" );
    }

    bool consecutive () const
    {
      DUNE_THROW( Dune::NotImplemented, "Method consecutive() not implemented yet" );
    }

    int oldOffSet ( const int block ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method oldOffSet() not implemented yet" );
    }
   
    int offSet ( const int block ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method offSet() not implemented yet" );
    }
  
    int numBlocks () const
    {
      DUNE_THROW( Dune::NotImplemented, "Method numBlocks() not implemented yet" );
    }


    // Non-interface methods
    // ---------------------

    MultiIndexType &order ( const ElementType &element ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method order() not implemented yet" );
    }

    // return whether polynomial order is constant
    // bool fixedOrder ( const ElementType &element ) const;
   
    // return maximum polyonmial order for given entity
    // typename MultiIndexType maxOrder ( const ElementType &element ) const;

  protected:
    const GridPartType &gridPart () const
    {
      return gridPart_;
    }

  private:
    // compute total number of DoFs
    static std::size_t computeSize ( const GridPartType &gridPart )
    {
      std::size_t size = 0;
      typedef typename GridPartType::template Codim< 0 >::Iterator IteratorType;
      const IteratorType end = gridPart.template end< 0 >();
      for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
        size += numDofs( *it );
      return size;
    }

    const GridPartType &gridPart_;
    size_t size_;
  };

} // namespace AnisotropicDG

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_DOFMAPPER_HH
