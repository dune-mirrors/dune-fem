#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_DOFMAPPER_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_DOFMAPPER_HH

// C++ includes
#include <cstddef>

// dune-fem includes
#include <dune/fem/space/mapper/dofmapper.hh>

// local includes
#include "multiindexset.hh"
#include "shapefunctionset.hh"

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
    //! \brief type of size integer
    typedef std::size_t SizeType;
  };



  // DofMapper
  // ---------

  template< class GridPart, int maxPolOrder >
  class DofMapper
  : public Dune::Fem::AdaptiveDofMapper< AnisotropicDG::DofMapperTraits< GridPart, maxPolOrder > >
  {
    typedef DofMapper< GridPart, maxPolOrder > ThisType;
    typedef Dune::Fem::AdaptiveDofMapper< AnisotropicDG::DofMapperTraits< GridPart, maxPolOrder > > BaseType;

  public:
    typedef typename BaseType::Traits Traits;
    typedef typename BaseType::ElementType ElementType;
    typedef typename BaseType::SizeType SizeType;

    typedef GridPart GridPartType;
    static const int dimension = GridPartType::dimension;

    typedef typename GridPartType::IndexSetType IndexSetType;
    typedef typename IndexSetType::IndexType GlobalKeyType;

    typedef typename MultiIndexSet< dimension, maxPolOrder >::MultiIndexType MultiIndexType;

    DofMapper ( const GridPartType &gridPart, const MultiIndexType multiIndex )
    : gridPart_( gridPart ),
      multiIndex_( multiIndex ),
      numDofs_( NumShapeFunctions< dimension, maxPolOrder >::count( multiIndex ) )
    {}


    // DofMapper interface methods
    // ---------------------------

    SizeType size () const
    {
      return indexSet().size( 0 )*numDofs_;
    }

    bool contains ( const int codim ) const
    {
      return ( codim == 0 );
    }

    bool fixedDataSize ( int codim ) const
    {
      return !contains( codim );
    }

    template< class Functor >
    void mapEach ( const ElementType &element, Functor f ) const
    {
      const int numDofs = ThisType::numDofs( element );
      for( int i = 0; i < numDofs; ++i )
        f( i, globalKey( element, i ) );
    }
    
    template< class Entity, class Functor >
    void mapEachEntityDof ( const Entity &entity, Functor f ) const
    {
      if( Entity::codimension != 0 )
        return;
      const int numDofs = ThisType::numDofs( entity );
      for( int i = 0; i < numDofs; ++i )
        f( i, globalKey( entity, i ) );
    }

    int maxNumDofs () const
    {
      return NumShapeFunctions< dimension, maxPolOrder >::max();
    }

    template< class Entity >
    int numDofs ( const Entity &entity ) const
    {
      assert( Entity::codimension == 0 );
      return numDofs_; 
    }

    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      return ( Entity::codimension == 0 ? numDofs( entity ) : 0 );
    }


    // Extended interface methods for AdaptviveDofMapper
    // -------------------------------------------------

    int numberOfHoles ( const int block ) const
    {
      return 0;
    }

    int oldIndex ( const int hole, const int block ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method oldIndex() not implemented yet" );
    }

    /** \brief return new index of hole for data block (with resprect to new offset) */
    int newIndex ( const int hole, const int block ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method newIndex() not implemented yet" );
    }

    /** \brief return true if compress will affect data */
    bool consecutive () const
    {
      return true;
    }

    /** \brief return old offsets for given block */
    int oldOffSet ( const int block ) const
    {
      return 0;
    }

    /** \brief return current offsets for given block */
    int offSet ( const int block ) const
    {
      return 0;
    }

    /** \brief return number of supported blocks */
    int numBlocks () const
    {
      return 1;
    }


    // Non-interface methods
    // ---------------------

    // return (multi index valued) polynomial order for element
    const MultiIndexType &order ( const ElementType &element ) const
    {
      return multiIndex_;
    }

    // return true, if only one shape function set is used 
    bool fixedOrder () const
    {
      return true;
    }
   
    // return maximum polyonmial order for given entity
    typename MultiIndexType::value_type maxOrder () const
    {
      return maxPolOrder;
    }

  protected:
    const GridPartType &gridPart () const
    {
      return gridPart_;
    }

    const IndexSetType &indexSet () const
    {
      return gridPart().indexSet();
    }

    template< class Entity >
    GlobalKeyType globalKey ( const Entity &entity, const int localDoF ) const
    {
      assert( Entity::codimension == 0 );
      return numDofs_*indexSet().index( entity ) + localDoF;
    }

  private:
    const GridPartType &gridPart_;
    MultiIndexType multiIndex_;
    SizeType numDofs_;
  };

} // namespace AnisotropicDG

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_DOFMAPPER_HH
