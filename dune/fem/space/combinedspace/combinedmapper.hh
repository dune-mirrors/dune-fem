#ifndef DUNE_COMBINEDMAPPER_HH
#define DUNE_COMBINEDMAPPER_HH

#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/common/dofmanager.hh>



namespace Dune
{
  namespace Fem
  {

    //! forward declaration
    template< class Grid,  class Mapper1, class Mapper2 >
    class CombinedMapper;

    //! Traits 
    template< class Grid, class Mapper1, class Mapper2 >
    struct CombinedMapperTraits
    {
      // we still need entities of codimension 0 here 
      typedef typename Mapper1 :: ElementType ElementType;

      static const int polynomialOrder1 = Mapper1 :: polynomialOrder;
      static const int polynomialOrder2 = Mapper2 :: polynomialOrder;
      static const int polynomialOrder = ( polynomialOrder1 > polynomialOrder2 ) ? polynomialOrder1 : polynomialOrder2;

      typedef CombinedMapper< Grid, Mapper1, Mapper2 >  DofMapperType;

      // Later on write an new dofIterator
      typedef DefaultDofMapIterator< ElementType, DofMapperType > DofMapIteratorType;
    };

    template< class Grid, class Mapper1, class Mapper2 >
    class CombinedMapper
    : public DofMapperDefault< CombinedMapperTraits< Grid, Mapper1, Mapper2 > >
    {
      typedef CombinedMapper< Grid, Mapper1, Mapper2 > ThisType;
      typedef DofMapperDefault< CombinedMapperTraits< Grid, Mapper1, Mapper2 > > BaseType;
      typedef Grid GridType;

      typedef Mapper1 MapperType1;
      typedef Mapper2 MapperType2;

      template< class Functor >
      struct FunctorWrapper
      {

        FunctorWrapper( Functor functor, int localOffset, int globalOffset ) :
          functor_( functor ) ,         
          localOffset_(  localOffset ) ,
          globalOffset_( globalOffset )
        {}

        void operator() ( int localDof, int globalDof ) 
        {
          functor_( localDof + localOffset_, globalDof + globalOffset_);
        }
        private:
          Functor functor_;
          int localOffset_;
          int globalOffset_;
      };

      public:
      typedef typename BaseType :: Traits Traits;
      typedef typename BaseType :: ElementType ElementType;

      typedef typename Traits :: DofMapIteratorType DofMapIteratorType;

      //! order of the Lagrange polynoms
      static const int polynomialOrder1 = Traits :: polynomialOrder1;
      static const int polynomialOrder2 = Traits :: polynomialOrder2;
      static const int polynomialOrder = Traits :: polynomialOrder;

      //! type of the DoF manager
      typedef DofManager< GridType > DofManagerType;

      public:
      //! constructor
      CombinedMapper( const GridType &grid,  const MapperType1 &mapper1, const MapperType2& mapper2 )
        : // BaseType( gridPart ), // hmmm
          dm_( DofManagerType :: instance( grid ) ), 
          mapper1_( mapper1 ),
          mapper2_( mapper2 ),
          globalOffset_( mapper1.size() ),
          oldGlobalOffset_( -1 )
      {
        dm_.addIndexSet( *this );
      }

      ~CombinedMapper() 
      {
        dm_.removeIndexSet( *this );
      }
      
      //! copy constructor
      CombinedMapper( const ThisType &other )
        : dm_( other.dm_ ),
          mapper1_( other.mapper1_ ),
          mapper2_( other.mapper2_ ),
          globalOffset_( other.globalOffset_ ),
          oldGlobalOffset_( other.oldGlobalOffset_ )
      {
        dm_.addIndexSet( *this );
      }

      /** \copydoc Dune::DofMapper::size() const */
      int size () const 
      {
        return mapper1_.size() + mapper2_.size();
      }

      /** \copydoc Dune::DofMapper::begin( const ElementType &element ) const */
      DofMapIteratorType begin ( const ElementType &element ) const   
      {
        return DofMapIteratorType( DofMapIteratorType::beginIterator, element, *this );
      }
      
      /** \copydoc Dune::DofMapper::end( const ElementType &element ) const */
      DofMapIteratorType end ( const ElementType &element ) const
      {
        return DofMapIteratorType( DofMapIteratorType::endIterator, element, *this );
      }

      /** \copydoc Dune::DofMapper::contains( const int codim ) const */
      bool contains ( const int codim ) const    
      {
        return ( mapper1_.contains( codim ) || mapper2_.contains( codim ) );
      }

      /** \copydoc Dune::DofMapper::mapTpGlobal( const ElementType &element, const int localDof) const */
      int mapToGlobal ( const ElementType &element, const int localDof ) const
      {
        assert( mapper1_.size() == globalOffset_ );
        const int localOffset = mapper1_.numDofs(element);

        int index;

        if( localDof - localOffset  < 0 )
          index = mapper1_.mapToGlobal( element, localDof );
        else
          index = mapper2_.mapToGlobal( element, localDof - localOffset ) + globalOffset_;
          
        assert( (0 <= index) && (index < size()) );
        return index;
      }
      
      /** \copydoc Dune::DofMapper::mapEntityDofToGlobal( const Entity &entity, const int localDof ) const */
      template< class Entity > 
      int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const
      {
        assert( mapper1_.size() == globalOffset_ );
        const int numEntityDofs =  mapper1_.numEntityDofs( entity);

        if( localDof - numEntityDofs < 0 )
          return mapper1_.mapEntityDofToGlobal( entity, localDof );
        else
          return mapper2_.mapEntityDofToGlobal( entity, localDof - numEntityDofs ) + globalOffset_;
      }
      
      /** \copydoc Dune::DofMapper::maxNumDofs() const */
      int maxNumDofs () const
      {
        return mapper1_.maxNumDofs() + mapper2_.maxNumDofs();
      }


      /** \copydoc Dune::DofMapper::mapEach( const ElementType &element, Functor f ) const */
      template< class Functor >
      void mapEach ( const ElementType &element, Functor f ) const
      {
        mapper1_.mapEach( element, FunctorWrapper< Functor > (f, 0, 0) );
        mapper2_.mapEach( element, FunctorWrapper< Functor > (f, mapper1_.numDofs( element ), globalOffset_) );
      }

      /** \copydoc Dune::DofMapper::numDofs( const ElementType &element ) const */
      int numDofs ( const ElementType &element ) const
      {
        int nDofs = mapper1_.numDofs( element ) + mapper2_.numDofs( element );
  //      assert( (nDofs >=0) && (nDofs < size()) );
        return nDofs;
      }

      /** \copydoc Dune::DofMapper::numEntityDofs( const Entity &entity ) const */
      template< class Entity >
      int numEntityDofs ( const Entity &entity ) const
      {
        return mapper1_.numEntityDofs( entity ) + mapper2_.numEntityDofs( entity );
      }

      /** \copydoc Dune::DofMapper::fixedSize( int codim ) const */
      bool fixedDataSize ( int codim ) const
      {
        return mapper1_.fixedDataSize( codim ) && mapper2_.fixedDataSize( codim );
      }

      /** \copydoc Dune::DofMapper::numerOfHoles( const int block ) const */
      int numberOfHoles(const int block) const 
      {
        const int numBlock1 = mapper1_.numBlocks(); 
        if ( block - numBlock1 < 0 ) 
          return mapper1_.numberOfHoles( block );
        else
          return mapper2_.numberOfHoles( block - numBlock1 );
      }
      
      /** \copydoc Dune::DofMapper::oldIndex( const int hole, const int block ) const */
      int oldIndex (const int hole, const int block) const 
      { 
        const int numBlock1 = mapper1_.numBlocks(); 
        if ( block - numBlock1 < 0 ) 
          return mapper1_.oldIndex( hole, block );
        else
          return mapper2_.oldIndex( hole, block - numBlock1 );
      }
        
      /** \copydoc Dune::DofMapper::newIndex( const int hole, const int block ) const */
      int newIndex (const int hole, const int block) const 
      { 
        const int numBlock1 = mapper1_.numBlocks(); 
        if ( block - numBlock1 < 0 ) 
          return mapper1_.newIndex( hole, block );
        else
          return mapper2_.newIndex( hole, block - numBlock1 );
      }

      /** \copydoc Dune::DofMapper::consecutive() const */
      bool consecutive () const 
      {
        return mapper1_.consecutive() && mapper2_.consecutive();
      }

      /** \copydoc Dune::DofMapper::oldOffSet( const int block ) const */
      int oldOffSet(const int block) const
      {
        const int numBlock1 = mapper1_.numBlocks();
        if ( block - numBlock1 < 0 )
          return mapper1_.oldOffSet( block );
        else
          return mapper2_.oldOffSet( block - numBlock1 ) + oldGlobalOffset_;
      }

      /** \copydoc Dune::DofMapper::OffSet( const int block ) const */
      int offSet(const int block) const
      {
        assert( globalOffset_ == mapper1_.size() );
        const int numBlock1 = mapper1_.numBlocks(); 
        if ( block - numBlock1 < 0 ) 
          return mapper1_.offSet( block );
        else
          return mapper2_.offSet( block - numBlock1 ) + globalOffset_;
      }

      /** \copydoc Dune::DofMapper::numBlocks() const */
      int numBlocks() const
      {
        return mapper1_.numBlocks() + mapper2_.numBlocks();
      }

      /** \copydoc Dune::DofMapper::resize() */
      void resize ()
      {
        oldGlobalOffset_ = globalOffset_;
        globalOffset_ = mapper1_.size();      
      }

      /** \copydoc Dune::DofMapper::compress() */
      bool compress ()
      {
        resize();
        return true;
      }

      void read_xdr ( const char *filename, int timestep )
      {
        resize();
      }

      void write_xdr ( const char *filename, int timestep )
      {}

      /** \copydoc Dune::DofMapper::insertEntity( const ElementType &element ) */
      void insertEntity ( const ElementType &element )
      {
        resize();
      }

      /** \copydoc Dune::DofMapper::removeEntity( const ElementType &element ) */
      void removeEntity ( const ElementType &element )
      {
        resize();
      }


    protected:
      DofManagerType &dm_;
      const MapperType1 &mapper1_;
      const MapperType2 &mapper2_;

      int globalOffset_;
      int oldGlobalOffset_;
    };

  }   // end namespace Fem  
} // end namespace Dune
#endif //header guards
