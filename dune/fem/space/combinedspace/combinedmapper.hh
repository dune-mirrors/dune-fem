#ifndef DUNE_FEM_COMBINEDSPACEMAPPER_HH
#define DUNE_FEM_COMBINEDSPACEMAPPER_HH

#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/common/dofmanager.hh>


namespace Dune
{

  namespace Fem
  {

    //! forward declaration
    template< class Grid,  class BlockMapper1, int blockSize1, class BlockMapper2, int blockSize2 >
    class CombinedSpaceMapper;

    //! Traits 
    template< class Grid,  class BlockMapper1, int blockSize1, class BlockMapper2, int blockSize2 >
    struct CombinedSpaceMapperTraits
    {
      // we still need entities of codimension 0 here 
      typedef typename BlockMapper1 :: ElementType ElementType;

      static const int polynomialOrder1 = BlockMapper1 :: polynomialOrder;
      static const int polynomialOrder2 = BlockMapper2 :: polynomialOrder;
      static const int polynomialOrder = ( polynomialOrder1 > polynomialOrder2 ) ? polynomialOrder1 : polynomialOrder2;

      typedef NonBlockMapper< BlockMapper1, blockSize1 > MapperType1; 
      typedef NonBlockMapper< BlockMapper2, blockSize2 > MapperType2; 

      typedef std::size_t SizeType;

      typedef CombinedSpaceMapper< Grid, BlockMapper1, blockSize1, BlockMapper2, blockSize2 >  DofMapperType;
    };

    template< class Grid,  class BlockMapper1, int blockSize1, class BlockMapper2, int blockSize2 >
    class CombinedSpaceMapper
    : public AdaptiveDofMapper< CombinedSpaceMapperTraits< Grid, BlockMapper1, blockSize1, BlockMapper2, blockSize2 > >
    {
      typedef AdaptiveDofMapper< CombinedSpaceMapperTraits< Grid, BlockMapper1, blockSize1, BlockMapper2, blockSize2 > > BaseType;
      typedef CombinedSpaceMapper< Grid, BlockMapper1, blockSize1, BlockMapper2, blockSize2 >  ThisType;
      typedef Grid GridType;

    public:
      typedef typename BaseType::Traits Traits;

    protected:
      typedef typename Traits::MapperType1 MapperType1;
      typedef typename Traits::MapperType2 MapperType2;

      template< class Functor >
      struct FunctorWrapper
      {
        FunctorWrapper ( Functor functor, int localOffset, int globalOffset )
        : functor_( functor ),         
          localOffset_(  localOffset ),
          globalOffset_( globalOffset )
        {}

        template< class GlobalKey >
        void operator() ( int localDof, const GlobalKey &globalKey ) 
        {
          functor_( localDof + localOffset_, globalKey + globalOffset_ );
        }

        template< class GlobalKey >
        void operator() ( const GlobalKey &globalKey ) 
        {
          functor_( globalKey + globalOffset_ );
        }

      private:
        Functor functor_;
        int localOffset_;
        int globalOffset_;
      };

      public:
      typedef typename BaseType :: ElementType ElementType;

      //! order of the Lagrange polynoms
      static const int polynomialOrder1 = Traits :: polynomialOrder1;
      static const int polynomialOrder2 = Traits :: polynomialOrder2;
      static const int polynomialOrder = Traits :: polynomialOrder;

      //! type of the DoF manager
      typedef DofManager< GridType > DofManagerType;

      public:
      //! constructor
      CombinedSpaceMapper( const GridType &grid, BlockMapper1 &blockMapper1, BlockMapper2& blockMapper2 )
      :  dm_( DofManagerType :: instance( grid ) ), 
         mapper1_( blockMapper1 ),
         mapper2_( blockMapper2 ),
         globalOffset_( mapper1_.size() ),
         oldGlobalOffset_( -1 )
      {
        dm_.addIndexSet( *this );
      }

      ~CombinedSpaceMapper() 
      {
        dm_.removeIndexSet( *this );
      }
      
      //! copy constructor
      CombinedSpaceMapper( const ThisType &other )
        : dm_( other.dm_ ),
          mapper1_( other.mapper1_ ),
          mapper2_( other.mapper2_ ),
          globalOffset_( other.globalOffset_ ),
          oldGlobalOffset_( other.oldGlobalOffset_ )
      {
        dm_.addIndexSet( *this );
      }

      /** \copydoc Dune::DofMapper::size const */
      int size () const 
      {
        return mapper1_.size() + mapper2_.size();
      }

      /** \copydoc Dune::DofMapper::contains(const int codim) const */
      bool contains ( const int codim ) const    
      {
        return ( mapper1_.contains( codim ) || mapper2_.contains( codim ) );
      }

      /** \copydoc Dune::DofMapper::mapToGlobal(const ElementType &entity, const int localDof) const */
      int mapToGlobal ( const ElementType &entity, const int localDof ) const
      {
        assert( mapper1_.size() == globalOffset_ );
        const int localOffset = mapper1_.numDofs(entity);

        int index;

        if( localDof - localOffset  < 0 )
          index = mapper1_.mapToGlobal( entity, localDof );
        else
          index = mapper2_.mapToGlobal( entity, localDof - localOffset ) + globalOffset_;
          
        assert( (0 <= index) && (index < size()) );
        return index;
      }
      
      /** \copydoc Dune::DofMapper::mapEachEntityDof(const Entity &entity,Functor f) const */
      template< class Entity, class Functor > 
      void mapEachEntityDof ( const Entity &entity, Functor f ) const
      {
        mapper1_.mapEachEntityDof( entity, FunctorWrapper< Functor >( f, 0, 0 ) );
        mapper2_.mapEachEntityDof( entity, FunctorWrapper< Functor >( f, mapper1_.numEntityDofs( entity ), globalOffset_) );
      }
      
      /** \copydoc Dune::DofMapper::maxNumDofs const */
      int maxNumDofs () const
      {
        return mapper1_.maxNumDofs() + mapper2_.maxNumDofs();
      }


      /** \copydoc Dune::DofMapper::mapEach(const ElementType &element, Functor f) const */
      template< class Functor >
      void mapEach ( const ElementType &element, Functor f ) const
      {
        mapper1_.mapEach( element, FunctorWrapper< Functor > (f, 0, 0) );
        mapper2_.mapEach( element, FunctorWrapper< Functor > (f, mapper1_.numDofs( element ), globalOffset_) );
      }

      /** \copydoc Dune::DofMapper::numDofs(const ElementType &element) const */
      int numDofs ( const ElementType &element ) const
      {
        int nDofs = mapper1_.numDofs( element ) + mapper2_.numDofs( element );
  //      assert( (nDofs >=0) && (nDofs < size()) );
        return nDofs;
      }

      /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
      template< class Entity >
      int numEntityDofs ( const Entity &entity ) const
      {
        return mapper1_.numEntityDofs( entity ) + mapper2_.numEntityDofs( entity );
      }

      /** \copydoc Dune::DofMapper::fixedDataSize(int codim) const */
      bool fixedDataSize ( int codim ) const
      {
        return mapper1_.fixedDataSize( codim ) && mapper2_.fixedDataSize( codim );
      }

      /** \copydoc Dune::DofMapper::numberOfHoles(const int block) const */
      int numberOfHoles(const int block) const 
      {
        const int numBlock1 = mapper1_.numBlocks(); 
        if ( block - numBlock1 < 0 ) 
          return mapper1_.numberOfHoles( block );
        else
          return mapper2_.numberOfHoles( block - numBlock1 );
      }
      
      /** \copydoc Dune::DofMapper::oldIndex(const int hole, const int block) const */
      int oldIndex (const int hole, const int block) const 
      { 
        const int numBlock1 = mapper1_.numBlocks(); 
        if ( block - numBlock1 < 0 ) 
          return mapper1_.oldIndex( hole, block );
        else
          return mapper2_.oldIndex( hole, block - numBlock1 );
      }
        
      /** \copydoc Dune::DofMapper::newIndex(const int hole, const int block) const */
      int newIndex (const int hole, const int block) const 
      { 
        const int numBlock1 = mapper1_.numBlocks(); 
        if ( block - numBlock1 < 0 ) 
          return mapper1_.newIndex( hole, block );
        else
          return mapper2_.newIndex( hole, block - numBlock1 );
      }

      /** \copydoc Dune::DofMapper::consecutive const */
      bool consecutive () const 
      {
        return mapper1_.consecutive() && mapper2_.consecutive();
      }

      /** \copydoc Dune::DofMapper::oldOffSet(const int block) const */
      int oldOffSet(const int block) const
      {
        const int numBlock1 = mapper1_.numBlocks();
        if ( block - numBlock1 < 0 )
          return mapper1_.oldOffSet( block );
        else
          return mapper2_.oldOffSet( block - numBlock1 ) + oldGlobalOffset_;
      }

      /** \copydoc Dune::DofMapper::offSet(const int block) const */
      int offSet(const int block) const
      {
        assert( globalOffset_ == mapper1_.size() );
        const int numBlock1 = mapper1_.numBlocks(); 
        if ( block - numBlock1 < 0 ) 
          return mapper1_.offSet( block );
        else
          return mapper2_.offSet( block - numBlock1 ) + globalOffset_;
      }

      int numBlocks() const
      {
        return mapper1_.numBlocks() + mapper2_.numBlocks();
      }

      void resize ()
      {
        oldGlobalOffset_ = globalOffset_;
        globalOffset_ = mapper1_.size();      
      }

      bool compress ()
      {
        resize();
        return true;
      }

      void backup() const
      {}

      void restore() 
      {
        resize();
      }

      template< class IOStream >
      void read ( IOStream &in )
      {
        resize();
      }

      template< class IOStream >
      void write ( IOStream &out ) const
      {}

      template< class Entity >
      void insertEntity ( const Entity &)
      {
        resize();
      }

      template< class Entity >
      void removeEntity ( const Entity & )
      {
        resize();
      }


    protected:
      DofManagerType &dm_;
      MapperType1 mapper1_;
      MapperType2 mapper2_;

      int globalOffset_;
      int oldGlobalOffset_;
    };

  }   // namespace Fem  

} // namespace Dune
 
#endif // #ifndef DUNE_FEM_COMBINEDSPACEMAPPER_HH
