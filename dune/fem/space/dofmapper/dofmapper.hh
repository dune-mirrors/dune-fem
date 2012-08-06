#ifndef DUNE_FEM_DOFMAPPER_DOFMAPPER_HH
#define DUNE_FEM_DOFMAPPER_DOFMAPPER_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/dofmapper/code.hh>
#include <dune/fem/space/dofmapper/exceptions.hh>
#include <dune/fem/space/mapper/dofmapper.hh>

namespace Dune
{

  namespace Fem
  {

    // DofMapper
    // ---------

    template< class GridPart >
    class DofMapper
    {
      typedef DofMapper< GridPart > ThisType;

      struct SubEntityInfo
      {
        SubEntityInfo ()
        : numDofs( 0 )
        {}

        unsigned int codim;
        unsigned int numDofs;
        std::size_t offset, oldOffset;
      };

      enum CodimType { CodimEmpty, CodimFixedSize, CodimVariableSize };

      struct BuildFunctor;

      template< class Functor >
      struct MapFunctor;

      static const int dimension = GridPart::dimension;

      typedef GenericReferenceElement< typename GridPart::ctype, dimension > RefElementType;
      typedef GenericReferenceElements< typename GridPart::ctype, dimension > RefElementsType;


      // forbid copying and assignment
      DofMapper ( const ThisType & );
      const ThisType &operator= ( const ThisType & );

    public:
      typedef GridPart GridPartType;

      typedef typename GridPartType::template Codim< 0 >::EntityType ElementType;

      typedef DefaultDofMapIterator< ElementType, ThisType > DofMapIteratorType;

      template< class CodeFactory >
      DofMapper ( const GridPartType &gridPart, const CodeFactory &codeFactory );

      ~DofMapper ();


      // mapping for DoFs
      /** \brief map each local DoF number to a global one
       *
       *  \param[in]  element  element, the DoFs belong to
       *  \param[in]  f        functor to call for each DoF
       *
       *  The functor has to be a copyable object satisfying the following
       *  interface:
       *  \code
       *  struct Functor
       *  {
       *    // application operator
       *    void operator() ( int localDoF, int globalDoF );
       *  };
       *  \endcode
       *
       *  For each DoF to be mapped, this method will call the application operator
       *  once.
       *  
       *  \note There is no guarantee on the order, in which the functor is applied.
       */
      template< class Functor >
      void mapEach ( const ElementType &element, Functor f ) const;

      void map ( const ElementType &element, std::vector< std::size_t > &indices ) const;
      std::size_t mapToGlobal ( const ElementType &element, int i ) const;

      DofMapIteratorType begin ( const ElementType &element ) const;
      DofMapIteratorType end ( const ElementType &element ) const;

      std::size_t maxNumDofs () const { return maxNumDofs_; }
      std::size_t numDofs ( const ElementType &element ) const { return code( element ).numDofs(); }


      // assignment of DoFs to entities

      template< class Entity, class Functor >
      void mapEachEntityDof ( const Entity &entity, Functor f ) const;

      template< class Entity >
      void mapEntityDofs ( const Entity &entity, std::vector< std::size_t > &indices ) const;

      template< class Entity >
      std::size_t mapEntityDofToGlobal ( const Entity &entity, int i ) const;

      template< class Entity >
      std::size_t numEntityDofs ( const Entity &entity ) const;


      // global information

      bool contains ( int codim ) const { return (codimType_[ codim ] != CodimEmpty); }

      bool fixedDataSize ( int codim ) const { return (codimType_[ codim ] == CodimFixedSize); }

      std::size_t size () const { return size_; }


      // adaptation interface

      int numBlocks () const { return blockMap_.size(); }

      std::size_t offSet ( int blk ) const;
      std::size_t oldOffSet ( int blk ) const;

      std::size_t numberOfHoles ( int blk ) const;

      std::size_t oldIndex ( std::size_t hole, int blk ) const;
      std::size_t newIndex ( std::size_t hole, int blk ) const;

      void update ();


      // adaptation methods (as for index sets)

      bool consecutive () const { return true; }

      template< class Entity >
      void insertEntity ( const Entity &entity ) { update(); }

      template< class Entity >
      void removeEntity ( const Entity &entity ) {}

      void resize () { update(); }

      bool compress () { update(); return true; }

      template <class StreamTraits> 
      void write( OutStreamInterface< StreamTraits >& out ) const {} 

      template <class StreamTraits> 
      void read( InStreamInterface< StreamTraits >& in )  
      {
        update();
      } 

      void backup () const {} 
      void restore () {}

    protected:
      typedef typename GridPartType::IndexSetType IndexSetType;
      typedef std::vector< GeometryType > BlockMapType;

      const DofMapperCode &code ( const GeometryType &gt ) const;
      const DofMapperCode &code ( const ElementType &element ) const { return code( element.type() ); }

      template< class Entity >
      const SubEntityInfo &subEntityInfo ( const Entity &entity ) const;

      const IndexSetType &indexSet () const { return gridPart_.indexSet(); }

      const GridPartType &gridPart_;
      std::vector< DofMapperCode > code_;
      unsigned int maxNumDofs_;
      std::size_t size_;
      std::vector< SubEntityInfo > subEntityInfo_;
      BlockMapType blockMap_;
      CodimType codimType_[ dimension+1 ];
    };



    // DofMapper::BuildFunctor
    // -----------------------

    template< class GridPart >
    struct DofMapper< GridPart >::BuildFunctor
    {
      explicit BuildFunctor ( std::vector< SubEntityInfo > &subEntityInfo )
      : subEntityInfo_( subEntityInfo )
      {}

      template< class Iterator >
      void operator() ( unsigned int gtIndex, unsigned int subEntity, Iterator it, Iterator end )
      {
        SubEntityInfo &info = subEntityInfo_[ gtIndex ];
        const unsigned int numDofs = end - it;
        if( info.numDofs == 0 )
          info.numDofs = numDofs;
        else if( info.numDofs != numDofs )
          DUNE_THROW( DofMapperError, "Inconsistent number of DoFs on subEntity (codim = " << info.codim << ")." );
      }

    private:
      std::vector< SubEntityInfo > &subEntityInfo_;
    };



    // DofMapper::MapFunctor
    // ---------------------

    template< class GridPart >
    template< class Functor >
    struct DofMapper< GridPart >::MapFunctor
    {
      static const bool isCartesian = Capabilities ::
        isCartesian< typename GridPart :: GridType > :: v ;

      MapFunctor ( const GridPart&  gridPart, const std::vector< SubEntityInfo > &subEntityInfo,
                   const ElementType &element, Functor functor )
      : gridPart_( gridPart ),
        indexSet_( gridPart_.indexSet() ),
        subEntityInfo_( subEntityInfo ),
        element_( element ),
        functor_( functor )
      {}

      template< class Iterator >
      void operator() ( unsigned int gtIndex, unsigned int subEntity, Iterator it, Iterator end )
      {
        enum { dimension = GridPart :: dimension };

        const SubEntityInfo &info = subEntityInfo_[ gtIndex ];
        const std::size_t subIndex = indexSet_.subIndex( element_, subEntity, info.codim );
        std::size_t index = info.offset + std::size_t( info.numDofs ) * subIndex;

        const std::size_t codim = info.codim ;

        const size_t numDofs = info.numDofs ;
        // for non-Cartesian grids check twist if on codim 1 noDofs > 1 
        // this should be the case for polOrder > 2  
        if( ! isCartesian && dimension == 2 && codim == 1 && numDofs > 1 ) 
        {
          typedef typename GridPart::ctype FieldType ;
          const GenericReferenceElement< FieldType, dimension > &refElem
              = GenericReferenceElements< FieldType, dimension >::general( element_.type() );

#ifndef NDEBUG
          const int vxSize = refElem.size( subEntity, codim, dimension );
          // two vertices per edge in 2d 
          assert( vxSize == 2 );
#endif
          const int vx[ 2 ] = { refElem.subEntity ( subEntity, codim, 0, dimension ),
                                refElem.subEntity ( subEntity, codim, 1, dimension ) };

          // flip index if face is twisted 
          if( gridPart_.grid().localIdSet().subId( gridEntity( element_ ), vx[ 0 ], dimension )
              > gridPart_.grid().localIdSet().subId( gridEntity( element_ ), vx[ 1 ], dimension ) )
          {
            std::size_t global[ numDofs ];
            std::size_t local[ numDofs ];

            std::size_t count = 0 ;
            while( it != end )
            {
              global[ count ] = index++;
              local [ count ] = *(it++);
              ++count ;
            }

            size_t reverse = numDofs - 1;
            for( std::size_t i=0; i<numDofs; ++ i, --reverse )
            {
              functor_( local[ i ], global[ reverse ] ); 
            }

            // already did mapping, then return  
            return ;
          }
        }

        // standard mapping  
        while( it != end )
        {
          functor_( *(it++), index++ );
        }
      }

    private:
      const GridPart& gridPart_;
      const IndexSetType &indexSet_;
      const std::vector< SubEntityInfo > &subEntityInfo_;
      const ElementType &element_;
      Functor functor_;
    };



    // Implementation of DofMapper
    // ---------------------------

    template< class GridPart >
    const int DofMapper< GridPart >::dimension;

    template< class GridPart >
    template< class CodeFactory >
    inline DofMapper< GridPart >
      ::DofMapper ( const GridPartType &gridPart, const CodeFactory &codeFactory )
    : gridPart_( gridPart ),
      code_( LocalGeometryTypeIndex::size( dimension ) ),
      maxNumDofs_( 0 ),
      subEntityInfo_( GlobalGeometryTypeIndex::size( dimension ) )
    {
      std::vector< GeometryType > gt( GlobalGeometryTypeIndex::size( dimension ) );

      const typename RefElementsType::Iterator end = RefElementsType::end();
      for( typename RefElementsType::Iterator it = RefElementsType::begin(); it != end; ++it )
      {
        const RefElementType &refElement = *it;

        for( int codim = 0; codim <= dimension; ++codim )
        {
          for( int i = 0; i < refElement.size( codim ); ++i )
          {
            const std::size_t gtIdx = GlobalGeometryTypeIndex::index( refElement.type( i, codim ) );
            gt[ gtIdx ] = refElement.type( i, codim );
            subEntityInfo_[ gtIdx ].codim = codim;
          }
        }

        DofMapperCode &code = code_[ LocalGeometryTypeIndex::index( refElement.type() ) ];
        code = codeFactory( refElement );
        maxNumDofs_ = std::max( code.numDofs(), maxNumDofs_ );
        code( BuildFunctor( subEntityInfo_ ) );
      }

      for( int codim = 0; codim <= dimension; ++codim )
        codimType_[ codim ] = CodimEmpty;

      unsigned int codimDofs[ dimension+1 ];
      for( unsigned int i = 0; i < subEntityInfo_.size(); ++i )
      {
        const SubEntityInfo &info = subEntityInfo_[ i ];
        if( info.numDofs == 0 )
          continue;

        if( codimType_[ info.codim ] == CodimEmpty )
          codimType_[ info.codim ] = CodimFixedSize;
        else if( codimDofs[ info.codim ] != info.numDofs )
          codimType_[ info.codim ] = CodimVariableSize;

        codimDofs[ info.codim ] = info.numDofs;
        blockMap_.push_back( gt[ i ] );
      }

      update();
      DofManager< typename GridPartType::GridType >::instance( gridPart_.grid() ).addIndexSet( *this );
    }


    template< class GridPart >
    inline DofMapper< GridPart >::~DofMapper ()
    {
      DofManager< typename GridPartType::GridType >::instance( gridPart_.grid() ).removeIndexSet( *this );
    }


    template< class GridPart >
    template< class Functor >
    inline void DofMapper< GridPart >
      ::mapEach ( const ElementType &element, Functor f ) const
    {
      code( element )( MapFunctor< Functor >( gridPart_, subEntityInfo_, element, f ) );
    }


    template< class GridPart >
    inline void DofMapper< GridPart >
      ::map ( const ElementType &element, std::vector< std::size_t > &indices ) const
    {
      indices.resize( numDofs( element ) );
      mapEach( element, AssignFunctor< std::vector< std::size_t > >( indices ) );
    }


    template< class GridPart >
    inline std::size_t DofMapper< GridPart >
      ::mapToGlobal ( const ElementType &element, int i ) const
    {
      std::size_t index = 0;
      mapEach( element, AssignSingleFunctor< std::size_t >( i, index ) );
      return index;
    }


    template< class GridPart >
    inline typename DofMapper< GridPart >::DofMapIteratorType
    DofMapper< GridPart >::begin ( const ElementType &element ) const
    {
      return DofMapIteratorType( DofMapIteratorType::beginIterator, element, *this );
    }


    template< class GridPart >
    inline typename DofMapper< GridPart >::DofMapIteratorType
    DofMapper< GridPart >::end ( const ElementType &element ) const
    {
      return DofMapIteratorType( DofMapIteratorType::endIterator, element, *this );
    }


    template< class GridPart >
    template< class Entity, class Functor >
    inline void DofMapper< GridPart >
      ::mapEachEntityDof ( const Entity &entity, Functor f ) const
    {
      const SubEntityInfo &info = subEntityInfo( entity );
      std::size_t numDofs = info.numDofs;
      std::size_t index = info.offset + numDofs * std::size_t( indexSet().index( entity ) );
      for( std::size_t i = 0; i < info.numDofs; ++i )
        f( i, index++ );
    }


    template< class GridPart >
    template< class Entity >
    inline void DofMapper< GridPart >
      ::mapEntityDofs ( const Entity &entity, std::vector< std::size_t > &indices ) const
    {
      indices.resize( numEntityDofs( entity ) );
      mapEachEntityDof( entity, AssignFunctor< std::vector< std::size_t > >( indices ) );
    }


    template< class GridPart >
    template< class Entity >
    inline std::size_t DofMapper< GridPart >
      ::mapEntityDofToGlobal ( const Entity &entity, const int i ) const
    {
      assert( (i >= 0) && (i < int( numEntityDofs( entity ) ) ) );
      const SubEntityInfo &info = subEntityInfo( entity );
      return info.offset + info.numDofs*indexSet().index( entity ) + i;
    }


    template< class GridPart >
    template< class Entity >
    inline std::size_t DofMapper< GridPart >
      ::numEntityDofs ( const Entity &entity ) const
    {
      return subEntityInfo( entity ).numDofs;
    }


    template< class GridPart >
    inline std::size_t DofMapper< GridPart >::offSet ( int blk ) const
    {
      assert( (blk >= 0) && (blk < numBlocks()) );
      const std::size_t gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
      return subEntityInfo_[ gtIdx ].offset;
    }


    template< class GridPart >
    inline std::size_t DofMapper< GridPart >::oldOffSet ( int blk ) const
    {
      assert( (blk >= 0) && (blk < numBlocks()) );
      const std::size_t gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
      return subEntityInfo_[ gtIdx ].oldOffset;
    }


    template< class GridPart >
    inline std::size_t DofMapper< GridPart >::numberOfHoles ( int blk ) const
    {
      assert( (blk >= 0) && (blk < numBlocks()) );
      const std::size_t gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
      const SubEntityInfo &info = subEntityInfo_[ gtIdx ];
      return std::size_t( info.numDofs ) * std::size_t( indexSet().numberOfHoles( info.codim ) );
    }


    template< class GridPart >
    inline std::size_t DofMapper< GridPart >::oldIndex ( std::size_t hole, int blk ) const
    {
      assert( (hole >= 0) && (hole < numberOfHoles( blk )) );
      const std::size_t gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
      const SubEntityInfo &info = subEntityInfo_[ gtIdx ];
      const std::size_t numDofs = info.numDofs;
      const std::size_t index = indexSet().oldIndex( hole / numDofs, info.codim );
      return info.offset + numDofs * index + (hole % numDofs);
    }


    template< class GridPart >
    inline std::size_t DofMapper< GridPart >::newIndex ( std::size_t hole, int blk ) const
    {
      assert( (hole >= 0) && (hole < numberOfHoles( blk )) );
      const std::size_t gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
      const SubEntityInfo &info = subEntityInfo_[ gtIdx ];
      const std::size_t numDofs = info.numDofs;
      const std::size_t index = indexSet().newIndex( hole / numDofs, info.codim );
      return info.offset + numDofs * index + (hole % numDofs);
    }


    template< class GridPart >
    inline void DofMapper< GridPart >::update ()
    {
      size_ = 0;
      for( typename BlockMapType::const_iterator it = blockMap_.begin(); it != blockMap_.end(); ++it )
      {
        SubEntityInfo &info = subEntityInfo_[ GlobalGeometryTypeIndex::index( *it ) ];
        info.oldOffset = info.offset;
        info.offset = size_;
        size_ += std::size_t( info.numDofs ) * std::size_t( indexSet().size( *it ) );
      }
    }


    template< class GridPart >
    inline const DofMapperCode &DofMapper< GridPart >
      ::code ( const GeometryType &gt ) const
    {
      return code_[ LocalGeometryTypeIndex::index( gt ) ];
    }


    template< class GridPart >
    template< class Entity >
    inline const typename DofMapper< GridPart >::SubEntityInfo &
    DofMapper< GridPart >::subEntityInfo ( const Entity &entity ) const
    {
      return subEntityInfo_[ GlobalGeometryTypeIndex::index( entity.type() ) ];
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DOFMAPPER_DOFMAPPER_HH
