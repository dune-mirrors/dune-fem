#ifndef DUNE_FEM_DOFMAPPER_INDEXSETDOFMAPPER_HH
#define DUNE_FEM_DOFMAPPER_INDEXSETDOFMAPPER_HH

#include <cassert>

#include <type_traits>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/common/indexset.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/code.hh>
#include <dune/fem/space/mapper/exceptions.hh>
#include <dune/fem/space/mapper/dofmapper.hh>

namespace Dune
{

  namespace Fem
  {

    namespace __IndexSetDofMapper
    {

      // DofMapper
      // ---------

      template< class GridPart >
      class DofMapper
      {
        typedef DofMapper< GridPart > ThisType;
      public:
        typedef std::size_t SizeType;
      protected:
        struct SubEntityInfo
        {
          SubEntityInfo ()
          : numDofs( 0 )
          {}

          unsigned int codim;
          unsigned int numDofs;
          SizeType offset, oldOffset;
        };

        enum CodimType { CodimEmpty, CodimFixedSize, CodimVariableSize };
        static const int dimension = GridPart::dimension;
        typedef Dune::ReferenceElement< typename GridPart::ctype, dimension > RefElementType;
        typedef Dune::ReferenceElements< typename GridPart::ctype, dimension > RefElementsType;

        struct BuildFunctor;

        struct SubEntityFilter
        {
          SubEntityFilter(const RefElementType &refElement, int subEntity, int codim)
          : active_(dimension+1), size_(0)
          {
            for (int c=0;c<=dimension;++c)
            {
              std::vector<bool> &a = active_[c];
              a.resize( refElement.size( c ), false );
              if (c<codim) continue;
              if (c==codim) { a[subEntity]=true; ++size_; continue; }
              for (int i=0;i<refElement.size(subEntity, codim, c);++i)
              {
                a[refElement.subEntity(subEntity, codim, i, c)] = true;
                ++size_;
              }
            }
          }
          bool operator()(int i,int c) const { return active_[c][i]; }
          private:
          std::vector< std::vector< bool > > active_;
          int size_;
        };
        template< class Functor >
        struct MapFunctor;
        struct SubEntityFilterFunctor;

      public:
        typedef SizeType GlobalKeyType;

        typedef GridPart GridPartType;

        typedef typename GridPartType::template Codim< 0 >::EntityType ElementType;

        template< class CodeFactory >
        DofMapper ( const GridPartType &gridPart, const CodeFactory &codeFactory );

        // mapping for DoFs
        /** \brief map each local DoF number to a global one

            \param[in]  element  element, the DoFs belong to
            \param[in]  f        functor to call for each DoF

            The functor has to be a copyable object satisfying the following
            interface:
            \code
            struct Functor
            {
              // application operator
              void operator() ( int localDoF, int globalDoF );
            };
            \endcode

            For each DoF to be mapped, this method will call the application operator
            once.

            \note There is no guarantee on the order, in which the functor is applied.
         */
        template< class Functor >
        void mapEach ( const ElementType &element, Functor f ) const;

        void map ( const ElementType &element, std::vector< GlobalKeyType > &indices ) const;

        /** \brief fills a vector of bools with true indicating that the corresponding 
         *  local degree of freedom is attached to the subentity specified by the (c,i)
         *  pair.
         *  A local dof is attached to a subentity S if it is attached either to that 
         *  subentity or to a subentity S'<S i.e. S' has codimension greater than c 
         *  which and lies within S. For example all dofs are attached to the element itself
         *  and dofs are attached to an edge also in the case where they are attached to 
         *  the vertices of that edge.
         **/ 
        void onSubEntity( const ElementType &element, int i, int c, std::vector< bool > &indices ) const;

        unsigned int maxNumDofs () const { return maxNumDofs_; }
        unsigned int numDofs ( const ElementType &element ) const { return code( element ).numDofs(); }

        // assignment of DoFs to entities
        template< class Entity, class Functor >
        void mapEachEntityDof ( const Entity &entity, Functor f ) const;

        template< class Entity >
        void mapEntityDofs ( const Entity &entity, std::vector< GlobalKeyType > &indices ) const;

        template< class Entity >
        unsigned int numEntityDofs ( const Entity &entity ) const;

        // global information

        bool contains ( int codim ) const { return (codimType_[ codim ] != CodimEmpty); }

        bool fixedDataSize ( int codim ) const { return (codimType_[ codim ] == CodimFixedSize); }

        SizeType size () const { return size_; }

        void update ();


        /* \name AdaptiveDofMapper interface methods
         * \{
         */

        /* Compatibility methods; users expect an AdaptiveDiscreteFunction to
         * compile over spaces built on top of a LeafGridPart or LevelGridPart.
         *
         * The AdaptiveDiscreteFunction requires the block mapper (i.e. this
         * type) to be adaptive. The CodimensionMapper however is truly
         * adaptive if and only if the underlying index set is adaptive. We
         * don't want to wrap the index set as 1) it hides the actual problem
         * (don't use the AdaptiveDiscreteFunction with non-adaptive index
         * sets), and 2) other dune-fem classes may make correct use of the
         * index set's capabilities.
         */

        static constexpr bool consecutive () noexcept { return false; }

        SizeType numBlocks () const
        {
          DUNE_THROW( NotImplemented, "Method numBlocks() called on non-adaptive block mapper" );
        }

        SizeType numberOfHoles ( int ) const
        {
          DUNE_THROW( NotImplemented, "Method numberOfHoles() called on non-adaptive block mapper" );
        }

        GlobalKeyType oldIndex ( int hole, int ) const
        {
          DUNE_THROW( NotImplemented, "Method oldIndex() called on non-adaptive block mapper" );
        }

        GlobalKeyType newIndex ( int hole, int ) const
        {
          DUNE_THROW( NotImplemented, "Method newIndex() called on non-adaptive block mapper" );
        }

        SizeType oldOffSet ( int ) const
        {
          DUNE_THROW( NotImplemented, "Method oldOffSet() called on non-adaptive block mapper" );
        }

        SizeType offSet ( int ) const
        {
          DUNE_THROW( NotImplemented, "Method offSet() called on non-adaptive block mapper" );
        }

        /* \} */

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
        SizeType size_;
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

      // The functor maps all DoFs for a given entity. Intentially, it
      // is passed as argument to DofMapperCode::operator() which then
      // calls the apply()-method for each sub-entity with DoFs in turn.
      template< class GridPart >
      template< class Functor >
      struct DofMapper< GridPart >::MapFunctor
      {
        static const bool isCartesian = Dune::Capabilities::
          isCartesian< typename GridPart :: GridType > :: v ;

        MapFunctor ( const GridPart&  gridPart, const std::vector< SubEntityInfo > &subEntityInfo,
                     const ElementType &element, Functor functor )
        : gridPart_( gridPart ),
          indexSet_( gridPart_.indexSet() ),
          subEntityInfo_( subEntityInfo ),
          element_( element ),
          functor_( functor )
        {}

        // subEntity is the sub-entity number, given codim, as returned
        // by refElem.subEntity(). The iterators iterate over all DoFs
        // attached to the given sub-entity.
        template< class Iterator >
        void operator() ( unsigned int gtIndex, unsigned int subEntity, Iterator it, Iterator end )
        {
          enum { dimension = GridPart :: dimension };

          const SubEntityInfo &info = subEntityInfo_[ gtIndex ];
          const SizeType subIndex = indexSet_.subIndex( element_, subEntity, info.codim );
          SizeType index = info.offset + SizeType( info.numDofs ) * subIndex;

          const unsigned int codim = info.codim ;

          const unsigned int numDofs = info.numDofs ;
          // for non-Cartesian grids check twist if on edges with noDofs > 1
          // this should be the case for polOrder > 2.
          // Note that in 3d this only solves the twist problem up to polOrder = 3
          if( ! isCartesian && codim == dimension-1 && numDofs > 1 )
          {
            typedef typename GridPart::ctype FieldType ;
            const Dune::ReferenceElement< FieldType, dimension > &refElem
                = Dune::ReferenceElements< FieldType, dimension >::general( element_.type() );

#ifndef NDEBUG
            const int vxSize = refElem.size( subEntity, codim, dimension );
            // two vertices per edge in 2d
            assert( vxSize == 2 );
#endif
            const int vx[ 2 ] = { refElem.subEntity ( subEntity, codim, 0, dimension ),
                                  refElem.subEntity ( subEntity, codim, 1, dimension ) };

            // flip index if face is twisted
            if( gridPart_.grid().localIdSet().subId( gridEntity( element_ ), vx[ 1 ], dimension )
                < gridPart_.grid().localIdSet().subId( gridEntity( element_ ), vx[ 0 ], dimension ) )
            {
              std::vector< unsigned int > global( numDofs );
              std::vector< unsigned int > local ( numDofs );

              unsigned int count = 0 ;
              while( it != end )
              {
                global[ count ] = index++;
                local [ count ] = *(it++);
                ++count ;
              }

              unsigned int reverse = numDofs - 1;
              for( unsigned int i=0; i<numDofs; ++ i, --reverse )
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

      template< class GridPart >
      struct DofMapper< GridPart >::SubEntityFilterFunctor
      {
        static const bool isCartesian = Dune::Capabilities::
          isCartesian< typename GridPart :: GridType > :: v ;

        SubEntityFilterFunctor( const GridPart &gridPart, const std::vector< SubEntityInfo > &subEntityInfo, 
                                const ElementType &element, int i, int c,
                                std::vector<bool> &vec )
        : gridPart_( gridPart ),
          subEntityInfo_( subEntityInfo ),
          element_( element ),
          vec_(vec),
          filter_(RefElementsType::general( element.type() ), i,c)  
        {}

        template< class Iterator >
        void operator() ( unsigned int gtIndex, unsigned int subEntity, Iterator it, Iterator end )
        {
          enum { dimension = GridPart :: dimension };

          const SubEntityInfo &info = subEntityInfo_[ gtIndex ];

          const unsigned int codim = info.codim ;

          const unsigned int numDofs = info.numDofs ;

          bool active = filter_(subEntity,codim);

          // for non-Cartesian grids check twist if on codim 1 noDofs > 1
          // this should be the case for polOrder > 2
          if( ! isCartesian && codim == dimension-1 && numDofs > 1 )
          {
            typedef typename GridPart::ctype FieldType ;
            const Dune::ReferenceElement< FieldType, dimension > &refElem
                = Dune::ReferenceElements< FieldType, dimension >::general( element_.type() );

#ifndef NDEBUG
            const int vxSize = refElem.size( subEntity, codim, dimension );
            // two vertices per edge in 2d
            assert( vxSize == 2 );
#endif
            const int vx[ 2 ] = { refElem.subEntity ( subEntity, codim, 0, dimension ),
                                  refElem.subEntity ( subEntity, codim, 1, dimension ) };

            // flip index if face is twisted
            if( gridPart_.grid().localIdSet().subId( gridEntity( element_ ), vx[ 1 ], dimension )
                < gridPart_.grid().localIdSet().subId( gridEntity( element_ ), vx[ 0 ], dimension ) )
            {
              std::vector< unsigned int > global( numDofs );
              std::vector< unsigned int > local ( numDofs );

              unsigned int count = 0 ;
              while( it != end )
              {
                local [ count ] = *(it++);
                ++count ;
              }

              unsigned int reverse = numDofs - 1;
              for( unsigned int i=0; i<numDofs; ++ i, --reverse )
              {
                vec_[ local[i] ] = active;
              }

              // already did mapping, then return
              return ;
            }
          }

          // standard mapping
          while( it != end )
          {
            vec_[ *(it++) ] = active;
          }
        }

      private:
        const GridPart& gridPart_;
        const std::vector< SubEntityInfo > &subEntityInfo_;
        const ElementType &element_;
        std::vector<bool> &vec_;
        SubEntityFilter filter_;
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
              const unsigned int gtIdx = GlobalGeometryTypeIndex::index( refElement.type( i, codim ) );
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
        ::map ( const ElementType &element, std::vector< SizeType > &indices ) const
      {
        indices.resize( numDofs( element ) );
        mapEach( element, AssignFunctor< std::vector< SizeType > >( indices ) );
      }

      template< class GridPart >
      template< class Entity, class Functor >
      inline void DofMapper< GridPart >
        ::mapEachEntityDof ( const Entity &entity, Functor f ) const
      {
        const SubEntityInfo &info = subEntityInfo( entity );
        const unsigned int numDofs = info.numDofs;
        SizeType index = info.offset + numDofs * SizeType( indexSet().index( entity ) );
        for( unsigned int i = 0; i < info.numDofs; ++i )
          f( i, index++ );
      }


      template< class GridPart >
      template< class Entity >
      inline void DofMapper< GridPart >
        ::mapEntityDofs ( const Entity &entity, std::vector< SizeType > &indices ) const
      {
        indices.resize( numEntityDofs( entity ) );
        mapEachEntityDof( entity, AssignFunctor< std::vector< SizeType > >( indices ) );
      }

      template< class GridPart >
      inline void DofMapper< GridPart >
        ::onSubEntity( const ElementType &element, int i, int c, std::vector< bool > &indices ) const
      {
        indices.resize( numEntityDofs( element ) );
        code( element )( SubEntityFilterFunctor( gridPart_, subEntityInfo_, element, i,c , indices ) );
      }

      template< class GridPart >
      template< class Entity >
      inline unsigned int
      DofMapper< GridPart >
        ::numEntityDofs ( const Entity &entity ) const
      {
        return subEntityInfo( entity ).numDofs;
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
          size_ += SizeType( info.numDofs ) * SizeType( indexSet().size( *it ) );
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



      // AdaptiveDofMapper
      // -----------------

      template< class GridPart >
      class AdaptiveDofMapper
        : public DofMapper< GridPart >
      {
        typedef AdaptiveDofMapper< GridPart > ThisType;
        typedef DofMapper< GridPart > BaseType;

      protected:
        typedef typename BaseType::SubEntityInfo SubEntityInfo;

      public:
        typedef typename BaseType::GridPartType GridPartType;
        typedef typename BaseType::SizeType SizeType;

        using BaseType::update;

        template< class CodeFactory >
        AdaptiveDofMapper ( const GridPartType &gridPart, const CodeFactory &codeFactory )
          : BaseType( gridPart, codeFactory )
        {
          DofManager< typename GridPartType::GridType >::instance( gridPart_.grid() ).addIndexSet( *this );
        }

        AdaptiveDofMapper ( const ThisType & ) = delete;

        ~AdaptiveDofMapper ()
        {
          DofManager< typename GridPartType::GridType >::instance( gridPart_.grid() ).removeIndexSet( *this );
        }

        ThisType &operator= ( const ThisType & ) = delete;

        // adaptation interface

        int numBlocks () const { return blockMap_.size(); }
        SizeType offSet ( int blk ) const;
        SizeType oldOffSet ( int blk ) const;

        SizeType numberOfHoles ( int blk ) const;

        SizeType oldIndex ( SizeType hole, int blk ) const;
        SizeType newIndex ( SizeType hole, int blk ) const;

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
        using BaseType::indexSet;

        using BaseType::blockMap_;
        using BaseType::gridPart_;
        using BaseType::subEntityInfo_;
      };



      // Implementation of AdaptiveDofMapper
      // -----------------------------------

      template< class GridPart >
      inline typename AdaptiveDofMapper< GridPart >::SizeType
      AdaptiveDofMapper< GridPart >::offSet ( int blk ) const
      {
        assert( (blk >= 0) && (blk < numBlocks()) );
        const unsigned int gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
        return subEntityInfo_[ gtIdx ].offset;
      }


      template< class GridPart >
      inline typename AdaptiveDofMapper< GridPart >::SizeType
      AdaptiveDofMapper< GridPart >::oldOffSet ( int blk ) const
      {
        assert( (blk >= 0) && (blk < numBlocks()) );
        const unsigned int gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
        return subEntityInfo_[ gtIdx ].oldOffset;
      }


      template< class GridPart >
      inline typename AdaptiveDofMapper< GridPart >::SizeType
      AdaptiveDofMapper< GridPart >::numberOfHoles ( int blk ) const
      {
        assert( (blk >= 0) && (blk < numBlocks()) );
        const unsigned int gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
        const SubEntityInfo &info = subEntityInfo_[ gtIdx ];
        return SizeType( info.numDofs ) * SizeType( indexSet().numberOfHoles( blockMap_[ blk ] ) );
      }


      template< class GridPart >
      inline typename AdaptiveDofMapper< GridPart >::SizeType
      AdaptiveDofMapper< GridPart >::oldIndex ( SizeType hole, int blk ) const
      {
        assert( (hole >= 0) && (hole < numberOfHoles( blk )) );
        const unsigned int gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
        const SubEntityInfo &info = subEntityInfo_[ gtIdx ];
        const unsigned int numDofs = info.numDofs;
        const SizeType index = indexSet().oldIndex( hole / numDofs, blockMap_[ blk ] );
        return info.offset + numDofs * index + (hole % numDofs);
      }


      template< class GridPart >
      inline typename AdaptiveDofMapper< GridPart >::SizeType
      AdaptiveDofMapper< GridPart >::newIndex ( SizeType hole, int blk ) const
      {
        assert( (hole >= 0) && (hole < numberOfHoles( blk )) );
        const unsigned int gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
        const SubEntityInfo &info = subEntityInfo_[ gtIdx ];
        const unsigned int numDofs = info.numDofs;
        const SizeType index = indexSet().newIndex( hole / numDofs, blockMap_[ blk ] );
        return info.offset + numDofs * index + (hole % numDofs);
      }



      // Implementation
      // --------------

      template< class GridPart, bool adaptive = Capabilities::isAdaptiveIndexSet< typename GridPart::IndexSetType >::v >
      struct Implementation
      {
        typedef typename std::conditional< adaptive, AdaptiveDofMapper< GridPart >, DofMapper< GridPart > >::type Type;
      };

    } // namespace __IndexSetDofMapper



    // IndexSetDofMapper
    // -----------------

    template< class GridPart >
    class IndexSetDofMapper
      : public __IndexSetDofMapper::template Implementation< GridPart >::Type
    {
      typedef typename __IndexSetDofMapper::template Implementation< GridPart >::Type BaseType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      template< class CodeFactory >
      IndexSetDofMapper ( const GridPartType &gridPart, const CodeFactory &codeFactory )
        : BaseType( gridPart, codeFactory )
      {}
    };


    // Capabilities
    // ------------

    namespace Capabilities
    {
      // isAdaptiveDofMapper
      // -------------------

      template< class GridPart >
      struct isAdaptiveDofMapper< IndexSetDofMapper< GridPart > >
      {
        static const bool v = Capabilities::isAdaptiveIndexSet< typename GridPart::IndexSetType >::v;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_DOFMAPPER_INDEXSETDOFMAPPER_HH
