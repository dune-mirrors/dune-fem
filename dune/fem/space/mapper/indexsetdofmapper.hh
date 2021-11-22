#ifndef DUNE_FEM_DOFMAPPER_INDEXSETDOFMAPPER_HH
#define DUNE_FEM_DOFMAPPER_INDEXSETDOFMAPPER_HH

#include <cassert>

#include <type_traits>
#include <utility>

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

    // DefaultLocalDofMapping
    // ----------------------

    template< class GridPart >
    class DefaultLocalDofMapping
    {
      struct Mapping
      {
        template< class Iterator, class Functor >
        void operator() ( std::size_t index, unsigned int numDofs, Iterator begin, Iterator end, Functor functor ) const
        {
          while( begin != end )
            functor( *(begin++), index++ );
        }
      };

    public:
      DefaultLocalDofMapping () {}
      DefaultLocalDofMapping ( const GridPart & ) {}

      Mapping operator() ( const typename GridPart::template Codim< 0 >::EntityType &element, unsigned int subEntity, unsigned int codim ) const { return {}; }
    };



    namespace __IndexSetDofMapper
    {

      // DofMapperBase
      // -------------

      template< class GridPart, class LocalDofMapping >
      class DofMapperBase
      {
        typedef DofMapperBase< GridPart, LocalDofMapping > ThisType;

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
        typedef Dune::ReferenceElements< typename GridPart::ctype, dimension > RefElementsType;
        typedef typename RefElementsType::ReferenceElement RefElementType;

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

      public:
        typedef SizeType GlobalKeyType;

        typedef GridPart GridPartType;
        typedef LocalDofMapping LocalDofMappingType;
        typedef typename GridPart::GridType GridType;

        typedef typename GridPartType::template Codim< 0 >::EntityType ElementType;

        template< class CodeFactory >
        DofMapperBase ( const GridPartType &gridPart, LocalDofMappingType localDofMapping, const CodeFactory &codeFactory );

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
         *  and lies within S. For example all dofs are attached to the element itself
         *  and dofs attached to a vertex of an edge are also attached to that edge.
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

        //! update mapper offsets
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

        template< class Entity >
        void insertEntity ( const Entity &entity )
        {
          DUNE_THROW( NotImplemented, "Method insertEntity(entity) called on non-adaptive block mapper" );
        }
        template< class Entity >
        void removeEntity ( const Entity &entity )
        {
          DUNE_THROW( NotImplemented, "Method removeEntity(entity) called on non-adaptive block mapper" );
        }

        void resize () { update(); }
        bool compress () { update(); return true;}
        void backup () const {}
        void restore () {}
        template <class StreamTraits>
        void write( OutStreamInterface< StreamTraits >& out ) const {}
        template <class StreamTraits>
        void read( InStreamInterface< StreamTraits >& in ) { update(); }

        /* \} */

      protected:
        //! submit request for codimensions used to index set
        void requestCodimensions ();

        typedef typename GridPartType::IndexSetType IndexSetType;
        typedef std::vector< GeometryType > BlockMapType;

        const DofMapperCode &code ( const GeometryType &gt ) const;
        const DofMapperCode &code ( const ElementType &element ) const { return code( element.type() ); }

        template< class Entity >
        const SubEntityInfo &subEntityInfo ( const Entity &entity ) const;

        const IndexSetType &indexSet () const { return indexSet_; }

        const IndexSetType &indexSet_;
        LocalDofMapping localDofMapping_;
        std::vector< DofMapperCode > code_;
        unsigned int maxNumDofs_;
        SizeType size_;
        std::vector< SubEntityInfo > subEntityInfo_;
        BlockMapType blockMap_;
        CodimType codimType_[ dimension+1 ];
      };



      // DofMapper::BuildFunctor
      // -----------------------

      template< class GridPart, class LocalDofMapping >
      struct DofMapperBase< GridPart, LocalDofMapping >::BuildFunctor
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




      // Implementation of DofMapper
      // ---------------------------

      template< class GridPart, class LocalDofMapping >
      const int DofMapperBase< GridPart, LocalDofMapping >::dimension;

      template< class GridPart, class LocalDofMapping >
      template< class CodeFactory >
      inline DofMapperBase< GridPart, LocalDofMapping >
        ::DofMapperBase ( const GridPartType &gridPart, LocalDofMappingType localDofMapping, const CodeFactory &codeFactory )
      // NOTE: Don't store gridPart in this class since the lifetime of gridPart
      // might be shorter than the lifetime of this class. The lifetime of
      // indexSet is guaranteed to be longer, so storage of that class is fine
      : indexSet_( gridPart.indexSet() ),
        localDofMapping_( std::move( localDofMapping ) ),
        code_( LocalGeometryTypeIndex::size( dimension ) ),
        maxNumDofs_( 0 ),
        subEntityInfo_( GlobalGeometryTypeIndex::size( dimension ) )
      {
        std::vector< GeometryType > gt( GlobalGeometryTypeIndex::size( dimension ) );

        const typename RefElementsType::Iterator end = RefElementsType::end();
        for( typename RefElementsType::Iterator it = RefElementsType::begin(); it != end; ++it )
        {
          const RefElementType refElement = *it;

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
          {
            continue;
          }

          // see commit message f86ab6e96a27fdecfa82de43fe9099f01e240e1b
          // Note: hasSingleGeometryType does not exist on all IndexSets
          static const bool hasSingleGeometryType = Dune::Capabilities::hasSingleGeometryType< typename GridPartType::GridType > :: v ;
          const auto & geomTypes = indexSet().types(info.codim);

          if (hasSingleGeometryType && geomTypes[0] != gt[i])
          {
            continue;
          }

          if( codimType_[ info.codim ] == CodimEmpty )
            codimType_[ info.codim ] = CodimFixedSize;
          else if( codimDofs[ info.codim ] != info.numDofs )
            codimType_[ info.codim ] = CodimVariableSize;

          codimDofs[ info.codim ] = info.numDofs;
          blockMap_.push_back( gt[ i ] );
        }

        // submit request for codimensions to index set
        requestCodimensions ();

        // update offsets
        update();
      }


      template< class GridPart, class LocalDofMapping >
      template< class Functor >
      inline void DofMapperBase< GridPart, LocalDofMapping >
        ::mapEach ( const ElementType &element, Functor f ) const
      {
        const auto &idxSet = indexSet();

        code( element )( [ this, &idxSet, &element, f ] ( unsigned int gtIndex, unsigned int subEntity, auto begin, auto end ) {
            const SubEntityInfo &info = subEntityInfo_[ gtIndex ];
            const SizeType subIndex = idxSet.subIndex( element, subEntity, info.codim );
            SizeType index = info.offset + SizeType( info.numDofs ) * subIndex;
            localDofMapping_( element, subEntity, info.codim )( index, info.numDofs, begin, end, f );
          } );
      }


      template< class GridPart, class LocalDofMapping >
      inline void DofMapperBase< GridPart, LocalDofMapping >
        ::map ( const ElementType &element, std::vector< SizeType > &indices ) const
      {
        indices.resize( numDofs( element ) );
        mapEach( element, AssignFunctor< std::vector< SizeType > >( indices ) );
      }

      template< class GridPart, class LocalDofMapping >
      template< class Entity, class Functor >
      inline void DofMapperBase< GridPart, LocalDofMapping >
        ::mapEachEntityDof ( const Entity &entity, Functor f ) const
      {
        const SubEntityInfo &info = subEntityInfo( entity );
        const unsigned int numDofs = info.numDofs;
        SizeType index = info.offset + numDofs * SizeType( indexSet().index( entity ) );
        for( unsigned int i = 0; i < info.numDofs; ++i )
          f( i, index++ );
      }


      template< class GridPart, class LocalDofMapping >
      template< class Entity >
      inline void DofMapperBase< GridPart, LocalDofMapping >
        ::mapEntityDofs ( const Entity &entity, std::vector< SizeType > &indices ) const
      {
        indices.resize( numEntityDofs( entity ) );
        mapEachEntityDof( entity, AssignFunctor< std::vector< SizeType > >( indices ) );
      }

      template< class GridPart, class LocalDofMapping >
      inline void DofMapperBase< GridPart, LocalDofMapping >
        ::onSubEntity( const ElementType &element, int i, int c, std::vector< bool > &indices ) const
      {
        const SubEntityFilter filter( RefElementsType::general( element.type() ), i, c );
        indices.resize( numDofs( element ) );
        code( element )( [ this, &indices, &filter ] ( unsigned int gtIndex, unsigned int subEntity, auto begin, auto end ) {
            const bool active = filter( subEntity, subEntityInfo_[ gtIndex ].codim );
            while( begin != end )
              indices[ *(begin++) ] = active;
          } );
      }

      template< class GridPart, class LocalDofMapping >
      template< class Entity >
      inline unsigned int
      DofMapperBase< GridPart, LocalDofMapping >
        ::numEntityDofs ( const Entity &entity ) const
      {
        return subEntityInfo( entity ).numDofs;
      }


      template< class GridPart, class LocalDofMapping >
      inline void DofMapperBase< GridPart, LocalDofMapping >::requestCodimensions ()
      {
        // this is only possible for index sets derived from Dune::Fem::IndexSet
        if constexpr ( Capabilities::isDuneFemIndexSet< IndexSetType >::v )
        {
          // collect all available codimensions
          std::vector< int > codimensions;
          codimensions.reserve( dimension+1 );

          for( typename BlockMapType::const_iterator it = blockMap_.begin(); it != blockMap_.end(); ++it )
          {
            SubEntityInfo &info = subEntityInfo_[ GlobalGeometryTypeIndex::index( *it ) ];
            codimensions.push_back( info.codim  );
          }

          // submit request for codimension to indexSet
          indexSet().requestCodimensions( codimensions );
        }
      }

      template< class GridPart, class LocalDofMapping >
      inline void DofMapperBase< GridPart, LocalDofMapping >::update ()
      {
        size_ = 0;
        for( const auto& geomType : blockMap_ )
        {
          SubEntityInfo &info = subEntityInfo_[ GlobalGeometryTypeIndex::index( geomType ) ];
          info.oldOffset = info.offset;
          info.offset = size_;
          size_ += SizeType( info.numDofs ) * SizeType( indexSet().size( geomType ) );
        }
      }


      template< class GridPart, class LocalDofMapping >
      inline const DofMapperCode &DofMapperBase< GridPart, LocalDofMapping >
        ::code ( const GeometryType &gt ) const
      {
        return code_[ LocalGeometryTypeIndex::index( gt ) ];
      }


      template< class GridPart, class LocalDofMapping >
      template< class Entity >
      inline const typename DofMapperBase< GridPart, LocalDofMapping >::SubEntityInfo &
      DofMapperBase< GridPart, LocalDofMapping >::subEntityInfo ( const Entity &entity ) const
      {
        return subEntityInfo_[ GlobalGeometryTypeIndex::index( entity.type() ) ];
      }



      // DofMapper
      // ---------

      template< class GridPart, class LocalDofMapping >
      class DofMapper
        : public DofMapperBase< GridPart, LocalDofMapping >
      {
        typedef DofMapper< GridPart, LocalDofMapping > ThisType;
        typedef DofMapperBase< GridPart, LocalDofMapping > BaseType;

      protected:
        typedef typename BaseType::SubEntityInfo SubEntityInfo;
        typedef typename BaseType::GridPartType::GridType GridType;
        typedef DofManager< GridType > DofManagerType;

      public:
        typedef typename BaseType::GridPartType GridPartType;
        typedef typename BaseType::LocalDofMappingType LocalDofMappingType;
        typedef typename BaseType::SizeType SizeType;

        template< class CodeFactory >
        DofMapper ( const GridPartType &gridPart, LocalDofMappingType localDofMapping, const CodeFactory &codeFactory )
          : BaseType( gridPart, std::move( localDofMapping ), codeFactory ),
            dofManager_( DofManagerType::instance( gridPart.grid() ) )
        {
          dofManager_.addIndexSet( *this );
        }

        DofMapper ( const ThisType & ) = delete;

        ~DofMapper ()
        {
          dofManager_.removeIndexSet( *this );
        }

        ThisType &operator= ( const ThisType & ) = delete;

      protected:
        DofManagerType& dofManager_;
      };


      // AdaptiveDofMapper
      // -----------------

      template< class GridPart, class LocalDofMapping >
      class AdaptiveDofMapper
        : public DofMapperBase< GridPart, LocalDofMapping >
      {
        typedef AdaptiveDofMapper< GridPart, LocalDofMapping > ThisType;
        typedef DofMapperBase< GridPart, LocalDofMapping > BaseType;

      protected:
        typedef typename BaseType::SubEntityInfo SubEntityInfo;
        typedef typename BaseType::GridPartType::GridType GridType;
        typedef DofManager< GridType > DofManagerType;

      public:
        typedef typename BaseType::GridPartType GridPartType;
        typedef typename BaseType::LocalDofMappingType LocalDofMappingType;
        typedef typename BaseType::SizeType SizeType;

        template< class CodeFactory >
        AdaptiveDofMapper ( const GridPartType &gridPart, LocalDofMappingType localDofMapping, const CodeFactory &codeFactory )
          : BaseType( gridPart, std::move( localDofMapping ), codeFactory ),
            dofManager_( DofManagerType::instance( gridPart.grid() ) )
        {
          dofManager_.addIndexSet( *this );
        }

        AdaptiveDofMapper ( const ThisType & ) = delete;

        ~AdaptiveDofMapper ()
        {
          dofManager_.removeIndexSet( *this );
        }

        ThisType &operator= ( const ThisType & ) = delete;

        // Adaptive DoF mappers are always up to date, so this method does nothing.
        // update is done several times during insertEntity
        void update () {}

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
        void insertEntity ( const Entity &entity ) { BaseType::update(); }

        template< class Entity >
        void removeEntity ( const Entity &entity ) {}

        void resize () { BaseType::update(); }

        bool compress () { BaseType::update(); return true; }

        template <class StreamTraits>
        void write( OutStreamInterface< StreamTraits >& out ) const {}

        template <class StreamTraits>
        void read( InStreamInterface< StreamTraits >& in )
        {
          BaseType::update();
        }

        void backup () const {}
        void restore () {}

      protected:
        using BaseType::indexSet;

        using BaseType::blockMap_;
        using BaseType::subEntityInfo_;

        DofManagerType& dofManager_;
      };



      // Implementation of AdaptiveDofMapper
      // -----------------------------------

      template< class GridPart, class LocalDofMapping >
      inline typename AdaptiveDofMapper< GridPart, LocalDofMapping >::SizeType
      AdaptiveDofMapper< GridPart, LocalDofMapping >::offSet ( int blk ) const
      {
        assert( (blk >= 0) && (blk < numBlocks()) );
        const unsigned int gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
        return subEntityInfo_[ gtIdx ].offset;
      }


      template< class GridPart, class LocalDofMapping >
      inline typename AdaptiveDofMapper< GridPart, LocalDofMapping >::SizeType
      AdaptiveDofMapper< GridPart, LocalDofMapping >::oldOffSet ( int blk ) const
      {
        assert( (blk >= 0) && (blk < numBlocks()) );
        const unsigned int gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
        return subEntityInfo_[ gtIdx ].oldOffset;
      }


      template< class GridPart, class LocalDofMapping >
      inline typename AdaptiveDofMapper< GridPart, LocalDofMapping >::SizeType
      AdaptiveDofMapper< GridPart, LocalDofMapping >::numberOfHoles ( int blk ) const
      {
        assert( (blk >= 0) && (blk < numBlocks()) );
        const unsigned int gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
        const SubEntityInfo &info = subEntityInfo_[ gtIdx ];
        return SizeType( info.numDofs ) * SizeType( indexSet().numberOfHoles( blockMap_[ blk ] ) );
      }


      template< class GridPart, class LocalDofMapping >
      inline typename AdaptiveDofMapper< GridPart, LocalDofMapping >::SizeType
      AdaptiveDofMapper< GridPart, LocalDofMapping >::oldIndex ( SizeType hole, int blk ) const
      {
        assert( (hole >= 0) && (hole < numberOfHoles( blk )) );
        const unsigned int gtIdx = GlobalGeometryTypeIndex::index( blockMap_[ blk ] );
        const SubEntityInfo &info = subEntityInfo_[ gtIdx ];
        const unsigned int numDofs = info.numDofs;
        const SizeType index = indexSet().oldIndex( hole / numDofs, blockMap_[ blk ] );
        return info.offset + numDofs * index + (hole % numDofs);
      }


      template< class GridPart, class LocalDofMapping >
      inline typename AdaptiveDofMapper< GridPart, LocalDofMapping >::SizeType
      AdaptiveDofMapper< GridPart, LocalDofMapping >::newIndex ( SizeType hole, int blk ) const
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

      template< class GridPart, class LocalDofMapping, bool adaptive = Capabilities::isAdaptiveIndexSet< typename GridPart::IndexSetType >::v >
      struct Implementation
      {
        typedef typename std::conditional< adaptive, AdaptiveDofMapper< GridPart, LocalDofMapping >, DofMapper< GridPart, LocalDofMapping > >::type Type;
      };

    } // namespace __IndexSetDofMapper



    // IndexSetDofMapper
    // -----------------

    template< class GridPart, class LocalDofMapping = DefaultLocalDofMapping< GridPart > >
    class IndexSetDofMapper
      : public __IndexSetDofMapper::template Implementation< GridPart, LocalDofMapping >::Type
    {
      typedef typename __IndexSetDofMapper::template Implementation< GridPart, LocalDofMapping >::Type BaseType;

    public:
      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::LocalDofMappingType LocalDofMappingType;

      template< class CodeFactory >
      IndexSetDofMapper ( const GridPartType &gridPart, LocalDofMappingType localDofMapping, const CodeFactory &codeFactory )
        : BaseType( gridPart, std::move( localDofMapping ), codeFactory )
      {}

      template< class CodeFactory >
      IndexSetDofMapper ( const GridPartType &gridPart, const CodeFactory &codeFactory )
        : BaseType( gridPart, LocalDofMappingType( gridPart ), codeFactory )
      {}
    };


    // Capabilities
    // ------------

    namespace Capabilities
    {
      // isAdaptiveDofMapper
      // -------------------

      template< class GridPart, class LocalDofMapping >
      struct isAdaptiveDofMapper< IndexSetDofMapper< GridPart, LocalDofMapping > >
      {
        static const bool v = Capabilities::isAdaptiveIndexSet< typename GridPart::IndexSetType >::v;
      };


      // isConsecutiveIndexSet
      // ---------------------

      template< class GridPart, class LocalDofMapping >
      struct isConsecutiveIndexSet< __IndexSetDofMapper::AdaptiveDofMapper< GridPart, LocalDofMapping > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_DOFMAPPER_INDEXSETDOFMAPPER_HH
