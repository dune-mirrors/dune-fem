#ifndef DUNE_FEM_SPACE_MAPPER_CODIMENSIONMAPPER_HH
#define DUNE_FEM_SPACE_MAPPER_CODIMENSIONMAPPER_HH

#include <cassert>

#include <algorithm>
#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/gridpart/common/indexset.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/mapper/dofmapper.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------

    template< class GridPart, int codim >
    class CodimensionMapper;



#ifndef DOXYGEN

    namespace __CodimensionMapper
    {

      // Traits
      // ------

      template< class GridPart, int codim >
      struct Traits
      {
        typedef CodimensionMapper< GridPart, codim > DofMapperType;

        static const int codimension = codim;

        typedef GridPart GridPartType;
        typedef typename GridPartType::IndexSetType IndexSetType;

        typedef typename GridPartType::template Codim< 0 >::EntityType ElementType;
        //typedef typename IndexSetType::IndexType SizeType;
        typedef std::size_t SizeType;
        typedef SizeType GlobalKeyType;
      };



      // DofMapper
      // ---------

      template< class T, template< class > class Base = Dune::Fem::DofMapper >
      class DofMapper
        : public Base< T >
      {
        typedef Base< T > BaseType;

      public:
        typedef typename BaseType::Traits Traits;

        static const int codimension = Traits::codimension;

        typedef typename Traits::GridPartType GridPartType;
        typedef typename Traits::IndexSetType IndexSetType;

        typedef typename BaseType::ElementType ElementType;
        typedef typename BaseType::SizeType SizeType;
        typedef typename Traits::GlobalKeyType GlobalKeyType;

        explicit DofMapper ( const GridPartType &gridPart )
          : DofMapper( gridPart.indexSet() )
        {}

        explicit DofMapper ( const IndexSetType &indexSet )
          : indexSet_( indexSet ),
            extension_( 0 ),
            maxNumDofs_( 0 )
        {
          AllGeomTypes< IndexSetType, typename GridPartType::GridType > types( indexSet );
          for( GeometryType type : types.geomTypes( 0 ) )
            maxNumDofs_ = std::max( maxNumDofs_, referenceElement( type ).size( codimension ) );

          if constexpr ( Capabilities::isDuneFemIndexSet< IndexSetType >:: v )
          {
            // submit codimension request to index set to enable support
            std::vector< int > codimensions( 1, int(Traits::codimension) );
            indexSet_.requestCodimensions( codimensions );
          }
        }

        /* \name DofMapper interface methods
         * \{
         */

        SizeType size () const
        {
          return indexSet().size( codimension ) + extension_;
        }

        static constexpr bool contains ( int codim ) noexcept
        {
          return codim == codimension;
        }

        static constexpr bool fixedDataSize ( int codim ) noexcept
        {
          return true;
        }

        template< class Function >
        void mapEach ( const ElementType &element, Function function ) const
        {
          if( codimension == 0 )
            function( 0, indexSet().index( element ) );
          else
          {
            const SizeType numDofs = this->numDofs( element );
            for( SizeType i = 0; i < numDofs; ++i )
              function( i, indexSet().subIndex( element, i, codimension ) );
          }
        }

        void map ( const ElementType &element, std::vector< GlobalKeyType > &indices ) const
        {
          indices.resize( numDofs( element ) );
          mapEach( element, [ &indices ] ( int local, GlobalKeyType global ) { indices[ local ] = global; } );
        }

        template< class Entity >
        void mapEntityDofs ( const Entity &entity, std::vector< GlobalKeyType > &indices ) const
        {
          indices.resize( numEntityDofs( entity ) );
          mapEachEntityDof( entity, AssignFunctor< std::vector< GlobalKeyType > >( indices ) );
        }

        template< class Entity, class Function >
        void mapEachEntityDof ( const Entity &entity, Function function ) const
        {
          assert( Entity::codimension == codimension );
          function( 0, indexSet().index( entity ) );
        }

        int maxNumDofs () const { return maxNumDofs_; }

        SizeType numDofs ( const ElementType &element ) const
        {
          return element.subEntities( codimension );
        }

        template< class Entity >
        static constexpr SizeType numEntityDofs ( const Entity &entity ) noexcept
        {
          return contains( Entity::codimension ) ? 1 : 0;
        }

        void update () {}

        void extendSize( const SizeType extension )
        {
          extension_ = extension;
        }

        void onSubEntity ( const ElementType &element, int i, int c, std::vector< bool > &indices ) const
        {
          indices.resize( numDofs(element) );
          std::fill(indices.begin(),indices.end(),false);
          indices[i] = true;
        }

        /* \} */

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
        const IndexSetType &indexSet () const { return indexSet_; }

      private:
        static auto referenceElement ( GeometryType type )
          -> decltype( ReferenceElements< typename GridPartType::ctype, GridPartType::dimension >::general( type ) )
        {
          return ReferenceElements< typename GridPartType::ctype, GridPartType::dimension >::general( type );
        }

        const IndexSetType &indexSet_;
        SizeType extension_;
        int maxNumDofs_;
      };



      // AdaptiveDofMapper
      // -----------------

      template< class Traits >
      class AdaptiveDofMapper
        : public DofMapper< Traits, Dune::Fem::AdaptiveDofMapper >
      {
        typedef DofMapper< Traits, Dune::Fem::AdaptiveDofMapper > BaseType;

      public:
        static const int codimension = BaseType::codimension;

        typedef typename BaseType::SizeType SizeType;
        typedef typename BaseType::GlobalKeyType GlobalKeyType;

      protected:
        using BaseType::indexSet;

      public:
        explicit AdaptiveDofMapper ( const typename BaseType::GridPartType &gridPart )
          : BaseType( gridPart )
        {}

        explicit AdaptiveDofMapper ( const typename BaseType::IndexSetType &indexSet )
          : BaseType( indexSet )
        {}

        static constexpr bool consecutive () noexcept { return true; }

        static constexpr SizeType numBlocks () noexcept
        {
          return 1;
        }

        SizeType numberOfHoles ( int ) const
        {
          return indexSet().numberOfHoles( codimension );
        }

        GlobalKeyType oldIndex ( int hole, int ) const
        {
          return indexSet().oldIndex( hole, codimension );
        }

        GlobalKeyType newIndex ( int hole, int ) const
        {
          return indexSet().newIndex( hole, codimension );
        }

        static constexpr SizeType oldOffSet ( int ) noexcept
        {
          return 0;
        }

        static constexpr SizeType offSet ( int ) noexcept
        {
          return 0;
        }
      };



      // Implementation
      // --------------

      template< class GridPart, int codim,
                bool adaptive = Capabilities::isAdaptiveIndexSet< typename GridPart::IndexSetType >::v >
      class Implementation
      {
        typedef __CodimensionMapper::Traits< GridPart, codim > Traits;

      public:
        typedef typename std::conditional< adaptive,
            AdaptiveDofMapper< Traits >,
            DofMapper< Traits >
          >::type Type;
      };

    } // namespace __CodimensionMapper

#endif // #ifndef DOXYGEN



    // CodimensionMapper
    // -----------------

    /** \brief mapper allocating one DoF per subentity of a given codimension
     *
     *  \tparam  GridPart  grid part type
     *  \tparam  codim  codimension
     *
     *  \note This mapper is adaptve (cf. AdaptiveDofMapper) if and only if the
     *        grid part's index set is adaptive, i.e. if
     *        Capabilities::isAdaptiveIndexSet< GridPart::IndexSetType >::v is \b true
     */
    template< class GridPart, int codim >
    class CodimensionMapper
      : public __CodimensionMapper::template Implementation< GridPart, codim >::Type
    {
      typedef typename __CodimensionMapper::template Implementation< GridPart, codim >::Type BaseType;

    public:
      explicit CodimensionMapper ( const typename BaseType::GridPartType &gridPart )
        : BaseType( gridPart )
      {}

      explicit CodimensionMapper ( const typename BaseType::IndexSetType *indexSet )
        : BaseType( *indexSet )
      {}

      explicit CodimensionMapper ( const typename BaseType::IndexSetType &indexSet )
        : BaseType( indexSet )
      {}
    };


    // Capabilities
    // ------------

    namespace Capabilities
    {
      // isAdaptiveDofMapper
      // -------------------

      template< class GridPart, int codim >
      struct isAdaptiveDofMapper< CodimensionMapper< GridPart, codim > >
      {
        static const bool v = Capabilities::isAdaptiveIndexSet< typename GridPart::IndexSetType >::v;
      };


      // isConsecutiveIndexSet
      // ---------------------

      template< class GridPart, int codim >
      struct isConsecutiveIndexSet< __CodimensionMapper::AdaptiveDofMapper< __CodimensionMapper::Traits< GridPart, codim > > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_MAPPER_CODIMENSIONMAPPER_HH
