#ifndef DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_LEGENDRE_HH
#define DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_LEGENDRE_HH

#include <cassert>
#include <cstddef>

#include <array>
#include <memory>
#include <type_traits>

#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

#include "basisfunctionsets.hh"

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // Internal forward declarations
      // -----------------------------

      template< class FunctionSpace, class GridPart, int maxOrder, bool hierarchicalOrdering, class Storage >
      class LegendreBasisFunctionSets;



#ifndef DOXYGEN

      // LegendreShapeFunctionSets
      // -------------------------

      template< class FunctionSpace, int order, bool hierarchicalOrdering, class Storage >
      class LegendreShapeFunctionSets
      {
        using ThisType = LegendreShapeFunctionSets< FunctionSpace, order, hierarchicalOrdering, Storage >;

      public:
        using ShapeFunctionSetType =
          Dune::Fem::SelectCachingShapeFunctionSet< LegendreShapeFunctionSet< FunctionSpace, hierarchicalOrdering >, Storage >;

      protected:
        LegendreShapeFunctionSets ()
        {
          for( int p = 0; p <= order; ++p )
            shapeFunctionSets_[ p ].reset( new ShapeFunctionSetType( type(), typename ShapeFunctionSetType::ImplementationType( p ) ) );
        }

        static const ThisType &instance ()
        {
          static ThisType instance;
          return instance;
        }

      public:
        static const ShapeFunctionSetType &get ( int p )
        {
          return *instance().shapeFunctionSets_[ p ];
        }

      private:
        static GeometryType type ()
        {
          return Dune::GeometryTypes::cube( FunctionSpace::dimDomain );
        }

        std::array< std::unique_ptr< ShapeFunctionSetType >, order+1 > shapeFunctionSets_;
      };



      // LegendreBasisFunctionSetsTraits
      // -------------------------------

      template< class FunctionSpace, class GridPart, int maxOrder, bool hierarchicalOrdering, class Storage >
      class LegendreBasisFunctionSetsTraits
      {
      public:
        using ImplementationType = LegendreBasisFunctionSets< FunctionSpace, GridPart, maxOrder, hierarchicalOrdering, Storage >;

        using GridPartType = GridPart;
        using Types = std::array< GeometryType, 1 >;

        using KeyType = int;
        using DataType = KeyType;

        using EntityType = typename GridPartType::template Codim< 0 >::EntityType;

        using ShapeFunctionSetsType = LegendreShapeFunctionSets< Dune::Fem::FunctionSpace< typename FunctionSpace::DomainFieldType, typename FunctionSpace::RangeFieldType, EntityType::mydimension, 1 >, maxOrder, hierarchicalOrdering, Storage >;
        using ShapeFunctionSetType = Dune::Fem::VectorialShapeFunctionSet< Dune::Fem::ShapeFunctionSetProxy< typename ShapeFunctionSetsType::ShapeFunctionSetType >, typename FunctionSpace::RangeType >;

        using BasisFunctionSetType = DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType >;

        static const int localBlockSize = BasisFunctionSetType::RangeType::dimension;
      };

#endif // #ifndef DOXYGEN



      // LegendreBasisFunctionSets
      // -------------------------

      /** \brief A family of local product basis function sets
       *
       *  \note This BasisFunctionSets object can be used only with cubic grids.
       *
       *  \tparam FunctionSpace  a Dune::Fem::FunctionSpace
       *  \tparam GridPart  a Dune::Fem::GridPart
       *  \tparam maxOrder  maximum order
       *  \tparam Storage  for certain caching features
       *
       *  \ingroup DiscreteFunctionSpace_Implementation_Legendre
       */
      template< class FunctionSpace, class GridPart, int maxOrder, bool hierarchicalOrdering, class Storage >
      class LegendreBasisFunctionSets
        : public BasisFunctionSets< LegendreBasisFunctionSetsTraits< FunctionSpace, GridPart, maxOrder, hierarchicalOrdering, Storage > >
      {
        using ThisType = LegendreBasisFunctionSets< FunctionSpace, GridPart, maxOrder, hierarchicalOrdering, Storage >;
        using BaseType = BasisFunctionSets< LegendreBasisFunctionSetsTraits< FunctionSpace, GridPart, maxOrder, hierarchicalOrdering, Storage > >;

      public:
        /** \copydoc Dune::Fem::BasisFunctionSets::GridPartType */
        using GridPartType = typename BaseType::GridPartType;
        /** \copydoc Dune::Fem::BasisFunctionSets::EntityType */
        using EntityType = typename BaseType::EntityType;

      private:
        using ShapeFunctionSetsType = typename BaseType::Traits::ShapeFunctionSetsType;
        using ShapeFunctionSetType = typename BaseType::Traits::ShapeFunctionSetType;

      public:
        /** \copydoc Dune::Fem::BasisFunctionSets::BasisFunctionSetType */
        using BasisFunctionSetType = typename BaseType::BasisFunctionSetType;

        /** \copydoc Dune::Fem::BasisFunctionSets::KeyType */
        using KeyType = typename BaseType::KeyType;
        /** \copydoc Dune::Fem::BasisFunctionSets::KeyType */
        using DataType = typename BaseType::DataType;

        /** \brief constructor */
        LegendreBasisFunctionSets () = default;

        /** \brief copy constructor */
        LegendreBasisFunctionSets ( const ThisType & ) = default;

        /** \brief move constructor */
        LegendreBasisFunctionSets ( ThisType && ) = default;

        /** \copydoc Dune::Fem::BasisFunctionSets::types */
        typename BaseType::Types types () const
        {
          return std::array< GeometryType, 1 >{{ Dune::GeometryTypes::cube( EntityType::mydimension ) }};
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::maxBlocks */
        std::size_t maxBlocks () const
        {
          return ShapeFunctionSetsType::get( maxOrder ).size();
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::maxBlocks */
        std::size_t maxBlocks ( GeometryType type ) const
        {
          return contains( type ) ? maxBlocks() : 0u;
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::blocks */
        std::size_t blocks ( GeometryType type, KeyType key ) const
        {
          return contains( type ) ? ShapeFunctionSetsType::get( key ).size() : 0u;
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::encode */
        static DataType encode ( const KeyType &key ) noexcept { return key; }

        /** \copydoc Dune::Fem::BasisFunctionSets::decode */
        static KeyType decode ( const DataType &data ) noexcept { return data; }

        /** \copydoc Dune::Fem::BasisFunctionSets::orthogonal */
        static constexpr bool orthogonal () noexcept
        {
          using GridType = typename GridPartType::GridType;
          return Dune::Capabilities::isCartesian< GridType >::v;
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::order */
        static constexpr int order () noexcept { return maxOrder; }

        /** \copydoc Dune::Fem::BasisFunctionSets::order */
        static constexpr int order ( GeometryType type ) noexcept { return order(); }

        /** \copydoc Dune::Fem::BasisFunctionSets::order */
        int order ( GeometryType type, KeyType key ) const
        {
          assert( contains( type ) );
          return static_cast< int >( key );
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::basisFunctionSet */
        BasisFunctionSetType basisFunctionSet ( const EntityType &entity, KeyType key ) const
        {
          assert( contains( entity.type() ) );
          return BasisFunctionSetType( entity, shapeFunctionSet( key ) );
        }

      private:
        static bool contains ( GeometryType type )
        {
          return (type.isCube() && type.dim() == EntityType::mydimension);
        }

        static ShapeFunctionSetType shapeFunctionSet ( KeyType key )
        {
          return ShapeFunctionSetType( &ShapeFunctionSetsType::get( static_cast< int >( key ) ) );
        }
      };

    } // namespace hpDG

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_LEGENDRE_HH
