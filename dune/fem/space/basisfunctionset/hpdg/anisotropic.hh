#ifndef DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_ANISOTROPIC_HH
#define DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_ANISOTROPIC_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <array>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/tensorproduct.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

#include "basisfunctionsets.hh"
#include "legendre.hh"

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // Internal forward declaration
      // ----------------------------

      template< class FunctionSpace, class GridPart, int maxOrder, class Storage >
      class AnisotropicBasisFunctionSets;



#ifndef DOXYGEN

      // LegendreShapeFunctionSetTuple
      // -----------------------------

      template< class FunctionSpace, int order, class Storage >
      class LegendreShapeFunctionSetTuple
      {
        // false == no hierarchical ordering
        using FactoryType = LegendreShapeFunctionSets< typename Dune::Fem::ToNewDimDomainFunctionSpace< FunctionSpace, 1 >::Type, order, false, Storage >;
        using ElementType = Dune::Fem::ShapeFunctionSetProxy< typename FactoryType::ShapeFunctionSetType >;

        template< int i, class MultiIndex >
        static ElementType get ( const MultiIndex &multiIndex )
        {
          return &FactoryType::get( multiIndex[ i ] );
        }

        template< class MultiIndex, int... i >
        static auto get ( const MultiIndex &multiIndex, std::integer_sequence< int,  i... > )
          -> decltype( std::make_tuple( get< i, MultiIndex >( multiIndex )... ) )
        {
          return std::make_tuple( get< i, MultiIndex >( multiIndex )... );
        }

      public:
        using Type = decltype( get( std::declval< Dune::FieldVector< int, FunctionSpace::dimDomain > >(), std::make_integer_sequence< int, FunctionSpace::dimDomain >() ) );

        template< class MultiIndex >
        static Type get ( const MultiIndex &multiIndex )
        {
          return get( multiIndex, std::make_integer_sequence< int, FunctionSpace::dimDomain >() );
        }
      };



      // AnisotropicShapeFunctionSet
      // ---------------------------

      template< class FunctionSpace, int order, class Storage >
      struct AnisotropicShapeFunctionSet
      : public Dune::Fem::TensorProductShapeFunctionSet< FunctionSpace, typename LegendreShapeFunctionSetTuple< FunctionSpace, order, Storage >::Type >
      {
        using BaseType = Dune::Fem::TensorProductShapeFunctionSet< FunctionSpace, typename LegendreShapeFunctionSetTuple< FunctionSpace, order, Storage >::Type >;

      public:
        AnisotropicShapeFunctionSet ()
          : AnisotropicShapeFunctionSet( multiIndex() )
        {}

        template< class MultiIndex >
        explicit AnisotropicShapeFunctionSet ( const MultiIndex &multiIndex )
          : BaseType( LegendreShapeFunctionSetTuple< FunctionSpace, order, Storage >::get( multiIndex ) )
        {}

      private:
        static Dune::FieldVector< int, FunctionSpace::dimDomain > multiIndex ()
        {
          return Dune::FieldVector< int, FunctionSpace::dimDomain >( order );
        }
      };



      // AnisotropicBasisFunctionSetsTraits
      // ----------------------------------

      template< class FunctionSpace, class GridPart, int maxOrder, class Storage >
      class AnisotropicBasisFunctionSetsTraits
      {
      public:
        using ImplementationType = AnisotropicBasisFunctionSets< FunctionSpace, GridPart, maxOrder, Storage >;

        using GridPartType = GridPart;
        using EntityType = typename GridPartType::template Codim< 0 >::EntityType;
        using Types = std::array< GeometryType, 1 >;

        using KeyType = Dune::FieldVector< int, FunctionSpace::dimDomain >;

        using ScalarFunctionSpaceType = typename Dune::Fem::ToNewDimRangeFunctionSpace< FunctionSpace, 1 >::Type;
        using ScalarShapeFunctionSetType = AnisotropicShapeFunctionSet< ScalarFunctionSpaceType, maxOrder, Storage >;
        using ShapeFunctionSetType = Dune::Fem::VectorialShapeFunctionSet< ScalarShapeFunctionSetType, typename FunctionSpace::RangeType >;

        using BasisFunctionSetType = DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType >;

        static const int localBlockSize = FunctionSpace::dimRange;

        using DataType = int;
      };

#endif // ifndef DOXYGEN



      // AnisotropicBasisFunctionSets
      // ----------------------------

      /** \brief A family of anisotropic local product basis function sets
       *
       *  \note This BasisFunctionSets object can be used only with cubic grids.
       *
       *  \tparam FunctionSpace  a Dune::Fem::FunctionSpace
       *  \tparam GridPart  a Dune::Fem::GridPart
       *  \tparam maxOrder  maximum order
       *  \tparam Storage  enable/disable Storage of quadratures
       *
       *  \ingroup DiscreteFunctionSpace_Implementation_Anisotropic
       */
      template< class FunctionSpace, class GridPart, int maxOrder, class Storage >
      class AnisotropicBasisFunctionSets
        : public BasisFunctionSets< AnisotropicBasisFunctionSetsTraits< FunctionSpace, GridPart, maxOrder, Storage > >
      {
        using BaseType = BasisFunctionSets< AnisotropicBasisFunctionSetsTraits< FunctionSpace, GridPart, maxOrder, Storage > >;

      public:
        /** \copydoc Dune::Fem::BasisFunctionSets::GridPartType */
        using GridPartType = typename BaseType::GridPartType;
        /** \copydoc Dune::Fem::BasisFunctionSets::EntityType */
        using EntityType = typename BaseType::EntityType;

        /** \copydoc Dune::Fem::BasisFunctionSets::BasisFunctionSetType */
        using BasisFunctionSetType = typename BaseType::BasisFunctionSetType;

      private:
        using ScalarShapeFunctionSetType = typename BaseType::Traits::ScalarShapeFunctionSetType;
        using ShapeFunctionSetType = typename BaseType::Traits::ShapeFunctionSetType;

      public:
        /** \copydoc Dune::Fem::BasisFunctionSets::KeyType */
        using KeyType = typename BaseType::KeyType;
        /** \copydoc Dune::Fem::BasisFunctionSets::DataType */
        using DataType = typename BaseType::DataType;

        /** \} */

        /** \copydoc Dune::Fem::BasisFunctionSets::types */
        typename BaseType::Types types () const noexcept
        {
          return std::array< GeometryType, 1 >{{ Dune::GeometryTypes::cube( EntityType::mydimension ) }};
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::maxBlocks */
        static constexpr std::size_t maxBlocks () noexcept
        {
          return Dune::power( int(maxOrder+1), int(FunctionSpace::dimDomain) );
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::maxBlocks */
        static std::size_t maxBlocks ( Dune::GeometryType type ) noexcept
        {
          assert( contains( type ) );
          return maxBlocks();
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::blocks */
        static std::size_t blocks ( GeometryType type, KeyType key ) noexcept
        {
          assert( contains( type ) );
          std::size_t blocks = 1;
          auto function = [&blocks]( int order ) -> void
          {
            assert( 0 <= order && order <= maxOrder );
            blocks *= order+1;
          };
          std::for_each( key.begin(), key.end(), function );
          assert( blocks == scalarShapeFunctionSet( key ).size() );
          return blocks;
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::encode */
        static DataType encode ( const KeyType &key ) noexcept
        {
          DataType data = 0, factor = 1;
          for( int i = 0; i < FunctionSpace::dimDomain-1; ++i )
          {
            data += key[ i ]*factor;
            factor *= maxOrder+1;
          }
          data += key[ FunctionSpace::dimDomain-1 ]*factor;
          assert( decode( data ) == key );
          return data;
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::decode */
        static KeyType decode ( DataType data ) noexcept
        {
          KeyType key;
          for( int i = 0; i < FunctionSpace::dimDomain-1; ++i )
          {
            key[ i ] = data % (maxOrder+1);
            data /= (maxOrder+1);
          }
          return std::move( key );
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::orthogonal */
        static constexpr bool orthogonal () noexcept
        {
          using GridType = typename GridPartType::GridType;
          return Dune::Capabilities::isCartesian< GridType >::v;
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::order */
        static constexpr int order () noexcept { return maxOrder; }

        /** \copydoc Dune::Fem::BasisFunctionSets::order */
        static constexpr int order ( Dune::GeometryType type ) noexcept
        {
          assert( type == GeometryType( GeometryType::cube, EntityType::mydimension ) );
          return order();
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::order */
        static int order ( Dune::GeometryType type, KeyType key ) noexcept
        {
          assert( type == GeometryType( GeometryType::cube, EntityType::mydimension ) );
          return *std::max_element( key.begin(), key.end() );
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::basisFunctionSet */
        static BasisFunctionSetType basisFunctionSet ( const EntityType &entity, const KeyType &key ) noexcept
        {
          assert( entity.type() == Dune::GeometryTypes::cube( EntityType::mydimension ) );
          return BasisFunctionSetType( entity, shapeFunctionSet( key ) );
        }

      private:
        static bool contains ( Dune::GeometryType type ) noexcept
        {
          return (type.isCube() && type.dim() == EntityType::mydimension);
        }

        static ScalarShapeFunctionSetType scalarShapeFunctionSet ( const KeyType &key ) noexcept
        {
          return ScalarShapeFunctionSetType( key );
        }

        static ShapeFunctionSetType shapeFunctionSet ( const KeyType &key ) noexcept
        {
          return ShapeFunctionSetType( scalarShapeFunctionSet( key ) );
        }
      };

    } // namespace hpDG

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_ANISOTROPIC_HH
