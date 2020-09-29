#ifndef DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_ORTHOGONAL_HH
#define DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_ORTHOGONAL_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <array>
#include <memory>
#include <type_traits>
#include <vector>

#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

#include <dune/fem/space/shapefunctionset/orthonormal.hh>

#include "typeindexset.hh"
#include "typemap.hh"
#include "basisfunctionsets.hh"

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // Internal forward declarations
      // -----------------------------

      template< class FunctionSpace, class GridPart, int maxOrder, class Storage >
      class OrthogonalBasisFunctionSets;



#ifndef DOXYGEN

      // OrthonormalShapeFunctionSets
      // ----------------------------

      template< class FunctionSpace, int order, class Storage >
      class OrthonormalShapeFunctionSets
      {
        using ThisType = OrthonormalShapeFunctionSets< FunctionSpace, order, Storage >;

      public:
        /** \brief type of shape function set */
        using ShapeFunctionSetType =
            Dune::Fem::SelectCachingShapeFunctionSet< OrthonormalShapeFunctionSet< FunctionSpace >, Storage >;

      private:
        static const int dimension = ShapeFunctionSetType::DomainType::dimension;

      protected:
        OrthonormalShapeFunctionSets ()
        {
          using GeometryTypeIndexSet = LocalGeometryTypeIndexSet< dimension, true >;
          const std::size_t types = GeometryTypeIndexSet::size();
          for( std::size_t i = 0u; i < types; ++i )
          {
            GeometryType type = GeometryTypeIndexSet::type( i );
            assert( !type.isNone() );
            types_.push_back( type );
            for( int p = 0; p <= order; ++p )
              shapeFunctionSets_[ type ][ p ].reset( new ShapeFunctionSetType( type, typename ShapeFunctionSetType::ImplementationType( type, p ) ) );
          }
        }

        static const ThisType &instance ()
        {
          static ThisType instance;
          return instance;
        }

      public:
        /** \brief vector of supported geometries */
        static const std::vector< GeometryType > &types ()
        {
          return instance().types();
        }

        /** \brief size of shape function set */
        static std::size_t maxSize ()
        {
          return OrthonormalShapeFunctions< dimension >::size( order );
        }

        /** \brief return size of shape function set for given geometry and polynomial order */
        static std::size_t size ( GeometryType type, int p )
        {
          return OrthonormalShapeFunctions< dimension >::size( p );
        }

        /** \brief get shape function set for given geometry and polynomial order */
        static const ShapeFunctionSetType &get ( GeometryType type, int p )
        {
          return *instance().shapeFunctionSets_[ type ][ p ];
        }

      private:
        std::vector< GeometryType > types_;
        LocalGeometryTypeMap< std::array< std::unique_ptr< ShapeFunctionSetType >, order+1 >, dimension > shapeFunctionSets_;
      };



      // OrthogonalBasisFunctionSetsTraits
      // ---------------------------------

      template< class FunctionSpace, class GridPart, int maxOrder, class Storage >
      class OrthogonalBasisFunctionSetsTraits
      {
      public:
        using ImplementationType = OrthogonalBasisFunctionSets< FunctionSpace, GridPart, maxOrder, Storage >;

        using GridPartType = GridPart;
        using Types = const std::vector< GeometryType > &;

        using KeyType = int;
        using DataType = KeyType;

        using EntityType = typename GridPartType::template Codim< 0 >::EntityType;

        using ShapeFunctionSetsType = OrthonormalShapeFunctionSets< Dune::Fem::FunctionSpace< typename FunctionSpace::DomainFieldType, typename FunctionSpace::RangeFieldType, EntityType::mydimension, 1 >, maxOrder, Storage >;
        using ShapeFunctionSetType = Dune::Fem::VectorialShapeFunctionSet< Dune::Fem::ShapeFunctionSetProxy< typename ShapeFunctionSetsType::ShapeFunctionSetType >, typename FunctionSpace::RangeType >;

        using BasisFunctionSetType = DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType >;

        static const int localBlockSize = BasisFunctionSetType::RangeType::dimension;
      };

#endif // #ifndef DOXYGEN



      // OrthogonalBasisFunctionSets
      // ---------------------------

      /** \brief A family of orthogonal local basis function sets
       *
       *  \tparam FunctionSpace  a Dune::Fem::FunctionSpace
       *  \tparam GridPart  a Dune::Fem::GridPart
       *  \tparam maxOrder  maximum order
       *  \tparam Storage  enable/disable Storage of quadratures
       *
       *  \ingroup DiscreteFunctionSpace_Implementation_Orthogonal
       */
      template< class FunctionSpace, class GridPart, int maxOrder, class Storage >
      class OrthogonalBasisFunctionSets
        : public BasisFunctionSets< OrthogonalBasisFunctionSetsTraits< FunctionSpace, GridPart, maxOrder, Storage > >
      {
        using ThisType = OrthogonalBasisFunctionSets< FunctionSpace, GridPart, maxOrder, Storage >;
        using BaseType = BasisFunctionSets< OrthogonalBasisFunctionSetsTraits< FunctionSpace, GridPart, maxOrder, Storage > >;

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

        /** \brief type of function space */
        using FunctionSpaceType = FunctionSpace;

        /** \brief range type of basis functions */
        using DomainType = typename FunctionSpaceType :: DomainType;

        /** \brief range type of basis functions */
        using RangeType = typename FunctionSpaceType :: RangeType;

        /** \name Construction
         *  \{
         */

        /** \brief constructor */
        OrthogonalBasisFunctionSets () = default;

        /** \brief copy constructor */
        OrthogonalBasisFunctionSets ( const ThisType & ) = default;

        /** \brief move constructor */
        OrthogonalBasisFunctionSets ( ThisType && ) = default;

        /** \} */

        /** \name Assignment
         *  \{
         */

        /** \brief assignment operator */
        ThisType &operator= ( const ThisType & ) = default;

        /** \brief move assignment operator */
        ThisType &operator= ( ThisType && ) = default;

        /** \} */

        /** \copydoc Dune::Fem::BasisFunctionSets::types */
        static const std::vector< GeometryType > &types ()
        {
          return ShapeFunctionSetsType::types();
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::maxBlocks */
        static std::size_t maxBlocks ()
        {
          return ShapeFunctionSetsType::maxSize();
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::maxBlocks */
        static std::size_t maxBlocks ( GeometryType type )
        {
          return contains( type ) ? maxBlocks() : 0u;
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::blocks */
        static std::size_t blocks ( GeometryType type, KeyType key )
        {
          return contains( type ) ? ShapeFunctionSetsType::size( type, key ) : 0u;
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::encode */
        static DataType encode ( const KeyType &key ) noexcept { return key; }

        /** \copydoc Dune::Fem::BasisFunctionSets::decode */
        static KeyType decode ( const DataType &data ) noexcept { return data; }

        /** \copydoc Dune::Fem::BasisFunctionSets::orthogonal */
        static constexpr bool orthogonal () noexcept
        {
          //using GridType = typename GridPartType::GridType;
          //return Dune::Capabilities::isCartesian< GridType >::v;
          return true;
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::order */
        static constexpr int order () noexcept { return maxOrder; }

        /** \copydoc Dune::Fem::BasisFunctionSets::order */
        static int order ( GeometryType type ) noexcept { return order(); }

        /** \copydoc Dune::Fem::BasisFunctionSets::order */
        static int order ( GeometryType type, KeyType key )
        {
          assert( contains( type ) );
          return static_cast< int >( key );
        }

        /** \copydoc Dune::Fem::BasisFunctionSets::basisFunctionSet */
        static BasisFunctionSetType basisFunctionSet ( const EntityType &entity, KeyType key )
        {
          assert( contains( entity.type() ) );
          return BasisFunctionSetType( entity, shapeFunctionSet( entity.type(), key ) );
        }

      private:
        static bool contains ( GeometryType type )
        {
          return (type.dim() <= 3 && !type.isNone());
        }

        static ShapeFunctionSetType shapeFunctionSet ( GeometryType type, KeyType key )
        {
          return ShapeFunctionSetType( &ShapeFunctionSetsType::get( type, static_cast< int >( key ) ) );
        }
      };

    } // namespace hpDG

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_BASISFUNCTIONSETS_ORTHOGONAL_HH
