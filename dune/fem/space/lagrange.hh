#ifndef HAVE_DUNE_FEM_SPACE_LAGRANGE
#define HAVE_DUNE_FEM_SPACE_LAGRANGE

#include <memory>

#include <dune/common/gmpfield.hh>

#include <dune/fem/space/lagrange/space.hh>

#if HAVE_DUNE_LOCALFUNCTIONS
#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#endif // #if HAVE_DUNE_LOCALFUNCTIONS

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/localfiniteelement/space.hh>
#include <dune/fem/space/localfiniteelement/dgspace.hh>
#include <dune/fem/space/localfiniteelement/quadratureinterpolation.hh>

namespace Dune
{

  namespace Fem
  {

    // LagrangeFiniteElementMap
    // ------------------------

#if HAVE_DUNE_LOCALFUNCTIONS
    template< class FunctionSpace, class GridPart, template< class, unsigned int > class PointSet = EquidistantPointSetDerived >
    class LagrangeFiniteElementMap
    {
      typedef LagrangeFiniteElementMap< FunctionSpace, GridPart, PointSet > ThisType;

    public:
      typedef GridPart GridPartType;

      typedef unsigned int KeyType;

      typedef typename FunctionSpace::DomainFieldType DomainFieldType;
      typedef typename FunctionSpace::RangeFieldType RangeFieldType;

      static const int dimLocal = GridPart::dimension;

      typedef LagrangeLocalFiniteElement< PointSet,dimLocal,double,double,
              // GMPField<64>, GMPField<256> > LocalFiniteElementType;
              double,
#if HAVE_GMP
              GMPField<256>
#else
              long double
#endif
                > LocalFiniteElementType;
      typedef typename LocalFiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType LocalCoefficientsType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolationType LocalInterpolationType;

      // -1 is default value if PointSet has no member pointSetId
      static const int pointSetId = detail::SelectPointSetId< PointSet<double, dimLocal> >::value;

      LagrangeFiniteElementMap ( const GridPart &gridPart, unsigned int order )
        : gridPart_( gridPart ), order_( order ), localFeVector_( size() )
      {}

      static std::size_t size () { return LocalGeometryTypeIndex::size(dimLocal); }

      int order () const { return order_; }

      template< class Entity >
      int order ( const Entity &entity ) const { return order(); }

      template< class Entity >
      std::tuple< std::size_t, const LocalBasisType &, const LocalInterpolationType & >
      operator() ( const Entity &e ) const
      {
        unsigned int index = localFiniteElement(e.type());
        const LocalFiniteElementType &lfe = *(localFeVector_[index]);
        return std::tuple< std::size_t, const LocalBasisType &, const LocalInterpolationType & >
          { index, lfe.localBasis(), lfe.localInterpolation() };
      }

      bool hasCoefficients ( const GeometryType &type ) const { return PointSet<double,0>::supports(type,order()); }

      const LocalCoefficientsType& localCoefficients ( const GeometryType &type ) const
      {
        unsigned int index = localFiniteElement(type);
        return localFeVector_[index]->localCoefficients();
      }

      const GridPartType &gridPart () const { return gridPart_; }

    private:
      std::size_t localFiniteElement ( const GeometryType &type ) const
      {
        std::size_t index = LocalGeometryTypeIndex::index(type);
        if ( !localFeVector_[ index ] )
          localFeVector_[ index ].reset( new LocalFiniteElementType( type, order_ ) );
        return index;
      }

      const GridPartType &gridPart_;
      unsigned int order_;
      mutable std::vector< std::unique_ptr< LocalFiniteElementType > > localFeVector_;
    };

    template< class FunctionSpace, class GridPart, unsigned int order,
      template< class, unsigned int > class PointSet = EquidistantPointSetDerived >
    struct FixedOrderLagrangeFiniteElementMap
    : public LagrangeFiniteElementMap<FunctionSpace,GridPart,PointSet>
    {
      typedef std::tuple<> KeyType; // no key
      static const unsigned int polynomialOrder = order;
      typedef typename GeometryWrapper<
          Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::topologyId,
                                                                  GridPart::dimension
        >::ImplType ImplType;
      typedef GenericLagrangeBaseFunction<
          typename FunctionSpace::ScalarFunctionSpaceType, ImplType, order
        > GenericBaseFunctionType;
      static const unsigned int numBasisFunctions = GenericBaseFunctionType::numBaseFunctions;
      FixedOrderLagrangeFiniteElementMap ( const GridPart &gridPart, const KeyType & )
        : LagrangeFiniteElementMap<FunctionSpace,GridPart,PointSet>(gridPart,order) {
        }
    };

    // LagrangeSpace
    // -------------

    template< class FunctionSpace, class GridPart,
              template< class, unsigned int > class PointSet = EquidistantPointSetDerived,
              class Storage = CachingStorage >
    using LagrangeSpace = LocalFiniteElementSpace<
            LagrangeFiniteElementMap< FunctionSpace, GridPart, PointSet >,
            FunctionSpace, Storage >;
    template< class FunctionSpace, class GridPart, unsigned int order,
              template< class, unsigned int > class PointSet = EquidistantPointSetDerived,
              class Storage = CachingStorage >
    using FixedOrderLagrangeSpace = LocalFiniteElementSpace<
            FixedOrderLagrangeFiniteElementMap< FunctionSpace, GridPart, order, PointSet >,
            FunctionSpace, Storage >;
    template< class FunctionSpace, class GridPart,
              template< class, unsigned int > class PointSet = EquidistantPointSetDerived,
              class Storage = CachingStorage >
    using DGLagrangeSpace = DiscontinuousLocalFiniteElementSpace<
            LagrangeFiniteElementMap< FunctionSpace, GridPart, PointSet >,
            FunctionSpace, Storage >;
    template< class FunctionSpace, class GridPart, unsigned int order,
              template< class, unsigned int > class PointSet = EquidistantPointSetDerived,
              class Storage = CachingStorage >
    using FixedOrderDGLagrangeSpace = DiscontinuousLocalFiniteElementSpace<
            FixedOrderLagrangeFiniteElementMap< FunctionSpace, GridPart, order, PointSet >,
            FunctionSpace, Storage >;
#else // #if HAVE_DUNE_LOCALFUNCTIONS
    // LagrangeSpace
    // -------------
    template< class FunctionSpace, class GridPart,
              class Storage = CachingStorage >
    using LagrangeSpace = DynamicLagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, Storage >;

#endif // #if HAVE_DUNE_LOCALFUNCTIONS

  } // namespace Fem

} // namespace Dune

#endif // #ifndef HAVE_DUNE_FEM_SPACE_LAGRANGE
