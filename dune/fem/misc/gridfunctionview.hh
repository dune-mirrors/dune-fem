#ifndef DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH
#define DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/common/intersectionside.hh>

namespace Dune
{

  namespace Fem
  {

    // GridFunctionView
    // ----------------
    struct IsGridFunctionView {};

    template< class GF, bool isDiscreteFunction = std::is_base_of< Fem::IsDiscreteFunction, GF >::value >
    struct GridFunctionView;

    template< class GF >
    struct GridFunctionView< GF, false >
    : public BindableGridFunction<typename GF::GridPartType, typename GF::RangeType>,
             public IsGridFunctionView
    {
      using Base = BindableGridFunction<typename GF::GridPartType, typename GF::RangeType>;
      typedef typename GF::EntityType Entity;
      typedef typename GF::RangeType Value;

      typedef typename Entity::Geometry::LocalCoordinate LocalCoordinate;

      GridFunctionView ( const GF &gf )
      : Base(gf.gridPart())
      , localFunction_( gf ) {}

      Value operator() ( const LocalCoordinate &x ) const
      {
        Value value;
        localFunction_.evaluate( x, value );
        return value;
      }

      void bind ( const Entity &entity ) { localFunction_.bind( entity ); }
      void unbind () {}
      template <class IntersectionType>
      void bind(const IntersectionType &intersection, IntersectionSide side)
      { localFunction_.bind(intersection, side); }

    private:
      Dune::Fem::ConstLocalFunction<GF> localFunction_;
    };

    template< class GF >
    struct GridFunctionView< GF, true >
    : public BindableGridFunctionWithSpace<typename GF::GridPartType, typename GF::RangeType>,
             public IsGridFunctionView
    {
      using Base = BindableGridFunctionWithSpace<typename GF::GridPartType, typename GF::RangeType>;
      typedef typename GF::EntityType Entity;
      typedef typename GF::RangeType Value;

      typedef typename Entity::Geometry::LocalCoordinate LocalCoordinate;

    private:
      typedef typename GF::DiscreteFunctionSpaceType DiscreteFunctionSpace;
      typedef typename DiscreteFunctionSpace::BasisFunctionSetType BasisFunctionSet;

    public:
      GridFunctionView ( const GF &gf )
        : Base(gf.gridPart(), gf.name(), gf.order())
        , gf_( gf )
      {
        localDofVector_.reserve( DiscreteFunctionSpace::localBlockSize * space().blockMapper().maxNumDofs() );
      }

      Value operator() ( const LocalCoordinate &x ) const
      {
        Value value;
        basisFunctionSet_.evaluateAll( x, localDofVector_, value );
        return value;
      }

      void bind ( const Entity &entity )
      {
        basisFunctionSet_ = space().basisFunctionSet( entity );
        localDofVector_.resize( basisFunctionSet_.size() );
        gf_.getLocalDofs( entity, localDofVector_ );
      }

      void unbind ()
      {
        basisFunctionSet_ = BasisFunctionSet();
        localDofVector_.resize( 0u );
      }

      template <class IntersectionType>
      void bind(const IntersectionType &intersection, IntersectionSide side)
      { defaultIntersectionBind(gf_, intersection, side); }

    private:
      const DiscreteFunctionSpace &space () const { return gf_.space(); }

      const GF &gf_;
      BasisFunctionSet basisFunctionSet_;
      DynamicVector< typename GF::DofType > localDofVector_;
    };



    // localFunction
    // -------------

    template< class GF,
      std::enable_if_t< !std::is_base_of< Fem::IsGridFunctionView, GF >::value &&
                         std::is_base_of< Fem::HasLocalFunction, GF >::value, int > = 0
      >
    inline static GridFunctionView< GF > localFunction ( const GF &gf )
    {
      return GridFunctionView< GF >( gf );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH
