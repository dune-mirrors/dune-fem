#ifndef DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH
#define DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH

#include <dune/fem/function/common/discretefunction.hh>

namespace Dune
{

  namespace FemPy
  {

    // GridFunctionView
    // ----------------

    template< class GF, bool isDiscreteFunction = std::is_base_of< Fem::IsDiscreteFunction, GF >::value >
    struct GridFunctionView;

    template< class GF >
    struct GridFunctionView< GF, false >
    {
      typedef typename GF::EntityType Entity;
      typedef typename GF::RangeType Value;

      typedef typename Entity::Geometry::LocalCoordinate LocalCoordinate;

      GridFunctionView ( const GF &gf ) : localFunction_( gf ) {}

      Value operator() ( const LocalCoordinate &x ) const
      {
        Value value;
        localFunction_.evaluate( x, value );
        return value;
      }

      void bind ( const Entity &entity ) { localFunction_.init( entity ); }
      void unbind () {}

    private:
      typename GF::LocalFunctionType localFunction_;
    };

    template< class GF >
    struct GridFunctionView< GF, true >
    {
      typedef typename GF::EntityType Entity;
      typedef typename GF::RangeType Value;

      typedef typename Entity::Geometry::LocalCoordinate LocalCoordinate;

    private:
      typedef typename GF::DiscreteFunctionSpaceType DiscreteFunctionSpace;
      typedef typename DiscreteFunctionSpace::BasisFunctionSetType BasisFunctionSet;

    public:
      GridFunctionView ( const GF &gf )
        : gf_( gf )
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

    private:
      const DiscreteFunctionSpace &space () const { return gf_.space(); }

      const GF &gf_;
      BasisFunctionSet basisFunctionSet_;
      DynamicVector< typename GF::DofType > localDofVector_;
    };



    // localFunction
    // -------------

    template< class GF, std::enable_if_t< std::is_base_of< Fem::HasLocalFunction, GF >::value, int > = 0 >
    inline static GridFunctionView< GF > localFunction ( const GF &gf )
    {
      return GridFunctionView< GF >( gf );
    }

  } // namespace FemPy


  namespace Fem
  {

    using FemPy::localFunction;

  } // namespace Fem

  namespace Python
  {

    using FemPy::localFunction;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH
