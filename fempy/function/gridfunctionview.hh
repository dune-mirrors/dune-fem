#ifndef DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH
#define DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH

#include <dune/fem/function/common/discretefunction.hh>

namespace Dune
{

  namespace FemPy
  {

    // GridFunctionView
    // ----------------

    template< class GF >
    struct GridFunctionView
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

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH
