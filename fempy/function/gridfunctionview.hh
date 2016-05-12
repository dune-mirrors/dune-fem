#ifndef DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH
#define DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH

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
      typename GF::LocalFunction localFunction_;
    };

  } // FemPy

} // Dune

#endif // #ifndef DUNE_FEMPY_FUNCTION_GRIDFUNCTIONVIEW_HH
