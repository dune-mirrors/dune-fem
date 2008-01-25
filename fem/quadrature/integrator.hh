#ifndef DUNE_FEM_INTEGRATOR_HH
#define DUNE_FEM_INTEGRATOR_HH

namespace Dune
{

  template< class Quadrature >
  class Integrator
  {
  public:
    typedef Quadrature QuadratureType;

    typedef typename QuadratureType :: EntityType EntityType;

  protected:
    typedef typename EntityType :: Geometry GeometryType;

  protected:
    const EntityType &entity_;
    const QuadratureType quadrature_;

  public:
    inline Integrator ( const EntityType &entity,
                        unsigned int order )
    : entity_( entity ),
      quadrature_( entity, order )
    {}

    template< class Function >
    inline void integrate ( const Function &function,
                            typename Function :: RangeType &ret ) const
    {
      typedef typename Function :: RangeType RangeType;
      typedef typename Function :: RangeFieldType RangeFieldType;

      const GeometryType &geometry = entity_.geometry();

      ret = 0;
      const unsigned int numQuadraturePoints = quadrature_.nop();
      for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
      {
        // evaluate function in quadrature point
        RangeType phi;
        function.evaluate( quadrature_[ pt ], phi );
       
        // calculate the weight of the quadrature point
        const RangeFieldType weight
          = geometry.integrationElement( quadrature_.point( pt ) )
            * quadrature_.weight( pt );

        ret.axpy( weight, phi );
      }
    }
  };
  
}

#endif
