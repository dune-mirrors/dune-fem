#ifndef DUNE_FEM_INTEGRATOR_HH
#define DUNE_FEM_INTEGRATOR_HH

namespace Dune
{

  /** \addtogroup Integrators
   * 
   *  Integrators are able to integrate a function.
   */



  /** \class   LocalIntegrator
   *  \ingroup Integrators
   *  \brief   integrator for a local function over a single entity
   *
   *  \param Quadrature  quadrature to use (either ElementQuadrature or
   *                     CachingQuadrature)
   */
  template< class Quadrature >
  class LocalIntegrator
  {
  public:
    //! type of quadrature to use
    typedef Quadrature QuadratureType;

    //! type of the entity
    typedef typename QuadratureType :: EntityType EntityType;

  protected:
    typedef typename EntityType :: Geometry GeometryType;

  protected:
    const EntityType &entity_;
    const QuadratureType quadrature_;

  public:
    /** \brief constructor
     * 
     *  \param[in]  entity  entity, for which the integrator is desired
     *  \param[in]  order   polynomial order for which the used quadrature
     *                      shall be exact
     */
    inline LocalIntegrator ( const EntityType &entity,
                             unsigned int order )
    : entity_( entity ),
      quadrature_( entity, order )
    {}

    /** \brief integrate a function
     *
     *  The function needs to have an evaluate method supporting
     *  \ref Dune::QuadraturePointWrapper "wrapped quadrature points".
     *  The declaration should look as follows:
     *  \code
     *  template< class Point >
     *  evaluate( Point &x, RangeType &ret );
     *  \endcode
     *
     *  \note The RangeType should be compatible with a Dune FieldVector.
     *
     *  \param[in]   function  function to integrate
     *  \param[out]  ret       value of the integral
     */
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
