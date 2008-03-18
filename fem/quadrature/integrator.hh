#ifndef DUNE_FEM_INTEGRATOR_HH
#define DUNE_FEM_INTEGRATOR_HH

//- Dune includes 
#include <dune/fem/quadrature/quadrature.hh>

namespace Dune
{

  /** \addtogroup Integrators
   * 
   *  Integrators are able to integrate a function.
   */



  /** \class   Integrator
   *  \ingroup Integrators
   *  \brief   integrator for arbitrary functions providing evaluate
   *
   *  \param Quadrature  quadrature to use (either ElementQuadrature or
   *                     CachingQuadrature)
   */
  template< class Quadrature >
  class Integrator
  {
  public:
    //! type of quadrature to use
    typedef Quadrature QuadratureType;

    //! type of the entity
    typedef typename QuadratureType :: EntityType EntityType;

  protected:
    typedef typename EntityType :: Geometry GeometryType;

  protected:
    const int order_;

  public:
    /** \brief constructor
     * 
     *  \param[in]  order   polynomial order for which the used quadrature
     *                      shall be exact
     */
    inline explicit Integrator ( unsigned int order )
    : order_( order )
    {}

    /** \brief add the integral over an entity to a variable
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
     *  \param[in]   entity    entity to integrate over
     *  \param[in]   function  function to integrate
     *  \param       ret       variable to which the value of the integral is
     *                         added
     */
    template< class Function >
    inline void integrateAdd ( const EntityType &entity,
                               const Function &function,
                               typename Function :: RangeType &ret ) const
    {
      typedef typename Function :: RangeType RangeType;
      typedef typename Function :: RangeFieldType RangeFieldType;

      const GeometryType &geometry = entity.geometry();
      const QuadratureType quadrature( entity, order_ );

      const unsigned int numQuadraturePoints = quadrature.nop();
      for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
      {
        // evaluate function in quadrature point
        RangeType phi;
        function.evaluate( quadrature[ pt ], phi );
       
        // calculate the weight of the quadrature point
        const RangeFieldType weight
          = geometry.integrationElement( quadrature.point( pt ) )
            * quadrature.weight( pt );

        ret.axpy( weight, phi );
      }
    }

    /** \brief integrate a function over an entity
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
     *  \param[in]   entity    entity to integrate over
     *  \param[in]   function  function to integrate
     *  \param[out]  ret       value of the integral
     */
    template< class Function >
    inline void integrate ( const EntityType &entity,
                            const Function &function,
                            typename Function :: RangeType &ret ) const
    {
      ret = 0;
      integrateAdd( entity, function, ret );
    }
  };
  
}

#endif
