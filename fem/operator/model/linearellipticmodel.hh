#ifndef DUNE_FEM_LINEARELLIPTICMODEL_HH
#define DUNE_FEM_LINEARELLIPTICMODEL_HH

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  class DefaultLinearEllipticModelProperties
  {
  public:
    enum { hasDirichletValues = true };
    enum { hasNeumannValues = true };
    enum { hasRobinValues = true };
    enum { hasGeneralizedNeumannValues = true };
    enum { hasConvectiveFlux = true };
    enum { hasMass = true };
    enum { hasSource = true };
  };



  template< class FunctionSpaceImp,
            class LinearEllipticModelImp,
            class PropertiesImp >
  class LinearEllipticModelInterface
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef PropertiesImp Properties;

  private:
    typedef LinearEllipticModelInterface
      < FunctionSpaceType, LinearEllipticModelImp, Properties >
      ThisType;

  public:
    typedef ThisType LinearEllipticModelInterfaceType;

    enum BoundaryType { Dirichlet, Neumann, Robin, GeneralizedNeumann };

    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    /** \brief obtain type of boundary intersection
     *
     *  \param[in] intersection intersection whose boundary type shall be returned
     *
     *  \returns boundary type of the intersection
     */
    template< class IntersectionIteratorType >
    inline BoundaryType boundaryType ( const IntersectionIteratorType &intersection ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION
        ( asImp().boundaryType( intersection ) );
      return asImp().boundaryType( intersection );
    }

    /** \brief evaluate Dirichlet boundary value
     *
     *  Given an intersection and a quadrature-point pair, evaluate the
     *  Dirichlet boundary data in the specified point.
     *
     *  \note The entity the boundary belongs to can be obtained from the
     *        intersection by the intersecion.inside method.
     *
     *  \note The quadrature is not guaranteed to be a codim-1-quadrature
     *
     *  \param[in]  intersection intersection on which the quadrature point lies
     *
     *  \param{in]  quadrature   quadrature to use
     *
     *  \param[in]  point        index of the evaluation point within the
     *                           quadrature
     *
     *  \param[out] phi          Dirichlet value in the quadrature point
     */
    template< class IntersectionIteratorType, class QuadratureType >
    inline void dirichletValues ( const IntersectionIteratorType &intersection,
                                  const QuadratureType &quadrature,
                                  int point,
                                  RangeType &phi ) const
    {
      assert( Properties :: hasDirichletValues );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().dirichletValues( intersection, quadrature, point, phi ) );
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline void neumannValues ( const IntersectionIteratorType &intersection,
                                const QuadratureType &quadrature,
                                int point,
                                RangeType &phi ) const
    {
      assert( Properties :: hasNeumannValues );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().neumannValues( intersection, quadrature, point, phi ) );
    }
    
    template< class IntersectionIteratorType, class QuadratureType >
    inline void generalizedNeumannValues
      ( const IntersectionIteratorType &intersection,
        const QuadratureType &quadrature,
        int point,
        RangeType &phi
      ) const
    {
      assert( Properties :: hasGeneralizedNeumannValues );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().generalizedNeumannValues( intersection, quadrature, point, phi ) );
    }


    template< class IntersectionIteratorType, class QuadratureType >
    inline RangeFieldType generalizedNeumannAlpha
      ( const IntersectionIteratorType &intersection,
        const QuadratureType &quadrature,
        int point
      ) const
    {
      assert( Properties :: hasGeneralizedNeumannValues );
      CHECK_INTERFACE_IMPLEMENTATION
        ( asImp().generalizedNeumannAlpha( intersection, quadrature, point ) );
      return asImp().generalizedNeumannAlpha( intersection, quadrature, point );
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline void robinValues ( const IntersectionIteratorType &intersection,
                              const QuadratureType &quadrature,
                              int point,
                              RangeType &phi ) const
    {
      assert( Properties :: hasRobinValues );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().robinValues( intersection, quadrature, point, phi ) );
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline RangeFieldType robinAlpha ( const IntersectionIteratorType &intersection,
                                       const QuadratureType &quadrature,
                                       int point ) const
    {
      assert( Properties :: hasRobinValues );
      CHECK_INTERFACE_IMPLEMENTATION
        ( asImp().robinAlpha( intersection, quadrature, point ) );
      return asImp().robinAlpha( intersection, quadrature, point );
    }

    template< class EntityType, class QuadratureType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const QuadratureType &quadrature,
                                int point,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().diffusiveFlux( entity, quadrature, point, gradient, flux ) );
    }
    
    template< class EntityType, class QuadratureType >
    inline void convectiveFlux ( const EntityType &entity,
                                 const QuadratureType &quadrature,
                                 int point,
                                 const RangeType &phi,
                                 DomainType &flux ) const
    {
      assert( Properties :: hasConvectiveFlux );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().convectiveFlux( entity, quadrature, point, phi, flux ) );
    }

    template< class EntityType, class QuadratureType >
    inline void mass ( const EntityType &entity,
                       const QuadratureType &quadrature,
                       int point,
                       RangeType &ret ) const
    {
      assert( Properties :: hasMass );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().mass( entity, quadrature, point, ret ) );
    }

    template< class EntityType, class QuadratureType >
    inline void source ( const EntityType &entity,
                         const QuadratureType &quadrature,
                         int point,
                         RangeType &ret ) const
    {
      assert( Properties :: hasSource );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().source( entity, quadrature, point, ret ) );
    }

  protected:
    inline const LinearEllipticModelImp &asImp () const
    {
      return static_cast< const LinearEllipticModelImp& >( *this );
    }

    inline LinearEllipticModelImp &asImp ()
    {
      return static_cast< LinearEllipticModelImp& >( *this );
    }
  };



  template< class FunctionSpaceImp,
            class LinearEllipticModelImp,
            class PropertiesImp
              = DefaultLinearEllipticModelProperties >
  class LinearEllipticModelDefault
  : public LinearEllipticModelInterface
    < FunctionSpaceImp, LinearEllipticModelImp, PropertiesImp >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef PropertiesImp Properties;

  private:
    typedef LinearEllipticModelDefault
      < FunctionSpaceType, LinearEllipticModelImp, Properties >
      ThisType;
    typedef LinearEllipticModelInterface
      < FunctionSpaceType, LinearEllipticModelImp, Properties >
      BaseType;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    
  public:
    template< class IntersectionIteratorType, class QuadratureType >
    inline void dirichletValues ( const IntersectionIteratorType &intersection,
                                  const QuadratureType &quadrature,
                                  int point,
                                  RangeType &phi ) const
    {
      assert( Properties :: hasDirichletValues );
      phi = 0;
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline void neumannValues ( const IntersectionIteratorType &intersection,
                                const QuadratureType &quadrature,
                                int point,
                                RangeType &phi ) const
    {
      assert( Properties :: hasNeumannValues );
      phi = 0;
    }
    
    template< class IntersectionIteratorType, class QuadratureType >
    inline void generalizedNeumannValues
      ( const IntersectionIteratorType &intersection,
        const QuadratureType &quadrature,
        int point,
        RangeType &phi
      ) const
    {
      assert( Properties :: hasGeneralizedNeumannValues );
      phi = 0;
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline RangeFieldType generalizedNeumannAlpha
      ( const IntersectionIteratorType &intersection,
        const QuadratureType &quadrature,
        int point
      ) const
    {
      assert( Properties :: hasGeneralizedNeumannValues );
      return 0;
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline void robinValues ( const IntersectionIteratorType &intersection,
                              const QuadratureType &quadrature,
                              int point,
                              RangeType &phi ) const
    {
      assert( Properties :: hasRobinValues );
      phi = 0;
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline RangeFieldType robinAlpha ( IntersectionIteratorType &intersection,
                                       QuadratureType &quadrature,
                                       int point ) const
    {
      assert( Properties :: hasRobinValues );
      return 0;
    }

    template< class EntityType, class QuadratureType >
    inline void convectiveFlux ( const EntityType &entity,
                                 const QuadratureType &quadrature,
                                 int point,
                                 const RangeType &phi,
                                 DomainType &flux ) const
    {
      assert( Properties :: hasConvectiveFlux );
      flux = 0;
    }

    template< class EntityType, class QuadratureType >
    inline void mass ( EntityType &entity,
                       QuadratureType &quadrature,
                       int point,
                       RangeType &ret ) const
    {
      assert( Properties :: hasMass );
      ret = 0;
    }

    template< class EntityType, class QuadratureType >
    inline void source ( const EntityType &entity,
                         const QuadratureType &quadrature,
                         int point,
                         RangeType &ret ) const
    {
      assert( Properties :: hasSource );
      source = 0;
    }
  };

}

#endif
