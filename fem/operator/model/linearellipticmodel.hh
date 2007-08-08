#ifndef DUNE_FEM_LINEARELLIPTICMODEL_HH
#define DUNE_FEM_LINEARELLIPTICMODEL_HH

#include <dune/common/bartonnackmanifcheck.hh>

#include <dune/fem/operator/model/diffusionmodel.hh>
#include <dune/fem/operator/model/boundarymodel.hh>

namespace Dune
{
  /*! @ingroup EllipticOperator
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
  : public DiffusionModelInterface< FunctionSpaceImp, LinearEllipticModelImp >,
    public BoundaryModelInterface< FunctionSpaceImp, LinearEllipticModelImp, PropertiesImp >
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

    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
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
  : public LinearEllipticModelInterface< FunctionSpaceImp, LinearEllipticModelImp, PropertiesImp >
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

    typedef BoundaryModelDefault< FunctionSpaceType, LinearEllipticModelImp, Properties >
      BoundaryModelDefaultType;

  private:
    const BoundaryModelDefaultType boundaryModelDefault_;

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
      boundaryModelDefault_.dirichletValues( intersection, quadrature, point, phi );
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline void neumannValues ( const IntersectionIteratorType &intersection,
                                const QuadratureType &quadrature,
                                int point,
                                RangeType &phi ) const
    {
      boundaryModelDefault_.neumannValues( intersection, quadrature, point, phi );
    }
    
    template< class IntersectionIteratorType, class QuadratureType >
    inline void robinValues ( const IntersectionIteratorType &intersection,
                              const QuadratureType &quadrature,
                              int point,
                              RangeType &phi ) const
    {
      boundaryModelDefault_.robinValues( intersection, quadrature, point, phi );
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline RangeFieldType robinAlpha ( IntersectionIteratorType &intersection,
                                       QuadratureType &quadrature,
                                       int point ) const
    {
      return boundaryModelDefault_.robinAlpha( intersection, quadrature, point );
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
