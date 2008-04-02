#ifndef DUNE_FEM_LINEARELLIPTICMODEL_HH
#define DUNE_FEM_LINEARELLIPTICMODEL_HH

#include <dune/common/bartonnackmanifcheck.hh>

#include <dune/fem/operator/model/diffusionmodel.hh>
#include <dune/fem/operator/model/boundarymodel.hh>

namespace Dune
{
  /** \ingroup EllipticOperator
   */

  struct DefaultLinearEllipticModelProperties
  {
    enum { hasDirichletValues = true };
    enum { hasNeumannValues = true };
    enum { hasRobinValues = true };
    enum { hasGeneralizedNeumannValues = true };
    enum { hasConvectiveFlux = true };
    enum { hasMass = true };
    enum { hasSource = true };
  };



  /** \class LinearEllipticModelInterface
   *  \brief interface for models of linear elliptic problems
   *
   *  This class models the data for equations of the following type:
   *  \f{displaymath}
   *    -\nabla \cdot (a(x) \nabla u) + \nabla \cdot (b( x ) u) + c( x ) u = f( x )
   *  \f}
   *  Here, the data functions are
   *  - \b a : the diffusive flux
   *  - \b b : the convective flux
   *  - \b c : the mass
   *  - \b f : the source
   *  .
   *  Additionally, boundary conditions are modelled by this class. For details,
   *  see BoundaryModelInterface.
   *
   *  Beside the properties of a BoundaryModelInterface, a LinearEllipticModel
   *  has the following properties:
   *  - \b hasConvectiveFlux : The convective flux is nontrivial
   *  - \b hasMass           : The mass in nontrivial
   *  - \b hasSource         : The source is nontrivial
   *  . 
   *  If one of these properties is \b false, the solver does not need to
   *  evaluate the corresponding terms. This may speed up the program.
   *
   *  \param  FunctionSpaceImp        function space to work in
   *  \param  LinearEllipticModelImp  actual implementation (Barton-Nackman)
   *  \param  PropertiesImp           properties of the implementation
   */
  template< class FunctionSpaceImp,
            class LinearEllipticModelImp,
            class PropertiesImp >
  class LinearEllipticModelInterface
  : public DiffusionModelInterface< FunctionSpaceImp, LinearEllipticModelImp >,
    public BoundaryModelInterface< FunctionSpaceImp, LinearEllipticModelImp, PropertiesImp >
  {
  public:
    //! type of the function space
    typedef FunctionSpaceImp FunctionSpaceType;

    //! type of the implementation (Barton-Nackman)
    typedef LinearEllipticModelImp LinearEllipticModelType;

    //! type properties the properties structure
    typedef PropertiesImp Properties;

  private:
    typedef LinearEllipticModelInterface
      < FunctionSpaceType, LinearEllipticModelType, Properties >
      ThisType;
    typedef DiffusionModelInterface< FunctionSpaceType, LinearEllipticModelType >
      BaseType;

  public:
    //! type of the interface
    typedef ThisType LinearEllipticModelInterfaceType;

    //! type of points within the domain
    typedef typename FunctionSpaceType :: DomainType DomainType;
    //! type of points within the range
    typedef typename FunctionSpaceType :: RangeType RangeType;
    //! type of the Jacobian (evaluated in some point)
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    //! field type of the domain
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    //! field type of the range
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  protected:
    using BaseType :: asImp;
    
  public:
    /** \brief evaluate the convective flux in a point
     *
     *  \param[in]   entity      entity to evaluate the flux on
     *  \param[in]   x           evaluaton point (in local coordinates)
     *  \param[in]   phi         value of the solution in the evaluation point
     *  \param[out]  flux        variable to receive the evaluated flux
     */
    template< class EntityType, class PointType >
    inline void convectiveFlux ( const EntityType &entity,
                                 const PointType &x,
                                 const RangeType &phi,
                                 JacobianRangeType &flux ) const
    {
      assert( Properties :: hasConvectiveFlux );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().convectiveFlux( entity, x, phi, flux ) );
    }
   
    /** \brief evaluate the mass in a point
     *
     *  \param[in]   entity      entity to evaluate the mass on
     *  \param[in]   x           evaluaton point (in local coordinates)
     *  \param[out]  ret         variable to receive the evaluated mass
     */
    template< class EntityType, class PointType >
    inline void mass ( const EntityType &entity,
                       const PointType &x,
                       RangeType &ret ) const
    {
      assert( Properties :: hasMass );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().mass( entity, x, ret ) );
    }

    /** \brief evaluate the source in a point
     *
     *  \param[in]   entity      entity to evaluate the source on
     *  \param[in]   x           evaluaton point (in local coordinates)
     *  \param[out]  ret         variable to receive the evaluated source
     */
    template< class EntityType, class PointType >
    inline void source ( const EntityType &entity,
                         const PointType &x,
                         RangeType &ret ) const
    {
      assert( Properties :: hasSource );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().source( entity, x, ret ) );
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

  public:
    using BaseType :: diffusiveFlux;
    using BaseType :: convectiveFlux;
    using BaseType :: mass;
    using BaseType :: source;

  protected:
    using BaseType :: asImp;

  private:
    const BoundaryModelDefaultType boundaryModelDefault_;

  public:
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    
  public:
    inline LinearEllipticModelDefault ()
    : boundaryModelDefault_()
    {
    }
    
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
    
    /** \copydoc Dune::DiffusionModelInterface::diffusiveFlux(const EntityType &entity,const PointType &x,const JacobianRangeType &gradient,JacobianRangeType &flux) const
     *
     *  The default implementation calls
     *  \code
     *  FieldVector< deriType, 0 > diffVar;
     *  diffusiveFlux( diffVar, entity, x, gradient, flux );
     *  \endcode
     */
    template< class EntityType, class PointType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const PointType &x,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux ) const
    {
      FieldVector< deriType, 0 > diffVar;
      asImp().diffusiveFlux( diffVar, entity, x, gradient, flux );
    }
   
    /** \copydoc Dune::LinearEllipticModelInterface::convectiveFlux(const EntityType &entity,const PointType &x,const RangeType &phi,JacobianRangeType &flux) const
     *
     *  The default implementation returns 0.
     */
    template< class EntityType, class PointType >
    inline void convectiveFlux ( const EntityType &entity,
                                 const PointType &x,
                                 const RangeType &phi,
                                 JacobianRangeType &flux ) const
    {
      assert( Properties :: hasConvectiveFlux );
      flux = 0;
    }

    /** \copydoc Dune::LinearEllipticModelInterface::mass(const EntityType &entity,const PointType &x,RangeType &ret) const
     *
     *  The default implementation returns 0.
     */
    template< class EntityType, class PointType >
    inline void mass ( const EntityType &entity,
                       const PointType &x,
                       RangeType &ret ) const
    {
      assert( Properties :: hasMass );
      ret = 0;
    }

    /** \copydoc Dune::LinearEllipticModelInterface::source(const EntityType &entity,const PointType &x,RangeType &ret) const
     *
     *  The default implementation returns 0.
     */
    template< class EntityType, class PointType >
    inline void source ( const EntityType &entity,
                         const PointType &x,
                         RangeType &ret ) const
    {
      assert( Properties :: hasSource );
      ret = 0;
    }
  };

}

#endif
