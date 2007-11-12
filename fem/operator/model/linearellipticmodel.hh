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

    //! type properties the properties structure
    typedef PropertiesImp Properties;

  private:
    typedef LinearEllipticModelInterface
      < FunctionSpaceType, LinearEllipticModelImp, Properties >
      ThisType;

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

  public:
#if 0
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
#endif

    /** \brief evaluate the convective flux in a point
     *
     *  \param[in]   entity      entity to evaluate the flux on
     *  \param[in]   x           evaluaton point (in local coordinates)
     *  \param[in]   phi         value of the solution in the evaluation point
     *  \param[out]  flux        variable to receive the evaluated flux
     */
    template< class EntityType >
    inline void convectiveFlux ( const EntityType &entity,
                                 const DomainType &x,
                                 const RangeType &phi,
                                 JacobianRangeType &flux ) const
    {
      assert( Properties :: hasConvectiveFlux );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().convectiveFlux( entity, x, phi, flux ) );
    }
    
    /** \brief evaluate the convective flux in a quadrature point
     *
     *  \param[in]   entity      entity to evaluate the flux on
     *  \param[in]   quadrature  quadrature to use
     *  \param[in]   quadPoint   number of the point within the quadrature
     *  \param[in]   phi         value of the solution in the evaluation point
     *  \param[out]  flux        variable to receive the evaluated flux
     */
    template< class EntityType, class QuadratureType >
    inline void convectiveFlux ( const EntityType &entity,
                                 const QuadratureType &quadrature,
                                 const int quadPoint,
                                 const RangeType &phi,
                                 JacobianRangeType &flux ) const
    {
      assert( Properties :: hasConvectiveFlux );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().convectiveFlux( entity, quadrature, quadPoint, phi, flux ) );
    }
    
    /** \brief evaluate the mass in a point
     *
     *  \param[in]   entity      entity to evaluate the mass on
     *  \param[in]   x           evaluaton point (in local coordinates)
     *  \param[out]  ret         variable to receive the evaluated mass
     */
    template< class EntityType >
    inline void mass ( const EntityType &entity,
                       const DomainType &x,
                       RangeType &ret ) const
    {
      assert( Properties :: hasMass );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().mass( entity, x, ret ) );
    }

    /** \brief evaluate the mass in a quadrature point
     *
     *  \param[in]   entity      entity to evaluate the mass on
     *  \param[in]   quadrature  quadrature to use
     *  \param[in]   quadPoint   number of the point within the quadrature
     *  \param[out]  ret         variable to receive the evaluated mass
     */
    template< class EntityType, class QuadratureType >
    inline void mass ( const EntityType &entity,
                       const QuadratureType &quadrature,
                       const int quadPoint,
                       RangeType &ret ) const
    {
      assert( Properties :: hasMass );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().mass( entity, quadrature, quadPoint, ret ) );
    }
    
    /** \brief evaluate the source in a point
     *
     *  \param[in]   entity      entity to evaluate the source on
     *  \param[in]   x           evaluaton point (in local coordinates)
     *  \param[out]  ret         variable to receive the evaluated source
     */
    template< class EntityType >
    inline void source ( const EntityType &entity,
                         const DomainType &x,
                         RangeType &ret ) const
    {
      assert( Properties :: hasSource );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().source( entity, x, ret ) );
    }

    /** \brief evaluate the source in a quadrature point
     *
     *  \param[in]   entity      entity to evaluate the source on
     *  \param[in]   quadrature  quadrature to use
     *  \param[in]   quadPoint   number of the point within the quadrature
     *  \param[out]  ret         variable to receive the evaluated source
     */
    template< class EntityType, class QuadratureType >
    inline void source ( const EntityType &entity,
                         const QuadratureType &quadrature,
                         const int quadPoint,
                         RangeType &ret ) const
    {
      assert( Properties :: hasSource );
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().source( entity, quadrature, quadPoint, ret ) );
    }

  protected:
    inline const LinearEllipticModelImp &asImp () const
    {
      return static_cast< const LinearEllipticModelImp & >( *this );
    }

    inline LinearEllipticModelImp &asImp ()
    {
      return static_cast< LinearEllipticModelImp & >( *this );
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


#if 0
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
#endif
    
    /** \copydoc Dune::DiffusionModelInterface::diffusiveFlux(const EntityType &entity,const QuadratureType &quadrature,const int quadPoint,const JacobianRangeType &gradient,JacobianRangeType &flux) const
     *
     *  The default implementation calls
     *  \code
     *  diffusiveFlux( entity, quadrature.point( quadpoint ), gradient, flux );
     *  \endcode
     */
    template< class EntityType, class QuadratureType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const QuadratureType &quadrature,
                                const int quadPoint,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux ) const
    {
      asImp().diffusiveFlux( entity, quadrature.point( quadPoint ), gradient, flux );
    }

    /** \copydoc Dune::LinearEllipticModelInterface::convectiveFlux(const EntityType &entity,const DomainType &x,const RangeType &phi,JacobianRangeType &flux) const
     *
     *  The default implementation returns 0.
     */
    template< class EntityType >
    inline void convectiveFlux ( const EntityType &entity,
                                 const DomainType &x,
                                 const RangeType &phi,
                                 JacobianRangeType &flux ) const
    {
      assert( Properties :: hasConvectiveFlux );
      flux = 0;
    }

    /** \copydoc Dune::LinearEllipticModelInterface::convectiveFlux(const EntityType &entity,const QuadratureType &quadrature,const int quadPoint,const RangeType &phi,JacobianRangeType &flux) const
     *
     *  The default implementation calls
     *  \code
     *  convectiveFlux( entity, quadrature.point( quadPoint ), phi, ret );
     *  \endcode
     */
    template< class EntityType, class QuadratureType >
    inline void convectiveFlux ( const EntityType &entity,
                                 const QuadratureType &quadrature,
                                 const int quadPoint,
                                 const RangeType &phi,
                                 JacobianRangeType &flux ) const
    {
      asImp().convectiveFlux( entity, quadrature.point( quadPoint ), phi, flux );
    }

    /** \copydoc Dune::LinearEllipticModelInterface::mass(const EntityType &entity,const DomainType &x,RangeType &ret) const
     *
     *  The default implementation returns 0.
     */
    template< class EntityType >
    inline void mass ( const EntityType &entity,
                       const DomainType &x,
                       RangeType &ret ) const
    {
      assert( Properties :: hasMass );
      ret = 0;
    }

    /** \copydoc Dune::LinearEllipticModelInterface::mass(const EntityType &entity,const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const
     *
     *  The default implementation calls
     *  \code
     *  mass( entity, quadrature.point( quadPoint ), ret );
     *  \endcode
     */
    template< class EntityType, class QuadratureType >
    inline void mass ( const EntityType &entity,
                       const QuadratureType &quadrature,
                       const int quadPoint,
                       RangeType &ret ) const
    {
      asImp().mass( entity, quadrature.point( quadPoint ), ret );
    }
    
    /** \copydoc Dune::LinearEllipticModelInterface::source(const EntityType &entity,const DomainType &x,RangeType &ret) const
     *
     *  The default implementation returns 0.
     */
    template< class EntityType >
    inline void source ( const EntityType &entity,
                         const DomainType &x,
                         RangeType &ret ) const
    {
      assert( Properties :: hasSource );
      ret = 0;
    }
    
    /** \copydoc Dune::LinearEllipticModelInterface::source(const EntityType &entity,const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const
     *
     *  The default implementation calls
     *  \code
     *  source( entity, quadrature.point( quadPoint ), ret );
     *  \endcode
     */
    template< class EntityType, class QuadratureType >
    inline void source ( const EntityType &entity,
                         const QuadratureType &quadrature,
                         const int quadPoint,
                         RangeType &ret ) const
    {
      asImp().source( entity, quadrature.point( quadPoint ), ret );
    }
  };

}

#endif
