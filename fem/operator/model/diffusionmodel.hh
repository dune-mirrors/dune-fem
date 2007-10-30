#ifndef DUNE_FEM_DIFFUSIONMODEL_HH
#define DUNE_FEM_DIFFUSIONMODEL_HH

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  /*! \ingroup FEM
   *  \class DiffusionModelInterface
   *  \brief Interface for a mathematical model of a DiffusionOperator
   *
   *  The DiffusionModelInterface specifies the way mathematical data can
   *  be incorporated into a DiffusionOperator, i.e. a linear operator
   *  performing
   *  \f[
   *  -\nabla \cdot \bigl( a( x ) \nabla u \bigr).
   *  \f]
   *  Here the matematical data is exactly the function \f$a\f$, which we call
   *  diffusive flux.
   *
   *  \param  FunctionSpaceImp   function space to work in
   *  \param  DiffusionModelImp  actual implementation (Barton-Nackman)
   */
  template< class FunctionSpaceImp, class DiffusionModelImp >
  class DiffusionModelInterface
  {
  public:
    //! type of the function space we are using
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef DiffusionModelInterface< FunctionSpaceType, DiffusionModelImp > ThisType;

  public:
    //! type of the interface
    typedef ThisType DiffusionModelInterfaceType;

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
    /** \brief calculate the diffusive flux \f$a( x ) \nabla u\f$ in a point
     *
     *  \param[in]  entity      entity to evaluate the flux on
     *  \param[in]  x           evaluation point (in local coordinates)
     *  \param[in]  gradient    \f$\nabla u\f$ in the evaluation point
     *  \param[out] flux        variable to receive the evaluated flux
     */
    template< class EntityType, class PointType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const PointType &x,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().diffusiveFlux( entity, x, gradient, flux ) );
    }

    /** \brief calculate the diffusive flux \f$a( x ) \nabla u\f$ in a quadrature
     *         point
     *
     *  \param[in]  entity      entity to evaluate the flux on
     *  \param[in]  quadrature  quadrature to use
     *  \param[in]  quadPoint   number of the point within the quadrature
     *  \param[in]  gradient    \f$\nabla u\f$ in the evaluation point
     *  \param[out] flux        variable to receive the evaluated flux
     */
    template< class EntityType, class QuadratureType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const QuadratureType &quadrature,
                                const int quadPoint,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().diffusiveFlux( entity, quadrature, quadPoint, gradient, flux ) );
    }

  protected:
    inline const DiffusionModelImp &asImp () const
    {
      return static_cast< const DiffusionModelImp & >( *this );
    }

    inline DiffusionModelImp &asImp ()
    {
      return static_cast< DiffusionModelImp & >( *this );
    }
  };


  
  template< class FunctionSpaceImp, class DiffusionModelImp >
  class DiffusionModelDefault
  : public DiffusionModelInterface< FunctionSpaceImp, DiffusionModelImp >
  {
  public:
    //! type of the function space we are using
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef DiffusionModelDefault< FunctionSpaceType, DiffusionModelImp > ThisType;
    typedef DiffusionModelInterface< FunctionSpaceType, DiffusionModelImp > BaseType;

  public:
    //! type of the interface
    typedef ThisType DiffusionModelInterfaceType;

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
    using BaseType :: diffusiveFlux;

  protected:
    using BaseType :: asImp;

  public:
    /** \copydoc Dune::DiffusionModelInterface::diffusiveFlux(const EntityType &entity,const PointType &x,const JacobianRangeType &gradient,JacobianRangeType &flux) const */
    template< class EntityType, class QuadratureType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const QuadraturePointWrapper< QuadratureType > &x,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux ) const
    {
      asImp().diffusiveFlux( entity, x.quadrature(), x.point(), gradient, flux );
    }

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
  };

}

#endif
