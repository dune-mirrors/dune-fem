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
   *  -\nabla \cdot \bigl( D( x ) \nabla u \bigr).
   *  \f]
   *  Here the matematical data is exactly the function \f$D\f$, which we call
   *  diffusive flux.
   *
   *  \param FunctionSpaceImp type of the function space we are using
   *  \param DiffusionModelImp implementation of the interface (Barton-Nackman)
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
    typedef ThisType DiffusionModelInterfaceType;

    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  public:
    /*! \brief calculate the diffusive flux in a quadrature point
     *
     *  \param[in] entity entity on which the flux shall be calculated
     *  \param[in] quadrature quadrature to be used
     *  \param[in] pt index of the current point within the quadrature
     *  \param[in] gradient \f$\nabla u\f$ evaluated in the quadrature point
     *  \param[out] flux calculated diffusive flux \f$D( x ) \nabla u\f$.
     */
    template< class EntityType, class QuadratureType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const QuadratureType &quadrature,
                                unsigned int pt,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().diffusiveFlux( entity, quadrature, pt, gradient, flux ) );
    }

  protected:
    inline const DiffusionModelImp& asImp () const
    {
      return static_cast< const DiffusionModelImp& >( *this );
    }

    inline DiffusionModelImp& asImp ()
    {
      return static_cast< DiffusionModelImp& >( *this );
    }
  };

}

#endif
