#ifndef DUNE_EXAMPLES_LAPLACE_MODEL_HH
#define DUNE_EXAMPLES_LAPLACE_MODEL_HH

namespace Dune
{

  template< class Problem >
  class LaplaceModel
  : public DiffusionModelDefault
    < typename Problem :: FunctionSpaceType, LaplaceModel< Problem > >,
    public BoundaryModelDefault
    < typename Problem :: FunctionSpaceType, LaplaceModel< Problem > >
  {
  public:
    typedef Problem ProblemType;
    
    typedef typename ProblemType :: FunctionSpaceType FunctionSpaceType;

  private:
    typedef DiffusionModelDefault< FunctionSpaceType, LaplaceModel >
      DiffusionModelBaseType;
    typedef BoundaryModelDefault< FunctionSpaceType, LaplaceModel >
      BoundaryModelBaseType;

  public:
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename BoundaryModelBaseType :: BoundaryType BoundaryType;

  public:
    using DiffusionModelBaseType :: diffusiveFlux;

  public:
    template< class EntityType, class PointType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const PointType &x,
                                const JacobianRangeType &gradient,
                                JacobianRangeType &flux ) const
    {
      // since we want to model the Laplace equation, just the identical flux.
      flux = gradient;
    }

    template< class IntersectionIteratorType >
    inline BoundaryType boundaryType ( const IntersectionIteratorType &intersection ) const
    {
      return BoundaryModelBaseType :: Dirichlet;
    }

    template< class IntersectionIteratorType, class QuadratureType >
    inline void dirichletValues ( const IntersectionIteratorType &intersection,
                                  const QuadratureType &quadrature,
                                  int point,
                                  RangeType &phi ) const
    {
      phi = 0;
    }
  };

}

#endif
