#ifndef DUNE_EXAMPLES_LAPLACE_SINEPROBLEM_HH
#define DUNE_EXAMPLES_LAPLACE_SINEPROBLEM_HH

namespace Dune
{

  template< class FunctionSpaceImp >
  class RightHandSide
  : public Function< FunctionSpaceImp, RightHandSide< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;
  
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    
  private:
    typedef RightHandSide< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    inline RightHandSide ( const FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }

    inline void evaluate ( const DomainType &x, RangeType &y ) const
    {
      enum { dimension = DomainType :: dimension };

      y = 0;
      for( int i = 0; i < dimension; ++i )
      {
        RangeType z = M_PI * M_PI;
        for( int j = 0; j < dimension; ++j )
        {
          const DomainFieldType &xj = x[ j ];
          z *= sin( M_PI * xj );
        }
        y += z;
      }
    }

    inline void evaluate ( const DomainType &x, RangeFieldType t, RangeType &y ) const
    {
      evaluate( x, y );
    }
  };



  template< class FunctionSpaceImp >
  class ExactSolution
  : public Function< FunctionSpaceImp, ExactSolution< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;
    
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  private:
    typedef ExactSolution< FunctionSpaceType > ThisType;
    typedef Function< FunctionSpaceType, ThisType > BaseType;

  public:
    inline ExactSolution ( FunctionSpaceType &functionSpace )
    : BaseType( functionSpace )
    {
    }

    inline void evaluate ( const DomainType &x, RangeType &y ) const
    {
      enum { dimension = DomainType :: dimension };
      
      y = 1;
      for( int i = 0; i < dimension; ++i )
      {
        const DomainFieldType &xi = x[ i ];
        y *= sin( M_PI * xi );
      }
    }

    inline void evaluate ( const DomainType &x, RangeFieldType t, RangeType &y ) const
    {
      evaluate( x, y );
    }
  };



  template< class FunctionSpaceImp >
  class LaplaceModel
  : public DiffusionModelDefault< FunctionSpaceImp, LaplaceModel< FunctionSpaceImp > >,
    public BoundaryModelDefault< FunctionSpaceImp, LaplaceModel< FunctionSpaceImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

  private:
    typedef LaplaceModel< FunctionSpaceType > ThisType;
    typedef DiffusionModelDefault< FunctionSpaceType, ThisType > DiffusionModelBaseType;
    typedef BoundaryModelDefault< FunctionSpaceType, ThisType > BoundaryModelBaseType;

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
    template< class EntityType >
    inline void diffusiveFlux ( const EntityType &entity,
                                const DomainType &x,
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
