#ifndef __DUNE_FEM_CC__
#define __DUNE_FEM_CC__

// assembling of the laplace operator using the 
// below defined getLocalMatrixMethod of the LaplaceOperator 

#include "laplace.hh" 

// where the quadratures are defined 
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/quadrature/quadrature.hh>

namespace Dune 
{

// L2 projection for RHS of laplace's equation
// note: this is not a real L2-Projection, though
template< class DiscreteFunctionImp >
class L2Projection
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
    
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType :: LocalFunctionType
    LocalFunctionType;
  
  typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
    BaseFunctionSetType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename GridPartType :: GridType GridType;

  enum { dimension = GridType :: dimension };
  
public:  
  template< int polOrd, class FunctionType >
  void project( FunctionType &function,
                DiscreteFunctionType &discreteFunction )
  {
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
    typedef typename EntityType :: Geometry GeometryType;
      
    const DiscreteFunctionSpaceType &discreteFunctionSpace
      = discreteFunction.space();
  
    discreteFunction.clear();

    const IteratorType endit = discreteFunctionSpace.end();
    for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
    {
      const GeometryType &geometry = (*it).geometry();
      
      LocalFunctionType localFunction = discreteFunction.localFunction( *it );
      const BaseFunctionSetType &baseFunctionSet
        = discreteFunctionSpace.baseFunctionSet( *it );

      CachingQuadrature< GridPartType, 0 > quadrature( *it, polOrd );
      const int numDofs = localFunction.numDofs();
      for( int i = 0; i < numDofs; ++i )
      {
        RangeType phi, psi;
        
        const int numQuadraturePoints = quadrature.nop();
        for( int qP = 0; qP < numQuadraturePoints; ++qP )
        {
          const double det
            = geometry.integrationElement( quadrature.point( qP ) );
          function.evaluate( geometry.global( quadrature.point( qP ) ), phi );
          baseFunctionSet.eval( i, quadrature, qP, psi );
          localFunction[ i ] += det * quadrature.weight( qP ) * (phi * psi);
        }
      }
    }
  }
};

#if 0
// used for calculation of the initial values 
template <class DiscreteFunctionType> 
class LagrangeInterpolation
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  
public:  
  template <class FunctionType> 
  void interpol (FunctionType &f, DiscreteFunctionType &discFunc)
  {
    const typename DiscreteFunctionType::FunctionSpaceType
        & space = discFunc.getFunctionSpace();  
  
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename FunctionSpaceType:: IteratorType IteratorType; 
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
      
    typedef typename FunctionSpaceType::RangeType RangeType;
    RangeType ret (0.0);
    //LocalFuncType lf = discFunc.newLocalFunction(); 

    IteratorType endit = space.end(); 
    for(IteratorType it = space.begin(); it != endit ; ++it)
    {
      //discFunc.localFunction(*it,lf); 
      LocalFuncType lf = discFunc.localFunction(*it); 
      int numDof = lf.numDofs ();  
      for(int i=0; i<numDof; i++)
      {
        f.evaluate(it->geometry()[i],ret);
        lf[i] = ret[0];
      }
    }
  }
};
#endif

} // end namespace 

#endif

