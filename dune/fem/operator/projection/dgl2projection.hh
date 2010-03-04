#ifndef DUNE_DGL2PROJECTION_HH
#define DUNE_DGL2PROJECTION_HH

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/discretefunctionadapter.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>

namespace Dune 
{

// implementation of L2 projection for discontinuous spaces 
class DGL2ProjectionImpl
{
  template <int dummy, bool hasLocalFunction> 
  struct ProjectChooser
  {
    template <class FunctionImp, class FunctionSpace> 
    class FunctionAdapter
    { 
      const FunctionImp& function_;
    public:  
      typedef FunctionSpace FunctionSpaceType;
      typedef typename FunctionSpaceType :: RangeType RangeType;
      typedef typename FunctionSpaceType :: DomainType DomainType;
      FunctionAdapter(const FunctionImp& f) : function_(f) {}

      void evaluate(const DomainType& local,
                    RangeType& ret) const 
      {
        function_.evaluate( local , ret );
      }
    };
    
    template <class FunctionImp, class DiscreteFunctionImp>
    static void project(const FunctionImp& f, 
                        DiscreteFunctionImp& discFunc,
                        int polOrd) 
    {
      // some typedefs 
      typedef typename DiscreteFunctionImp :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
      typedef FunctionAdapter<FunctionImp, typename  DiscreteFunctionSpaceType :: FunctionSpaceType> FunctionAdapterType;
      // create function adapter in case of incorrect implementation 
      FunctionAdapterType af( f );
      // create discrete function adapter 
      DiscreteFunctionAdapter< FunctionAdapterType, GridPartType> adapter(
          "DGL2projection::adapter" , f , discFunc.space().gridPart());
      
      DGL2ProjectionImpl::projectFunction(adapter, discFunc, polOrd);
    }
  };

  template <int dummy> 
  struct ProjectChooser<dummy,true>
  {
    template <class FunctionImp, class DiscreteFunctionImp>
    static void project(const FunctionImp& f, 
                        DiscreteFunctionImp& discFunc,
                        int polOrd) 
    {
      DGL2ProjectionImpl::projectFunction(f, discFunc, polOrd);
    }
  };

public:  
  //! project function onto discrete function space  
  template <class FunctionImp, class DiscreteFunctionImp>
  static void project(const FunctionImp& f, DiscreteFunctionImp& discFunc, int polOrd = -1) 
  {
    ProjectChooser<0, Conversion<FunctionImp, HasLocalFunction> ::exists > :: project(f,discFunc,polOrd);
  }

protected:  
  template <class FunctionImp, class DiscreteFunctionImp>
  static void projectFunction(const FunctionImp& func, 
                              DiscreteFunctionImp& discFunc, 
                              int polOrd = -1) 
  {
    typedef typename DiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionImp::LocalFunctionType LocalFuncType;
    typedef typename DiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::Traits::IteratorType Iterator;
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType ; 
    typedef typename GridPartType::GridType GridType;

    typedef typename FunctionImp::LocalFunctionType LocalFType;
    
    const DiscreteFunctionSpaceType& space =  discFunc.space();

    // type of quadrature 
    typedef CachingQuadrature<GridPartType,0> QuadratureType; 
    // type of local mass matrix 
    typedef LocalDGMassMatrix< DiscreteFunctionSpaceType, QuadratureType > LocalMassMatrixType;

    const int quadOrd = (polOrd == -1) ? (2 * space.order()) : polOrd;
    
    // create local mass matrix object
    LocalMassMatrixType massMatrix( space, quadOrd );

    // check whether geometry mappings are affine or not 
    const bool affineMapping = massMatrix.affine();

    // clear destination
    discFunc.clear();

    const Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      // get entity 
      const typename GridType::template Codim<0>::Entity& en = *it; 
      // get geometry 
      const typename GridType::template Codim<0>::Geometry& geo = en.geometry(); 
      
      // get quadrature 
      QuadratureType quad(en, quadOrd);
      
      // get local function of destination 
      LocalFuncType lf = discFunc.localFunction(en);
      // get local function of argument 
      const LocalFType f = func.localFunction(en);

      const int quadNop = quad.nop();

      typename DiscreteFunctionSpaceType :: RangeType value ;

      for(int qP = 0; qP < quadNop ; ++qP) 
      {
        const double intel = (affineMapping) ? 
             quad.weight(qP) : // affine case 
             quad.weight(qP) * geo.integrationElement( quad.point(qP) ); // general case 

        // evaluate function 
        f.evaluate(quad[ qP ], value );
        
        // apply weight 
        value *= intel;

        // add to local function 
        lf.axpy( quad[ qP ], value );
      }

      // in case of non-linear mapping apply inverse 
      if ( ! affineMapping ) 
      {
        massMatrix.applyInverse( en, lf );
      } 
    }
  }
};

} // end namespace Dune 
#endif
