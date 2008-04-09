#ifndef DUNE_L2PROJECTION_HH
#define DUNE_L2PROJECTION_HH

#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/discretefunctionadapter.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>

namespace Dune 
{

struct L2ProjectionImpl
{
  template <int dummy, bool isDiscreteFunction> 
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
          "L2projection::adapter" , f , discFunc.space().gridPart());
      
      L2ProjectionImpl::projectFunction(adapter, discFunc, polOrd);
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
      L2ProjectionImpl::projectFunction(f, discFunc, polOrd);
    }
  };

  template <class FunctionImp, class DiscreteFunctionImp>
  static void project(const FunctionImp& f, DiscreteFunctionImp& discFunc, int polOrd_ = -1) 
  {
    ProjectChooser<0, Conversion<FunctionImp, IsDiscreteFunction> ::exists > :: project(f,discFunc,polOrd_);
  }

protected:  
  template <class FunctionImp, class DiscreteFunctionImp>
  static void projectFunction(const FunctionImp& func, 
                              DiscreteFunctionImp& discFunc, 
                              int polOrd_ = -1) 
  {
    typedef typename DiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionImp::LocalFunctionType LocalFuncType;
    typedef typename DiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::Traits::IteratorType Iterator;
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType ; 
    typedef typename GridPartType::GridType GridType;

    typedef typename FunctionImp::LocalFunctionType LocalFType;
    
    typename DiscreteFunctionSpaceType::RangeType ret (0.0);
    typename DiscreteFunctionSpaceType::RangeType phi (0.0);
    const DiscreteFunctionSpaceType& space =  discFunc.space();

    // type of quadrature 
    typedef CachingQuadrature<GridPartType,0> QuadratureType; 
    // type of local mass matrix 
    typedef LocalDGMassMatrix< DiscreteFunctionSpaceType, QuadratureType > LocalMassMatrixType;

    const int quadOrd = (polOrd_ == -1) ? (2 * space.order()) : polOrd_;
    
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

      // get base function set 
      const BaseFunctionSetType & baseset = lf.baseFunctionSet();

      const int quadNop = quad.nop();
      const int numDofs = lf.numDofs();

      for(int qP = 0; qP < quadNop ; ++qP) 
      {
        const double intel = (affineMapping) ? 
             quad.weight(qP) : // affine case 
             quad.weight(qP) * geo.integrationElement( quad.point(qP) ); // general case 

        // evaluate function 
        f.evaluate(quad[qP], ret);
        
        // do projection 
        for(int i=0; i<numDofs; ++i) 
        {
          baseset.evaluate(i, quad[qP], phi);
          lf[i] += intel * (ret * phi) ;
        }
      }

      // in case of non-linear mapping apply inverse 
      if ( ! affineMapping ) 
      {
        massMatrix.applyInverse( en, lf );
      } 
    }
  }
};

/*======================================================================*/
/*! @ingroup L2ProjectionOperator
 *  \class L2Projection
 *  \brief The L2Projection class provides methods for projecting a function
 *         unto a given discrete function space. Note that this implementation
 *         assumes orthorgonal base functions!
 */
/*======================================================================*/
template <typename DFieldType, typename RFieldType,
          typename DType , typename RType>
class L2Projection : public Operator<DFieldType, RFieldType,DType , RType> {
 public:
  typedef DType DomainType;
  typedef RType  RangeType;
  typedef DFieldType DomainFieldType;
  typedef RFieldType RangeFieldType;

  //! Constructor taking degree for quadrature rule
  //! if no argument given a default value is chosen depending on the order
  //! in the discrete function space
  L2Projection(int polOrd = -1) : polOrd_(polOrd) {}

  //! apply L2 projection
  virtual void operator() (const DomainType& f, RangeType& discFunc) const 
  {
    L2ProjectionImpl::project(f,discFunc,polOrd_);
  }

private:
  const int polOrd_;
};

}
#endif
