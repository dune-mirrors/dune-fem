#ifndef DUNE_L2PROJECTION_HH
#define DUNE_L2PROJECTION_HH

#include <dune/fem/quadrature/cachequad.hh>
// #include <dune/grid/utility/twistutility.hh>
#include <dune/fem/operator/common/operator.hh>
namespace Dune 
{

struct L2ProjectionImpl
{
  template <class FunctionImp, class DiscreteFunctionImp>
  static void project(const FunctionImp& f, DiscreteFunctionImp& discFunc, int polOrd_ = -1) 
  {
    projectFunction(f,discFunc);
  }

  template <class FunctionImp, class DiscreteFunctionImp>
  static void projectDiscreteFunction(const FunctionImp& func, 
          DiscreteFunctionImp& discFunc, int polOrd_ = -1) 
  {
    typedef typename DiscreteFunctionImp::FunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionImp::LocalFunctionType LocalFuncType;
    typedef typename DiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::Traits::IteratorType Iterator;
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType ; 
    typedef typename GridPartType::GridType GridType;

    typedef typename FunctionImp::LocalFunctionType LocalFType;
    
    typename DiscreteFunctionSpaceType::RangeType ret (0.0);
    typename DiscreteFunctionSpaceType::RangeType phi (0.0);
    const DiscreteFunctionSpaceType& space =  discFunc.space();

    int polOrd = polOrd_;
    if (polOrd == -1) 
      polOrd = 2 * space.order();

    discFunc.clear();

    Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      typename GridType::template Codim<0>::Entity& en = *it; 
      
      // Get quadrature rule
      CachingQuadrature<GridPartType,0> quad(en, polOrd);
      LocalFuncType lf = discFunc.localFunction(en);
      LocalFType f = func.localFunction(en);

      //! Note: BaseFunctions must be ortho-normal!!!!
      const BaseFunctionSetType & baseset = lf.baseFunctionSet();

      const int quadNop = quad.nop();
      const int numDofs = lf.numDofs();

      for(int qP = 0; qP < quadNop ; ++qP) 
      {
        f.evaluate(quad,qP, ret);
        for(int i=0; i<numDofs; ++i) 
        {
          baseset.evaluate(i,quad,qP,phi);
          lf[i] += quad.weight(qP) * (ret * phi) ;
        }
      }
    }
  }
  
  template <class FunctionImp, class DiscreteFunctionImp>
  static void projectFunction(const FunctionImp& f, DiscreteFunctionImp& discFunc, int polOrd_ = -1) 
  {
    typedef typename DiscreteFunctionImp::FunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionImp::LocalFunctionType LocalFuncType;
    typedef typename DiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::Traits::IteratorType Iterator;
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType ; 
    typedef typename GridPartType::GridType GridType;

    typename DiscreteFunctionSpaceType::RangeType ret (0.0);
    typename DiscreteFunctionSpaceType::RangeType phi (0.0);
    const DiscreteFunctionSpaceType& space =  discFunc.space();

    int polOrd = polOrd_;
    if (polOrd == -1) 
      polOrd = 2 * space.order();

    discFunc.clear();

    Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      typename GridType::template Codim<0>::Entity& en = *it; 
      
      // Get quadrature rule
      CachingQuadrature<GridPartType,0> quad(en, polOrd);
      LocalFuncType lf = discFunc.localFunction(en);

      //! Note: BaseFunctions must be ortho-normal!!!!
      const BaseFunctionSetType & baseset = lf.baseFunctionSet();

      const typename GridType::template Codim<0>::Entity::Geometry& 
        itGeom = en.geometry();
     
      const int quadNop = quad.nop();
      const int numDofs = lf.numDofs();

      for(int qP = 0; qP < quadNop ; ++qP) 
      {
        f.evaluate(itGeom.global(quad.point(qP)), ret);
        for(int i=0; i<numDofs; ++i) 
        {
          baseset.evaluate(i,quad,qP,phi);
          lf[i] += quad.weight(qP) * (ret * phi) ;
        }
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
