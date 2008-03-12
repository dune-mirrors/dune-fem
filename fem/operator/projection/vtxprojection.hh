#ifndef DUNE_VTXPROJECTION_HH
#define DUNE_VTXPROJECTION_HH

#include <dune/fem/space/lagrangespace/lagrangespace.hh>
// #include <dune/grid/utility/twistutility.hh>
#include <dune/fem/operator/common/operator.hh>
namespace Dune 
{

template <class GridPartType>
struct WeightDefault {
  typedef FieldVector<double,GridPartType::GridType::dimension> DomainType;
  template <class EntityType>
  void setEntity(const EntityType& en) {
    volume_ = en.geometry().volume();
  }
  double operator()(const DomainType& point) {
    return volume_;
  }
  double volume_;
};
  
struct VtxProjectionImpl
{
  template <class ArgFunctionImp, class DiscreteFunctionImp,class WeightType>
  static void project(const ArgFunctionImp& f, DiscreteFunctionImp& discFunc,
                      WeightType& weight) 
  {
    typedef typename ArgFunctionImp::FunctionSpaceType ArgFunctionSpaceType;
    typedef typename ArgFunctionImp::LocalFunctionType ArgLocalFuncType;
    typedef typename DiscreteFunctionImp::FunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionImp::LocalFunctionType LocalFuncType;
    typedef typename DiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::Traits::IteratorType Iterator;
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType ; 
    // typedef typename GridPartType::GridType GridType;

    typename ArgFunctionSpaceType::RangeType val;
    const unsigned int dimRange = ArgFunctionSpaceType :: dimRange;
    const DiscreteFunctionSpaceType& space =  discFunc.space();
    typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
            LagrangePointSetType;
    typedef typename LagrangePointSetType :: template Codim< 0 >
                     :: SubEntityIteratorType
            EntityDofIteratorType;

    discFunc.clear();
    DiscreteFunctionImp weightDF("weight",space);
    weightDF.clear();

    Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      typename Iterator::Entity& en = *it;
      weight.setEntity(en);
      LocalFuncType lw  = weightDF.localFunction(en);
      LocalFuncType ldf = discFunc.localFunction(en);
      const ArgLocalFuncType larg = f.localFunction(en);
      const LagrangePointSetType &lagrangePointSet
                = space.lagrangePointSet( en );
      EntityDofIteratorType itPoint
             = lagrangePointSet.template beginSubEntity< 0 >( 0 );
      const EntityDofIteratorType enditPoint
            = lagrangePointSet.template endSubEntity< 0 >( 0 );
      for( ; itPoint != enditPoint; ++itPoint ) {
        const unsigned int dof = *itPoint;
        const typename ArgFunctionSpaceType::DomainType &point = lagrangePointSet.point( dof );
        larg.evaluate(point, val);
        double w = weight(point);
        val *= w;
        for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate ) {
          ldf[ dimRange * dof + coordinate ] += val[coordinate];
          lw[  dimRange * dof + coordinate ] += w;
        }
      }
    }
    typename DiscreteFunctionImp::DofIteratorType
      itdof = discFunc.dbegin();
    const typename DiscreteFunctionImp::DofIteratorType
      enddof = discFunc.dend();
    typename DiscreteFunctionImp::DofIteratorType
      itwdof = weightDF.dbegin();
    // const typename DiscreteFunctionImp::DofIteratorType
    //   endwdof = weightDF.dend();
    for (;itdof != enddof;++itdof,++itwdof) {
      if (*itwdof>0) {
        *itdof /= *itwdof;
      }
      assert(*itwdof>0 || *itdof == 0);
    }

  }
};

/*======================================================================*/
/*! \ingroup VtxProjectionOperator
 *  \class VtxProjection
 *  \brief The Projection class which average
 *         discontinuous function in the Lagrangepoints
 */
/*======================================================================*/
template <typename DFieldType, typename RFieldType,
          typename DType , typename RType>
class VtxProjection : public Operator<DFieldType, RFieldType,DType , RType> {
 public:
  typedef DType DomainType;
  typedef RType  RangeType;
  typedef DFieldType DomainFieldType;
  typedef RFieldType RangeFieldType;
  typedef typename RType::DiscreteFunctionSpaceType::GridPartType GridPartType;
  //! Constructor 
  VtxProjection() {}

  //! apply projection
  template <class WeightType>
  void operator() (const DomainType& f, RangeType& discFunc,
                   WeightType& weight) const 
  {
    VtxProjectionImpl::project(f,discFunc,weight);
  }
  //! apply projection with default weight
  void operator() (const DomainType& f, RangeType& discFunc) const
  {
    WeightDefault<GridPartType> weight; 
    VtxProjectionImpl::project(f,discFunc, weight);
  }

private:
};

}
#endif
