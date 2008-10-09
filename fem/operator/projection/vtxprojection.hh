#ifndef DUNE_VTXPROJECTION_HH
#define DUNE_VTXPROJECTION_HH

#include <dune/fem/space/lagrangespace/lagrangespace.hh>
// #include <dune/grid/utility/twistutility.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/space/common/communicationmanager.hh>
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
    typedef CommunicationManager< DiscreteFunctionSpaceType > CommunicationManagerType;
    // typedef typename GridPartType::GridType GridType;

    typename ArgFunctionSpaceType::RangeType val;
    typedef typename ArgFunctionSpaceType::DomainType DomainType;
    const unsigned int dimRange = ArgFunctionSpaceType :: dimRange;
    const DiscreteFunctionSpaceType& space =  discFunc.space();
    typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
            LagrangePointSetType;
    typedef typename LagrangePointSetType :: template Codim< 0 >
                     :: SubEntityIteratorType
            EntityDofIteratorType;
    typedef typename LagrangePointSetType :: template Codim< 1 >
                     :: SubEntityIteratorType
            FaceDofIteratorType;

    typedef typename Iterator::Entity EntityType;

    discFunc.clear();
    DiscreteFunctionImp weightDF("weight",space);
    weightDF.clear();

    const Iterator endit = space.end();
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
    
    discFunc.space().communicate( discFunc );
    weightDF.space().communicate( weightDF );

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

    // make function continuous over hanging nodes

    if (!GridPartType::conforming) {
      const GridPartType& gridPart =  space.gridPart();
      for(Iterator it = space.begin(); it != endit ; ++it) {
        const EntityType& en = *it;
        typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
        const IntersectionIteratorType endnit = gridPart.iend(en);
        for (IntersectionIteratorType nit = gridPart.ibegin(en); nit != endnit; ++nit) 
        {
          typedef typename IntersectionIteratorType::Intersection IntersectionType;
          const IntersectionType& inter=*nit;
          if (inter.neighbor()) 
          {
            // get neighbor 
            typename EntityType::EntityPointer ep = inter.outside();
            EntityType & nb = *ep;
            if (en.level()>nb.level()) {
              const int numInSelf = inter.numberInSelf();
              const LagrangePointSetType &lagrangePointSet
                    = space.lagrangePointSet( en );
              FaceDofIteratorType itPoint
                 = lagrangePointSet.template beginSubEntity< 1 >( numInSelf );
              const FaceDofIteratorType enditPoint
                 = lagrangePointSet.template endSubEntity< 1 >( numInSelf );
              const typename IntersectionType::LocalGeometry& geoIn  = inter.intersectionSelfLocal();
              const typename IntersectionType::LocalGeometry& geoOut = inter.intersectionNeighborLocal();
              LocalFuncType ldfIn  = discFunc.localFunction(en);
              LocalFuncType ldfOut = discFunc.localFunction(nb);
              for( ; itPoint != enditPoint; ++itPoint ) {
                const unsigned int dof = *itPoint;
                const DomainType &point = lagrangePointSet.point( dof );
                DomainType x = geoOut.global(geoIn.local(point));
                ldfOut.evaluate(x, val);
                for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate ) {
                  ldfIn[ dimRange * dof + coordinate ] = val[coordinate];
                }
              }
            }
          }
        }
      }
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
