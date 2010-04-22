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

    typedef typename DiscreteFunctionImp::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionImp::LocalFunctionType LocalFunctionType;

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

    typedef typename Iterator::Entity              EntityType;
    typedef typename EntityType :: EntityPointer   EntityPointer;
    typedef typename EntityType :: Geometry        Geometry;
    typedef typename Geometry :: LocalCoordinate   LocalCoordinate;

    discFunc.clear();
    DiscreteFunctionImp weightDF("weight",space);
    weightDF.clear();

    const Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      const EntityType& en = *it;
      weight.setEntity(en);

      LocalFunctionType lw  = weightDF.localFunction(en);
      LocalFunctionType ldf = discFunc.localFunction(en);

      const ArgLocalFuncType larg = f.localFunction(en);
      const LagrangePointSetType &lagrangePointSet = space.lagrangePointSet( en );

      const unsigned int numPoints = lagrangePointSet.nop();
      for( unsigned int pt = 0; pt < numPoints; ++pt )
      {
        larg.evaluate( lagrangePointSet[ pt ], val );
        double w = weight( lagrangePointSet.point( pt ) );
        for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
        {
          ldf[ dimRange * pt + coordinate ] += w * val[ coordinate];
          lw[ dimRange * pt + coordinate ] += w;
        }
      }
    }
    
    discFunc.communicate();
    weightDF.communicate();

    typedef typename DiscreteFunctionImp :: DofIteratorType  DofIteratorType;

    DofIteratorType        itdof = discFunc.dbegin();
    const DofIteratorType enddof = discFunc.dend();
    DofIteratorType       itwdof = weightDF.dbegin();

    for (;itdof != enddof;++itdof,++itwdof) 
    {
      if (*itwdof>0) 
      {
        *itdof /= *itwdof;
      }
      assert(*itwdof>0 || *itdof == 0);
    }

    // make function continuous over hanging nodes

    if( ! GridPartType::conforming )
    {
      const GridPartType &gridPart =  space.gridPart();
      for( Iterator it = space.begin(); it != endit ; ++it )
      {
        const EntityType &entity = *it;
        typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
        const IntersectionIteratorType iend = gridPart.iend( entity );
        for( IntersectionIteratorType iit = gridPart.ibegin( entity ); iit != iend; ++iit )
        {
          typedef typename IntersectionIteratorType::Intersection IntersectionType;
          const IntersectionType &intersection = *iit;
          if( intersection.neighbor() )
          {
            // get neighbor 
            EntityPointer ep = intersection.outside();
            const EntityType&  neighbor = *ep;
            // if non-conforming situation 
            if( entity.level() > neighbor.level() )
            {
              const int indexInInside = intersection.indexInInside();
              const LagrangePointSetType &lagrangePointSet = space.lagrangePointSet( entity );
              FaceDofIteratorType itPoint
                 = lagrangePointSet.template beginSubEntity< 1 >( indexInInside );
              const FaceDofIteratorType enditPoint
                 = lagrangePointSet.template endSubEntity< 1 >( indexInInside );

              typedef typename IntersectionType :: LocalGeometry   LocalGeometry;
              const LocalGeometry &geoIn  = intersection.geometryInInside();
              const LocalGeometry &geoOut = intersection.geometryInOutside();

              LocalFunctionType ldfIn  = discFunc.localFunction( entity );
              LocalFunctionType ldfOut = discFunc.localFunction( neighbor );

              for( ; itPoint != enditPoint; ++itPoint )
              {
                const unsigned int dof = *itPoint;
                const LocalCoordinate& point = lagrangePointSet.point( dof );
                const LocalCoordinate x = geoOut.global( geoIn.local( point ) );

                ldfOut.evaluate( x, val );
                for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
                  ldfIn[ dimRange * dof + coordinate ] = val[ coordinate ];
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
