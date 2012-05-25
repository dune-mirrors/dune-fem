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
  template< class Function, class DiscreteFunction, class Weight >
  static void project ( const Function &f, DiscreteFunction &u, Weight &weight )
  {
    typedef typename Function::FunctionSpaceType FunctionSpaceType;
    typedef typename Function::LocalFunctionType LocalFunctionType;

    typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunction::LocalFunctionType LocalDiscreteFunctionType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename DiscreteFunctionSpaceType::LagrangePointSetType LagrangePointSetType;

    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridPartType::IntersectionType IntersectionType;
    typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
    typedef typename GridPartType::template Codim< 0 >::EntityPointerType EntityPointerType;
    typedef typename GridPartType::template Codim< 0 >::GeometryType GeometryType;;

    typedef typename LagrangePointSetType::template Codim< 0 >::SubEntityIteratorType EntityDofIteratorType;
    typedef typename LagrangePointSetType::template Codim< 1 >::SubEntityIteratorType FaceDofIteratorType;

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename GeometryType::LocalCoordinate LocalCoordinateType;

    const unsigned int dimRange = FunctionSpaceType::dimRange;
    const DiscreteFunctionSpaceType &space = u.space();

    u.clear();
    DiscreteFunction w ( "weight", space );
    w.clear();

    const IteratorType end = space.end();
    for ( IteratorType it = space.begin(); it != end; ++it )
    {
      const EntityType &entity = *it;
      weight.setEntity( entity );

      LocalDiscreteFunctionType lw = w.localFunction( entity );
      LocalDiscreteFunctionType lu = u.localFunction( entity );

      const LocalFunctionType lf = f.localFunction( entity );
      const LagrangePointSetType &lagrangePointSet = space.lagrangePointSet( entity );

      const unsigned int numPoints = lagrangePointSet.nop();
      for( unsigned int pt = 0; pt < numPoints; ++pt )
      {
        RangeType val;
        lf.evaluate( lagrangePointSet[ pt ], val );

        double wght = weight( lagrangePointSet.point( pt ) );

        for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
        {
          lu[ dimRange*pt + coordinate ] += wght*val[ coordinate ];
          lw[ dimRange*pt + coordinate ] += wght;
        }
      }
    }
    
    u.communicate();
    w.communicate();

    typedef typename DiscreteFunction::DofIteratorType DofIteratorType;

    const DofIteratorType udend = u.dend();
    DofIteratorType udit = u.dbegin();
    DofIteratorType wdit = w.dbegin();
    for( ; udit != udend; ++udit, ++wdit )
    {
      assert( (*wdit > 0) || (*udit == 0) );
      *udit /= (*wdit > 0 ? *wdit : RangeFieldType( 1 ));
    }

    // make function continuous over hanging nodes

    if( !GridPartType::conforming )
    {
      const GridPartType &gridPart =  space.gridPart();
      for( IteratorType it = space.begin(); it != end; ++it )
      {
        const EntityType &entity = *it;

        const LagrangePointSetType &lagrangePointSet = space.lagrangePointSet( entity );

        const IntersectionIteratorType iend = gridPart.iend( entity );
        for( IntersectionIteratorType iit = gridPart.ibegin( entity ); iit != iend; ++iit )
        {
          const IntersectionType &intersection = *iit;

          if( intersection.neighbor() )
          {
            // get neighbor 
            EntityPointerType outside = intersection.outside();
            const EntityType &neighbor = *outside;

            // if non-conforming situation 
            if( entity.level() > neighbor.level() )
            {
              const int indexInInside = intersection.indexInInside();

              typedef typename IntersectionType::LocalGeometry LocalGeometryType;
              const LocalGeometryType &geoIn  = intersection.geometryInInside();
              const LocalGeometryType &geoOut = intersection.geometryInOutside();

              LocalDiscreteFunctionType uIn  = u.localFunction( entity );
              LocalDiscreteFunctionType uOut = u.localFunction( neighbor );

              const FaceDofIteratorType fdend = lagrangePointSet.template endSubEntity< 1 >( indexInInside );
              FaceDofIteratorType fdit = lagrangePointSet.template beginSubEntity< 1 >( indexInInside );
              for( ; fdit != fdend; ++fdit )
              {
                const LocalCoordinateType &xIn = lagrangePointSet.point( *fdit );
                const LocalCoordinateType xOut = geoOut.global( geoIn.local( xIn ) );

                RangeType val;
                uOut.evaluate( xOut, val );

                for( unsigned int coordinate = 0; coordinate < dimRange; ++coordinate )
                  uIn[ dimRange*(*fdit) + coordinate ] = val[ coordinate ];
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
