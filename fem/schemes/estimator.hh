/**************************************************************************

  The dune-fem module is a module of DUNE (see www.dune-project.org).
  It is based on the dune-grid interface library
  extending the grid interface by a number of discretization algorithms
  for solving non-linear systems of partial differential equations.

  Copyright (C) 2003 - 2015 Robert Kloefkorn
  Copyright (C) 2003 - 2010 Mario Ohlberger
  Copyright (C) 2004 - 2015 Andreas Dedner
  Copyright (C) 2005        Adrian Burri
  Copyright (C) 2005 - 2015 Mirko Kraenkel
  Copyright (C) 2006 - 2015 Christoph Gersbacher
  Copyright (C) 2006 - 2015 Martin Nolte
  Copyright (C) 2011 - 2015 Tobias Malkmus
  Copyright (C) 2012 - 2015 Stefan Girke
  Copyright (C) 2013 - 2015 Claus-Justus Heine
  Copyright (C) 2013 - 2014 Janick Gerstenberger
  Copyright (C) 2013        Sven Kaulman
  Copyright (C) 2013        Tom Ranner
  Copyright (C) 2015        Marco Agnese
  Copyright (C) 2015        Martin Alkaemper


  The dune-fem module is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The dune-fem module is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

**************************************************************************/
/**Residual estimator for Poisson's equation. Only for Dirichlet/periodic boundaries ATM. */
#ifndef ESTIMATOR_HH
#define ESTIMATOR_HH

//- Dune-fem includes
#include <dune/fem/misc/compatibility.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/operator/matrix/blockmatrix.hh>

#include "diffusionmodel.hh"

// Estimator
// ---------
template< class DiscreteFunction, class Model >
class Estimator
{
  typedef Estimator< DiscreteFunction, Model > ThisType;

public:
  typedef DiscreteFunction DiscreteFunctionType;

  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

  typedef typename GridPartType :: GridType GridType;
  typedef typename GridPartType :: IndexSetType IndexSetType;
  typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;

  typedef typename IntersectionIteratorType :: Intersection IntersectionType;

  typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
  typedef typename GridType :: template Codim< 0 > :: EntityPointer
    ElementPointerType;
  typedef typename ElementType::Geometry GeometryType;
  static const int dimension = GridType :: dimension;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > ElementQuadratureType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 1 > FaceQuadratureType;

  //typedef Dune::ElementQuadrature< GridPartType, 0 > ElementQuadratureType;
  //typedef Dune::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;

  typedef Dune::FieldMatrix<double,dimension,dimension> JacobianInverseType;
  typedef std :: vector< double > ErrorIndicatorType;

  // include the model
  typedef Model ModelType;

private:
  const DiscreteFunctionType &uh_;
  const DiscreteFunctionSpaceType &dfSpace_;
  GridPartType &gridPart_;
  const IndexSetType &indexSet_;
  GridType &grid_;
  ErrorIndicatorType indicator_;
  const ModelType& model_;

public:
  static_assert( static_cast< unsigned int >( GridType::dimension ) == static_cast< unsigned int >( GridType::dimensionworld ),
                 "the estimator is not implemented for surfaces problems" );

  explicit Estimator ( const DiscreteFunctionType &uh, const ModelType& model )
  : uh_( uh ),
    dfSpace_( uh.space() ),
    gridPart_( dfSpace_.gridPart() ),
    indexSet_( gridPart_.indexSet() ),
    grid_( gridPart_.grid() ),
    indicator_( indexSet_.size( 0 ) ),
    model_(model)
  {}

  //! calculate estimator
  template< class RHSFunctionType >
  double estimate ( const RHSFunctionType &rhs )
  {
    // clear all local estimators
    clear();

    //! [Error estimator]
    // calculate local estimator
    const IteratorType end = dfSpace_.end();
    for( IteratorType it = dfSpace_.begin(); it != end; ++it )
      estimateLocal( rhs, *it );

    double error = 0.0;

    // sum up local estimators
    const typename ErrorIndicatorType :: const_iterator endEstimator = indicator_.end();
    for( typename ErrorIndicatorType :: const_iterator it = indicator_.begin();
          it != endEstimator; ++it )
      error += *it;

    // obtain global sum
    error = grid_.comm().sum( error );

    return sqrt( error );
    //! [Error estimator]
  }

  //! mark all elements due to given tolerance
  bool mark ( const double tolerance ) const
  {
    // possible strategies
    enum Strategy { none=0, maximum=1, equiv=2, uniform=3 };
    static const std::string strategyNames []
          = { "none", "maximum", "equidistribution", "uniform" };
    Strategy strategy = (Strategy) Dune::Fem::Parameter :: getEnum("adaptation.strategy", strategyNames );

    double localTol2 = 0;
    switch( strategy )
    {
    case maximum:
      {
        double maxError = 0 ;
        // sum up local estimators
        const typename ErrorIndicatorType :: const_iterator endEstimator = indicator_.end();
        for( typename ErrorIndicatorType :: const_iterator it = indicator_.begin();
              it != endEstimator; ++it )
        {
          maxError = std::max( *it, maxError );
        }
        // get global maxError
        maxError = grid_.comm().max( maxError );
        localTol2 = 0.25 * maxError ;
      }

      break ;
    case equiv:
      {
        double sumError = 0 ;
        // sum up local estimators
        const size_t indSize = indicator_.size();
        for( size_t i=0; i<indSize; ++i)
        {
          sumError += indicator_[ i ];
        }

        // get global sum of number of elements and local error sum
        double buffer[ 2 ] = { (double)indexSet_.size( 0 ) , sumError };
        grid_.comm().sum( buffer, 2 );
        localTol2 = 0.95 * tolerance*tolerance * buffer[ 1 ] / buffer[ 0 ];
      }
      break ;
    default:
      break ;
    }

    //! [Mark entities]
    int marked = 0;
    // loop over all elements
    const IteratorType end = dfSpace_.end();
    for( IteratorType it = dfSpace_.begin(); it != end; ++it )
    {
      const ElementType &entity = *it;
      // check local error indicator
      if( indicator_[ indexSet_.index( entity ) ] > localTol2 )
      {
        // mark entity for refinement
        grid_.mark( 1, entity );
        // grid was marked
        marked = 1;
      }
    }

    // get global max
    marked = grid_.comm().max( marked );
    //! [Mark entities]
    return bool(marked);
  }

protected:
  void clear ()
  {
    // resize and clear
    indicator_.resize( indexSet_.size( 0 ) );
    typedef typename ErrorIndicatorType :: iterator IteratorType;
    const IteratorType end = indicator_.end();
    for( IteratorType it = indicator_.begin(); it != end; ++it )
      *it = 0.0;
  }

  //! caclulate error on element
  template< class RHSFunctionType >
  void estimateLocal ( const RHSFunctionType &rhs, const ElementType &entity )
  {
    const typename ElementType :: Geometry &geometry = entity.geometry();

    const double volume = geometry.volume();
    const double h2 = (dimension == 2 ? volume : std :: pow( volume, 2.0 / (double)dimension ));
    const int index = indexSet_.index( entity );
    const LocalFunctionType uLocal = uh_.localFunction( entity );
    typename RHSFunctionType::LocalFunctionType localRhs = rhs.localFunction( entity );

    ElementQuadratureType quad( entity, 2*(dfSpace_.order() + 1) );

    model_.init(entity);

    // compute element residual
    const int numQuadraturePoints = quad.nop();
    for( int qp = 0; qp < numQuadraturePoints; ++qp )
    {
      const typename ElementQuadratureType :: CoordinateType &x = quad.point( qp );
      const double weight = quad.weight( qp ) * geometry.integrationElement( x );

      typename LocalFunctionType::RangeType values;
      uLocal.evaluate(x, values);

      typename LocalFunctionType::JacobianRangeType jacobian;
      uLocal.jacobian(x, jacobian);

      typename LocalFunctionType::HessianRangeType hessian;
      uLocal.hessian(x, hessian);

      RangeType y;
      localRhs.evaluate( quad[qp], y );

      RangeType tmp;
      model_.fluxDivergence(quad[qp], values, jacobian, hessian, tmp);
      //laplacianLocal( uLocal, quad[ qp ], tmp );
      y -= tmp;

      indicator_[ index ] += h2 * weight * y.two_norm2();
    }

    // calculate face contribution
    IntersectionIteratorType end = gridPart_.iend( entity );
    for( IntersectionIteratorType it = gridPart_.ibegin( entity ); it != end; ++it )
    {
      const IntersectionType &intersection = *it;
      // if we got an element neighbor
      if( intersection.neighbor() )
        estimateIntersection( intersection, entity, uLocal );
    }
  }

  //! caclulate error on boundary intersections
  void estimateIntersection ( const IntersectionType &intersection,
                              const ElementType &inside,
                              const LocalFunctionType &uInside )
  {
    const ElementType outside = Dune::Fem::make_entity( intersection.outside() );

    const int insideIndex = indexSet_.index( inside );
    const int outsideIndex = indexSet_.index( outside );

    // only compute intersection estimator once - either if the
    // outside entity has a larger index than the inside entity or
    // the intersection is on the boundary (parallel case...)
    const bool isOutsideInterior = (outside.partitionType() == Dune::InteriorEntity);
    if( !isOutsideInterior || (insideIndex < outsideIndex) )
    {
      const LocalFunctionType uOutside = uh_.localFunction( outside );

      double error;

      if( intersection.conforming() )
        error = estimateIntersection< true >( intersection, uInside, uOutside );
      else
        error = estimateIntersection< false >( intersection, uInside, uOutside );

      if( error > 0.0 )
      {
        const double volume
          = 0.5 * (inside.geometry().volume() + outside.geometry().volume());
        const double h = (dimension == 1 ? volume : std::pow( volume, 1.0 / (double)dimension ));
        indicator_[ insideIndex ] += h * error;
        if( isOutsideInterior )
          indicator_[ outsideIndex ] += h * error;
      }
    }
  }

  //! caclulate error on element intersections
  template< bool conforming  >
  double estimateIntersection ( const IntersectionType &intersection,
                                const LocalFunctionType &uInside,
                                const LocalFunctionType &uOutside ) const
  {
    // make sure correct method is called
    assert( intersection.conforming() == conforming );

    // use IntersectionQuadrature to create appropriate face quadratures
    typedef Dune :: Fem :: IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType;
    typedef typename IntersectionQuadratureType :: FaceQuadratureType Quadrature ;

    const int quadOrder = 2 * (dfSpace_.order() - 1);
    // create intersection quadrature
    IntersectionQuadratureType intersectionQuad( gridPart_, intersection, quadOrder );

    // get appropriate quadrature references
    const Quadrature &quadInside  = intersectionQuad.inside();
    const Quadrature &quadOutside = intersectionQuad.outside();

    double error = 0.0;
    const int numQuadraturePoints = quadInside.nop();
    for( int qp = 0; qp < numQuadraturePoints; ++qp )
    {
      DomainType integrationNormal
        = intersection.integrationOuterNormal( quadInside.localPoint( qp ) );
      const double integrationElement = integrationNormal.two_norm();

      // evaluate | (d u_l * n_l) + (d u_r * n_r) | = | (d u_l - d u_r) * n_l |
      RangeType valueInside, valueOutside;
      uInside.evaluate( quadInside[ qp ], valueInside );
      uOutside.evaluate( quadOutside[ qp ], valueOutside );
      JacobianRangeType jacobianInside, jacobianOutside;
      uInside.jacobian( quadInside[ qp ], jacobianInside );
      uOutside.jacobian( quadOutside[ qp ], jacobianOutside );
      JacobianRangeType fluxInside, fluxOutside;
      model_.diffusiveFlux(quadInside[qp], valueInside, jacobianInside, fluxInside);
      model_.diffusiveFlux(quadOutside[qp],valueOutside,jacobianOutside,fluxOutside);

      // apply normal
      RangeType jump;
      fluxInside -= fluxOutside;
      fluxInside.mv( integrationNormal, jump );

      error += quadInside.weight( qp ) * jump.two_norm2() / integrationElement;
    }
    return error;
  }

  template< class PointType >
  void laplacianLocal ( const LocalFunctionType &u_h,
                        const PointType &x, RangeType &result ) const
  {
    typename LocalFunctionType::HessianRangeType hessian;
    u_h.hessian( x, hessian );

    result = 0;
    for( int r = 0; r < LocalFunctionType::dimRange; ++r )
    {
      for( int i = 0; i < dimension; ++i )
        result[ r ] += hessian[ r ][ i ][ i ];
    }
  }
};

#endif // #ifndef ESTIMATOR_HH
