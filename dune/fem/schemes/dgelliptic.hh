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
#ifndef DGELLIPTIC_HH
#define DGELLIPTIC_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>

#include <dune/fem/schemes/elliptic.hh>

// EllipticOperator
// ----------------

template <class DFSpace>
struct DefaultPenalty
{
  DefaultPenalty(const DFSpace &space, double penalty)
  : space_(space)
  , penalty_(penalty)
  {}
  template <class Intersection>
  double operator()(const Intersection &intersection,
                    double intersectionArea, double area, double nbArea) const
  {
    const double hInv = intersectionArea / std::min( area, nbArea );
    return penalty_ * hInv;
  }
  const double &factor() const { return penalty_; }
  private:
  const DFSpace &space_;
  double penalty_;
};

template< class DiscreteFunction, class Model,
  class Penalty = DefaultPenalty<typename DiscreteFunction::DiscreteFunctionSpaceType > >
struct DGEllipticOperator
: public virtual Dune::Fem::Operator< DiscreteFunction >
{
  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model            ModelType;

  typedef DiscreteFunction DomainDiscreteFunctionType;
  typedef DiscreteFunction RangeDiscreteFunctionType;
protected:
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;

  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  static const int dimDomain = LocalFunctionType::dimDomain;
  static const int dimRange = LocalFunctionType::dimRange;

public:
  //! constructor
  DGEllipticOperator ( const DiscreteFunctionSpaceType &space,
                       const ModelType &model,
                       const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
  : model_( model ),
    penalty_( space, parameter.getValue< double >( "penalty", 40 ) )
  {
    // std::cout << "dg operator with penalty:" << penalty_.factor() << std::endl;
  }

  DGEllipticOperator ( const DiscreteFunctionSpaceType &dSpace,
                       const DiscreteFunctionSpaceType &rSpace,
                       const ModelType &model,
                       const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
  : model_( model ),
    penalty_( dSpace, parameter.getValue< double >( "penalty", 40 ) )
  {
    // std::cout << "dg operator with penalty:" << penalty_.factor() << std::endl;
  }

  //! application operator
  virtual void operator() ( const DomainDiscreteFunctionType &u, RangeDiscreteFunctionType &w ) const
  { apply(u,w); }
  template <class GF>
  void apply( const GF &u, RangeDiscreteFunctionType &w ) const;

protected:
  const ModelType &model () const { return model_; }
  Penalty penalty() const { return penalty_; }

private:
  const ModelType &model_;
  Penalty penalty_;
};

// DifferentiableDGEllipticOperator
// ------------------------------

template< class JacobianOperator, class Model,
  class Penalty = DefaultPenalty<typename JacobianOperator::DomainFunctionType::DiscreteFunctionSpaceType > >
struct DifferentiableDGEllipticOperator
: public DGEllipticOperator< typename JacobianOperator::DomainFunctionType, Model, Penalty >,
  public Dune::Fem::DifferentiableOperator< JacobianOperator >
{
  typedef DGEllipticOperator< typename JacobianOperator::DomainFunctionType, Model, Penalty > BaseType;

  typedef JacobianOperator JacobianOperatorType;

  typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename BaseType::ModelType ModelType;
  typedef typename BaseType::DomainDiscreteFunctionType DomainDiscreteFunctionType;
  typedef typename BaseType::RangeDiscreteFunctionType RangeDiscreteFunctionType;

protected:
  typedef typename DomainDiscreteFunctionType::DiscreteFunctionSpaceType DomainDiscreteFunctionSpaceType;
  typedef typename DomainDiscreteFunctionType::LocalFunctionType         DomainLocalFunctionType;
  typedef typename DomainLocalFunctionType::RangeType                    DomainRangeType;
  typedef typename DomainLocalFunctionType::JacobianRangeType            DomainJacobianRangeType;
  typedef typename RangeDiscreteFunctionType::DiscreteFunctionSpaceType RangeDiscreteFunctionSpaceType;
  typedef typename RangeDiscreteFunctionType::LocalFunctionType         RangeLocalFunctionType;
  typedef typename RangeLocalFunctionType::RangeType                    RangeRangeType;
  typedef typename RangeLocalFunctionType::JacobianRangeType            RangeJacobianRangeType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::RangeFieldType RangeFieldType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;

  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  static const int dimDomain = LocalFunctionType::dimDomain;
  static const int dimRange = LocalFunctionType::dimRange;

public:
  //! constructor
  DifferentiableDGEllipticOperator ( const DiscreteFunctionSpaceType &space,
                         const ModelType &model,
                         const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
  : BaseType( space, model, parameter )
  {}
  DifferentiableDGEllipticOperator ( const DiscreteFunctionSpaceType &dSpace,
                         const DiscreteFunctionSpaceType &rSpace,
                         const ModelType &model,
                         const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
  : BaseType( dSpace, rSpace, model, parameter )
  {}

  //! method to setup the jacobian of the operator for storage in a matrix
  void jacobian ( const DomainDiscreteFunctionType &u, JacobianOperatorType &jOp ) const
  { apply(u,jOp); }
  template <class GridFunctionType>
  void apply ( const GridFunctionType &u, JacobianOperatorType &jOp ) const;
  using BaseType::apply;

protected:
  using BaseType::model;
  using BaseType::penalty;
};

// Implementation of DGEllipticOperator
// ----------------------------------

template< class RangeDiscreteFunction, class Model, class Penalty >
template<class GF>
void DGEllipticOperator< RangeDiscreteFunction, Model, Penalty >
  ::apply( const GF &u, RangeDiscreteFunctionType &w ) const
{
  // clear destination
  w.clear();

  // get discrete function space
  const DiscreteFunctionSpaceType &dfSpace = w.space();
  Dune::Fem::ConstLocalFunction< GF > uLocal( u );
  Dune::Fem::ConstLocalFunction< GF > uOutLocal( u );
  Dune::Fem::AddLocalContribution< RangeDiscreteFunctionType > wLocal( w );

  // iterate over grid
  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    // get entity (here element)
    const EntityType &entity = *it;

    bool needsCalculation = model().init( entity );
    if (! needsCalculation )
      continue;

    // get elements geometry
    const GeometryType &geometry = entity.geometry();

    // get local representation of the discrete functions
    auto uGuard = Dune::Fem::bindGuard( uLocal, entity );
    auto wGuard = Dune::Fem::bindGuard( wLocal, entity );

    // obtain quadrature order
    const int quadOrder = uLocal.order() + wLocal.order();

    { // element integral
      QuadratureType quadrature( entity, quadOrder );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
        const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

        RangeType vu;
        uLocal.evaluate( quadrature[ pt ], vu );
        JacobianRangeType du;
        uLocal.jacobian( quadrature[ pt ], du );

        // compute mass contribution (studying linear case so linearizing around zero)
        RangeType avu( 0 );
        model().source( quadrature[ pt ], vu, du, avu );
        avu *= weight;
        // add to local functional wLocal.axpy( quadrature[ pt ], avu );

        JacobianRangeType adu( 0 );
        // apply diffusive flux
        model().flux( quadrature[ pt ], vu, du, adu );
        adu *= weight;

        // add to local function
        wLocal.axpy( quadrature[ pt ], avu, adu );
      }
    }
    if ( ! dfSpace.continuous() )
    {
      const double area = entity.geometry().volume();
      const IntersectionIteratorType iitend = dfSpace.gridPart().iend( entity );
      for( IntersectionIteratorType iit = dfSpace.gridPart().ibegin( entity ); iit != iitend; ++iit ) // looping over intersections
      {
        //! [Compute skeleton terms: iterate over intersections]
        const IntersectionType &intersection = *iit;
        if ( intersection.neighbor() )
        {
          const EntityType outside = intersection.outside() ;

          typedef typename IntersectionType::Geometry  IntersectionGeometryType;
          const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

          // compute penalty factor
          const double intersectionArea = intersectionGeometry.volume();
          const double beta = penalty()(intersection,intersectionArea,area,outside.geometry().volume());

          auto uOutGuard = Dune::Fem::bindGuard( uOutLocal, outside );

          FaceQuadratureType quadInside( dfSpace.gridPart(), intersection, quadOrder, FaceQuadratureType::INSIDE );
          FaceQuadratureType quadOutside( dfSpace.gridPart(), intersection, quadOrder, FaceQuadratureType::OUTSIDE );
          const size_t numQuadraturePoints = quadInside.nop();
          //! [Compute skeleton terms: iterate over intersections]

          for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
          {
            //! [Compute skeleton terms: obtain required values on the intersection]
            // get coordinate of quadrature point on the reference element of the intersection
            const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
            DomainType normal = intersection.integrationOuterNormal( x );
            double faceVol = normal.two_norm();
            normal /= faceVol; // make it into a unit normal
            const double weight = quadInside.weight( pt ) * faceVol;

            RangeType vuIn,vuOut;
            JacobianRangeType duIn, duOut;
            uLocal.evaluate( quadInside[ pt ], vuIn );
            uLocal.jacobian( quadInside[ pt ], duIn );
            uOutLocal.evaluate( quadOutside[ pt ], vuOut );
            uOutLocal.jacobian( quadOutside[ pt ], duOut );
            RangeType jump = vuIn - vuOut;
            // compute -0.5 * [u] x normal
            JacobianRangeType dvalue;
            for (int r=0;r<dimRange;++r)
              for (int d=0;d<dimDomain;++d)
                dvalue[r][d] = -0.5 * normal[d] * jump[r];
            JacobianRangeType aduIn,aduOut;
            model().init( outside );
            model().flux( quadOutside[ pt ], vuOut, duOut, aduOut );
            auto pFactor = model().penalty( quadOutside[ pt ], vuOut );
            model().init( entity );
            model().flux( quadInside[ pt ], vuIn, duIn, aduIn );
            JacobianRangeType affine;
            model().flux( quadInside[ pt ], RangeType(0), JacobianRangeType(0), affine);
            pFactor += model().penalty( quadInside[ pt ], vuIn );
            //! [Compute skeleton terms: obtain required values on the intersection]

            //! [Compute skeleton terms: compute factors for axpy method]
            RangeType value;
            JacobianRangeType advalue;
            // penalty term : beta [u] [phi] = beta (u+ - u-)(phi+ - phi-)=beta (u+ - u-)phi+
            value = jump;
            for (unsigned int r=0;r<dimRange;++r)
              value[r] *= beta * pFactor[r]/2.;
            // {A grad u}.[phi] = {A grad u}.phi+ n_+ = 0.5*(grad u+ + grad u-).n_+ phi+
            aduIn += aduOut;
            aduIn *= -0.5;
            aduIn.umv(normal,value);
            //  [ u ] * { {A} grad phi_en } = -normal(u+ - u-) * 0.5 {A}    grad phi_en
            //  we actually need  G(u)tau with G(x,u) = d/sigma_j D_i(x,u,sigma)
            //  - we might need to assume D(x,u,sigma) = G(x,u)sigma + affine(x)
            model().flux( quadInside[ pt ], vuIn, dvalue, advalue );
            // advalue -= affine;

            value *= weight;
            advalue *= weight;
            wLocal.axpy( quadInside[ pt ], value, advalue );
            //! [Compute skeleton terms: compute factors for axpy method]
          }
        }
        else if( intersection.boundary() )
        {
          std::array<int,dimRange> components;
          components.fill(0);
          model().isDirichletIntersection( intersection, components);

          typedef typename IntersectionType::Geometry IntersectionGeometryType;
          const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

          // compute penalty factor
          const double intersectionArea = intersectionGeometry.volume();
          const double beta = penalty()(intersection,intersectionArea,area,area);

          FaceQuadratureType quadInside( dfSpace.gridPart(), intersection, quadOrder, FaceQuadratureType::INSIDE );
          const size_t numQuadraturePoints = quadInside.nop();

          for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
          {
            const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
            const DomainType normal = intersection.integrationOuterNormal( x );
            const double weight = quadInside.weight( pt );

            RangeType bndValue;
            model().dirichlet(1, quadInside[pt], bndValue);

            RangeType value;
            JacobianRangeType dvalue,advalue;

            RangeType vuIn,jump;
            JacobianRangeType duIn, aduIn;
            uLocal.evaluate( quadInside[ pt ], vuIn );
            uLocal.jacobian( quadInside[ pt ], duIn );
            jump = vuIn;
            jump -= bndValue;

            // penalty term : beta [u] [phi] = beta (u+ - u-)(phi+ - phi-)=beta (u+ - u-)phi+
            auto pFactor = model().penalty( quadInside[ pt ], vuIn );
            value = jump;
            for (unsigned int r=0;r<dimRange;++r)
              value[r] *= beta * pFactor[r] * intersectionGeometry.integrationElement( x );

            //  [ u ] * { grad phi_en } = -normal(u+ - u-) * 0.5 grad phi_en
            // here we need a diadic product of u x n
            for (int r=0;r<dimRange;++r)
              for (int d=0;d<dimDomain;++d)
                dvalue[r][d] = -0.5 * normal[d] * jump[r];

            model().flux( quadInside[ pt ], jump, dvalue, advalue );

            // consistency term
            // {A grad u}.[phi] = {A grad u}.phi+ n_+ = 0.5*(grad u+ + grad u-).n_+ phi+
            model().flux( quadInside[ pt ], vuIn, duIn, aduIn );
            aduIn.mmv(normal,value);

            for (int r=0;r<dimRange;++r)
              if (!components[r]) // do not use dirichlet constraints here
              {
                value[r] = 0;
                advalue[r] = 0;
              }

            value *= weight;
            advalue *= weight;
            wLocal.axpy( quadInside[ pt ], value, advalue );
          }
        }
      }
    }

  }

  // communicate data (in parallel runs)
  w.communicate();
}

// Implementation of DifferentiableDGEllipticOperator
// ------------------------------------------------
template< class JacobianOperator, class Model, class Penalty >
template<class GF>
void DifferentiableDGEllipticOperator< JacobianOperator, Model, Penalty >
  ::apply ( const GF &u, JacobianOperator &jOp ) const
{
  // std::cout << "starting assembly\n";
  // Dune::Timer timer;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  typedef Dune::Fem::TemporaryLocalMatrix< DomainDiscreteFunctionSpaceType, RangeDiscreteFunctionSpaceType > TmpLocalMatrixType;

  const DomainDiscreteFunctionSpaceType &domainSpace = jOp.domainSpace();
  const RangeDiscreteFunctionSpaceType  &rangeSpace = jOp.rangeSpace();

  Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> stencil(domainSpace,rangeSpace);
  jOp.reserve(stencil);
  jOp.clear();

  Dune::Fem::ConstLocalFunction< GF > uLocal( u );
  TmpLocalMatrixType jLocal( domainSpace, rangeSpace );
  TmpLocalMatrixType localOpNb( domainSpace, rangeSpace );

  const GridPartType& gridPart = rangeSpace.gridPart();

  const unsigned int numDofs = rangeSpace.blockMapper().maxNumDofs() *
                               DiscreteFunctionSpaceType :: localBlockSize ;

  std::vector< RangeType > phi( numDofs );
  std::vector< JacobianRangeType > dphi( numDofs );

  std::vector< RangeType > phiNb( numDofs );
  std::vector< JacobianRangeType > dphiNb( numDofs );

  // std::cout << "   in assembly: start element loop size=" << rangeSpace.gridPart().grid().size(0) << " time=  " << timer.elapsed() << std::endl;;
  const IteratorType end = rangeSpace.end();
  for( IteratorType it = rangeSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;

    bool needsCalculation = model().init( entity );
    if (! needsCalculation )
      continue;

    const GeometryType geometry = entity.geometry();

    auto uGuard = Dune::Fem::bindGuard( uLocal, entity );

    jLocal.bind( entity, entity );
    jLocal.clear();

    const BasisFunctionSetType &baseSet = jLocal.domainBasisFunctionSet();
    const unsigned int numBaseFunctions = baseSet.size();

    QuadratureType quadrature( entity, 2*rangeSpace.order() );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
      const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

      // evaluate all basis functions at given quadrature point
      baseSet.evaluateAll( quadrature[ pt ], phi );

      // evaluate jacobians of all basis functions at given quadrature point
      baseSet.jacobianAll( quadrature[ pt ], dphi );

      // get value for linearization
      RangeType u0;
      JacobianRangeType jacU0;
      uLocal.evaluate( quadrature[ pt ], u0 );
      uLocal.jacobian( quadrature[ pt ], jacU0 );

      RangeType aphi( 0 );
      JacobianRangeType adphi( 0 );
      for( unsigned int localCol = 0; localCol < numBaseFunctions; ++localCol )
      {
        // if mass terms or right hand side is present
        model().linSource( u0, jacU0, quadrature[ pt ], phi[ localCol ], dphi[ localCol ], aphi );

        // if gradient term is present
        model().linFlux( u0, jacU0, quadrature[ pt ], phi[ localCol ], dphi[ localCol ], adphi );

        // get column object and call axpy method
        jLocal.column( localCol ).axpy( phi, dphi, aphi, adphi, weight );
      }
    }
    if ( rangeSpace.continuous() )
    {
      jOp.addLocalMatrix( entity, entity, jLocal );
      jLocal.unbind();
      continue;
    }

    double area = geometry.volume();
    const IntersectionIteratorType endiit = gridPart.iend( entity );
    for ( IntersectionIteratorType iit = gridPart.ibegin( entity );
          iit != endiit ; ++ iit )
    {
      const IntersectionType& intersection = *iit ;

      if( intersection.neighbor() )
      {
        const EntityType neighbor = intersection.outside() ;

        typedef typename IntersectionType::Geometry  IntersectionGeometryType;
        const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

        //! [Assemble skeleton terms: get contributions on off diagonal block]
        // get local matrix for face entries
        localOpNb.bind( neighbor, entity );
        localOpNb.clear();

        // get neighbor's base function set
        const BasisFunctionSetType &baseSetNb = localOpNb.domainBasisFunctionSet();
        //! [Assemble skeleton terms: get contributions on off diagonal block]

        // compute penalty factor
        const double intersectionArea = intersectionGeometry.volume();
        const double beta = penalty()(intersection,intersectionArea,area,neighbor.geometry().volume());

        // here we assume that the intersection is conforming
        FaceQuadratureType faceQuadInside(gridPart, intersection, 2*rangeSpace.order() + 1,
                                          FaceQuadratureType::INSIDE);
        FaceQuadratureType faceQuadOutside(gridPart, intersection, 2*rangeSpace.order() + 1,
                                           FaceQuadratureType::OUTSIDE);

        const size_t numFaceQuadPoints = faceQuadInside.nop();
        for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
        {
          const typename FaceQuadratureType::LocalCoordinateType &x = faceQuadInside.localPoint( pt );
          DomainType normal = intersection.integrationOuterNormal( x );
          double faceVol = normal.two_norm();
          normal /= faceVol; // make it into a unit normal

          const double quadWeight = faceQuadInside.weight( pt );
          const double weight = quadWeight * faceVol;

          //! [Assemble skeleton terms: obtain values om quadrature point]
          RangeType u0En;
          JacobianRangeType u0EnJac;
          uLocal.evaluate( faceQuadInside[ pt ], u0En );
          uLocal.jacobian( faceQuadInside[ pt ], u0EnJac );

          /////////////////////////////////////////////////////////////
          // evaluate basis function of face inside E^- (entity)
          /////////////////////////////////////////////////////////////

          // evaluate all basis functions for quadrature point pt
          baseSet.evaluateAll( faceQuadInside[ pt ], phi );

          // evaluate the jacobians of all basis functions
          baseSet.jacobianAll( faceQuadInside[ pt ], dphi );

          /////////////////////////////////////////////////////////////
          // evaluate basis function of face inside E^+ (neighbor)
          /////////////////////////////////////////////////////////////

          // evaluate all basis functions for quadrature point pt on neighbor
          baseSetNb.evaluateAll( faceQuadOutside[ pt ], phiNb );

          // evaluate the jacobians of all basis functions on neighbor
          baseSetNb.jacobianAll( faceQuadOutside[ pt ], dphiNb );

          model().init( entity );
          for( unsigned int i = 0; i < numBaseFunctions; ++i )
          {
            JacobianRangeType adphiEn = dphi[ i ];
            model().linFlux( u0En, u0EnJac,   faceQuadInside[ pt ], phi[i], adphiEn, dphi[ i ] );
          }
          auto pFactor = model().penalty( faceQuadInside[ pt ], u0En );
          model().init( neighbor );
          for( unsigned int i = 0; i < numBaseFunctions; ++i )
          {
            JacobianRangeType adphiNb = dphiNb[ i ];
            model().linFlux( u0En, u0EnJac,   faceQuadOutside[ pt ], phiNb[i], adphiNb, dphiNb[ i ] );
          }
          pFactor += model().penalty( faceQuadOutside[ pt ], u0En );
          model().init( entity );
          //! [Assemble skeleton terms: obtain values om quadrature point]

          //! [Assemble skeleton terms: compute factors for axpy method]
          for( unsigned int localCol = 0; localCol < numBaseFunctions; ++localCol )
          {
            RangeType valueEn(0), valueNb(0);
            JacobianRangeType dvalueEn(0), dvalueNb(0);

            //  -{ A grad u } * [ phi_en ]
            dphi[localCol].usmv( -0.5, normal, valueEn );

            //  -{ A grad u } * [ phi_en ]
            dphiNb[localCol].usmv( -0.5, normal, valueNb );

            //  [ u ] * [ phi_en ] = u^- * phi_en^-
            for (unsigned int r=0;r<dimRange;++r)
            {
              valueEn[r] += beta*pFactor[r]/2.*phi[localCol][r];
              valueNb[r] -= beta*pFactor[r]/2.*phiNb[localCol][r];
            }
            // here we need a diadic product of u x n
            for ( int r=0; r< dimRange; ++r )
              for ( int d=0; d< dimDomain; ++d )
              {
                //  [ u ] * { grad phi_en }
                dvalueEn[r][d] = - 0.5 * normal[d] * phi[localCol][r];

                //  [ u ] * { grad phi_en }
                dvalueNb[r][d] = 0.5 * normal[d] * phiNb[localCol][r];
              }

            jLocal.column( localCol ).axpy( phi, dphi, valueEn, dvalueEn, weight );
            localOpNb.column( localCol ).axpy( phi, dphi, valueNb, dvalueNb, weight );
          }
        }
        jOp.addLocalMatrix( neighbor, entity, localOpNb );
        localOpNb.unbind();
      }
      else if( intersection.boundary() )
      {
        std::array<int,dimRange> components;
        components.fill(0);
        model().isDirichletIntersection( intersection, components);

        typedef typename IntersectionType::Geometry  IntersectionGeometryType;
        const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

        // compute penalty factor
        const double intersectionArea = intersectionGeometry.volume();
        const double beta = penalty()(intersection,intersectionArea,area,area);

        // here we assume that the intersection is conforming
        FaceQuadratureType faceQuadInside(gridPart, intersection, 2*rangeSpace.order() + 1,
                                          FaceQuadratureType::INSIDE);

        const size_t numFaceQuadPoints = faceQuadInside.nop();
        for( size_t pt = 0; pt < numFaceQuadPoints; ++pt )
        {
          const typename FaceQuadratureType::LocalCoordinateType &x = faceQuadInside.localPoint( pt );
          DomainType normal = intersection.integrationOuterNormal( x );
          double faceVol = normal.two_norm();
          normal /= faceVol; // make it into a unit normal

          const double quadWeight = faceQuadInside.weight( pt );
          const double weight = quadWeight * faceVol;

          RangeType u0En;
          JacobianRangeType u0EnJac;
          uLocal.evaluate( faceQuadInside[ pt ], u0En );
          uLocal.jacobian( faceQuadInside[ pt ], u0EnJac );

          /////////////////////////////////////////////////////////////
          // evaluate basis function of face inside E^- (entity)
          /////////////////////////////////////////////////////////////

          // evaluate all basis functions for quadrature point pt
          baseSet.evaluateAll( faceQuadInside[ pt ], phi );

          // evaluate the jacobians of all basis functions
          baseSet.jacobianAll( faceQuadInside[ pt ], dphi );

          for( unsigned int i = 0; i < numBaseFunctions; ++i )
          {
            JacobianRangeType adphiEn = dphi[ i ];
            model().linFlux( u0En, u0EnJac, faceQuadInside[ pt ], phi[i], adphiEn, dphi[ i ] );
          }

          auto pFactor = model().penalty( faceQuadInside[ pt ], u0En );

          for( unsigned int localCol = 0; localCol < numBaseFunctions; ++localCol )
          {
            RangeType valueEn(0);
            JacobianRangeType dvalueEn(0);

            //  -{ A grad u } * [ phi_en ]
            dphi[localCol].usmv( -1.0, normal, valueEn );

            //  [ u ] * [ phi_en ] = u^- * phi_en^-
            for (unsigned int r=0;r<dimRange;++r)
              valueEn[r] += beta*pFactor[r]*phi[ localCol ][r];

            // here we need a diadic product of u x n
            for ( int r=0; r< dimRange; ++r )
              for ( int d=0; d< dimDomain; ++d )
              {
                //  [ u ] * { grad phi_en }
                dvalueEn[r][d] = - 0.5 * normal[d] * phi[localCol][r];
              }

            for (int r=0;r<dimRange;++r)
              if (!components[r]) // do not use dirichlet constraints here
              {
                valueEn[r] = 0;
                dvalueEn[r] = 0;
              }

            jLocal.column( localCol ).axpy( phi, dphi, valueEn, dvalueEn, weight );
          }
        }
      } // is boundary
    } // end of intersection loop
    jOp.addLocalMatrix( entity, entity, jLocal );
    jLocal.unbind();
  } // end grid traversal
  jOp.flushAssembly();
  // std::cout << "   in assembly: final    " << timer.elapsed() << std::endl;;
}

#endif // ELLIPTIC_HH
