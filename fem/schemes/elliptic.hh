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
#ifndef DUNE_FEM_SCHEMES_ELLIPTIC_HH
#define DUNE_FEM_SCHEMES_ELLIPTIC_HH

#include <cstddef>

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/common/bindguard.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/function/localfunction/const.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/fem/schemes/dirichletconstraints.hh>

// EllipticOperator
// ----------------

struct NoConstraints
{
  template <class ModelType, class DiscreteFunctionSpaceType>
  NoConstraints( const ModelType&, const DiscreteFunctionSpaceType& )
  {}

  template < class DiscreteFunctionType >
  void operator ()( const DiscreteFunctionType& u, DiscreteFunctionType& w ) const
  {}

  template < class GridFunctionType, class DiscreteFunctionType >
  void operator ()( const GridFunctionType& u, DiscreteFunctionType& w ) const
  {}

  template <class LinearOperator>
  void applyToOperator( LinearOperator& linearOperator ) const
  {}
};



//! [Class for elliptic operator]
template< class DomainDiscreteFunction, class RangeDiscreteFunction, class Model,
          class Constraints = Dune::DirichletConstraints< Model, typename RangeDiscreteFunction::DiscreteFunctionSpaceType > >
struct EllipticOperator
: public virtual Dune::Fem::Operator< DomainDiscreteFunction, RangeDiscreteFunction >
//! [Class for elliptic operator]
{
  typedef DomainDiscreteFunction DomainDiscreteFunctionType;
  typedef RangeDiscreteFunction  RangeDiscreteFunctionType;
  typedef Model                  ModelType;
  typedef Constraints      ConstraintsType;           // the class taking care of boundary constraints e.g. dirichlet bc

protected:
  typedef typename DomainDiscreteFunctionType::DiscreteFunctionSpaceType DomainDiscreteFunctionSpaceType;
  typedef typename DomainDiscreteFunctionType::LocalFunctionType         DomainLocalFunctionType;
  typedef typename DomainLocalFunctionType::RangeType                    DomainRangeType;
  typedef typename DomainLocalFunctionType::JacobianRangeType            DomainJacobianRangeType;
  typedef typename RangeDiscreteFunctionType::DiscreteFunctionSpaceType RangeDiscreteFunctionSpaceType;
  typedef typename RangeDiscreteFunctionType::LocalFunctionType         RangeLocalFunctionType;
  typedef typename RangeLocalFunctionType::RangeType                    RangeRangeType;
  typedef typename RangeLocalFunctionType::JacobianRangeType            RangeJacobianRangeType;

  // the following types must be identical for domain and range
  typedef typename RangeDiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;
  typedef typename RangeDiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename RangeDiscreteFunctionSpaceType::GridPartType  GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 1 > FaceQuadratureType;

public:
  EllipticOperator ( const ModelType &model, const RangeDiscreteFunctionSpaceType &rangeSpace )
    : model_( model ) , constraints_( model, rangeSpace )
  {}

  // prepare the solution vector
  void prepare( RangeDiscreteFunctionType &u ) const
  {
    // set boundary values for solution
    constraints()( u );
  }
  void prepare( const RangeDiscreteFunctionType &u, RangeDiscreteFunctionType &w ) const
  {
    // set boundary values for solution
    constraints()( u, w );
  }

  //! application operator
  virtual void operator() ( const DomainDiscreteFunctionType &u, RangeDiscreteFunctionType &w ) const
  { apply(u,w); }
  template <class GF>
  void apply( const GF &u, RangeDiscreteFunctionType &w ) const;

protected:
  const ModelType &model () const { return model_; }
  const ConstraintsType &constraints () const { return constraints_; }

private:
  const ModelType &model_;
  ConstraintsType constraints_;
};



// DifferentiableEllipticOperator
// ------------------------------

//! [Class for linearizable elliptic operator]
template< class JacobianOperator, class Model,
          class Constraints = Dune::DirichletConstraints< Model, typename JacobianOperator::RangeFunctionType::DiscreteFunctionSpaceType > >
struct DifferentiableEllipticOperator
: public EllipticOperator< typename JacobianOperator::DomainFunctionType, typename JacobianOperator::RangeFunctionType, Model, Constraints >,
  public Dune::Fem::DifferentiableOperator< JacobianOperator >
//! [Class for linearizable elliptic operator]
{
  typedef EllipticOperator< typename JacobianOperator::DomainFunctionType, typename JacobianOperator::RangeFunctionType, Model, Constraints > BaseType;

  typedef JacobianOperator JacobianOperatorType;

  typedef typename BaseType::DomainDiscreteFunctionType DomainDiscreteFunctionType;
  typedef typename BaseType::RangeDiscreteFunctionType  RangeDiscreteFunctionType;
  typedef typename BaseType::ModelType ModelType;

protected:
  typedef typename DomainDiscreteFunctionType::DiscreteFunctionSpaceType DomainDiscreteFunctionSpaceType;
  typedef typename DomainDiscreteFunctionType::LocalFunctionType         DomainLocalFunctionType;
  typedef typename DomainLocalFunctionType::RangeType                    DomainRangeType;
  typedef typename DomainLocalFunctionType::JacobianRangeType            DomainJacobianRangeType;
  typedef typename RangeDiscreteFunctionType::DiscreteFunctionSpaceType RangeDiscreteFunctionSpaceType;
  typedef typename RangeDiscreteFunctionType::LocalFunctionType         RangeLocalFunctionType;
  typedef typename RangeLocalFunctionType::RangeType                    RangeRangeType;
  typedef typename RangeLocalFunctionType::JacobianRangeType            RangeJacobianRangeType;

  // the following types must be identical for domain and range
  typedef typename RangeDiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;
  typedef typename RangeDiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename RangeDiscreteFunctionSpaceType::GridPartType  GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  typedef typename BaseType::QuadratureType QuadratureType;
  // quadrature for faces - used for Neuman b.c.
  typedef typename BaseType::FaceQuadratureType FaceQuadratureType;

public:
  //! contructor
  DifferentiableEllipticOperator ( const ModelType &model, const RangeDiscreteFunctionSpaceType &space )
  : BaseType( model, space )
  {}

  //! method to setup the jacobian of the operator for storage in a matrix
  void jacobian ( const DomainDiscreteFunctionType &u, JacobianOperatorType &jOp ) const
  { apply(u,jOp); }
  template <class GridFunctionType>
  void apply ( const GridFunctionType &u, JacobianOperatorType &jOp ) const;
  using BaseType::apply;

protected:
  using BaseType::model;
  using BaseType::constraints;
};



// Implementation of EllipticOperator
// ----------------------------------

template< class DomainDiscreteFunction, class RangeDiscreteFunction, class Model, class Constraints >
template<class GF>
void EllipticOperator< DomainDiscreteFunction, RangeDiscreteFunction, Model, Constraints >
  ::apply( const GF &u, RangeDiscreteFunctionType &w ) const
{
  w.clear();
  // get discrete function space
  const RangeDiscreteFunctionSpaceType &dfSpace = w.space();

  // get local representation of the discrete functions
  Dune::Fem::ConstLocalFunction< GF > uLocal( u );
  Dune::Fem::AddLocalContribution< RangeDiscreteFunctionType > wLocal( w );

  // iterate over grid
  for( const EntityType &entity : dfSpace )
  {
    if( !model().init( entity ) )
      continue;

    // get elements geometry
    const GeometryType &geometry = entity.geometry();

    auto uGuard = Dune::Fem::bindGuard( uLocal, entity );
    auto wGuard = Dune::Fem::bindGuard( wLocal, entity );

    // obtain quadrature order
    const int quadOrder = uLocal.order() + wLocal.order();

    { // element integral
      QuadratureType quadrature( entity, quadOrder );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        //! [Compute local contribution of operator]
        const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
        const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

        DomainRangeType vu;
        uLocal.evaluate( quadrature[ pt ], vu );
        DomainJacobianRangeType du;
        uLocal.jacobian( quadrature[ pt ], du );

        // compute mass contribution (studying linear case so linearizing around zero)
        RangeRangeType avu( 0 );
        model().source( quadrature[ pt ], vu, du, avu );
        avu *= weight;
        // add to local functional wLocal.axpy( quadrature[ pt ], avu );

        RangeJacobianRangeType adu( 0 );
        // apply diffusive flux
        model().diffusiveFlux( quadrature[ pt ], vu, du, adu );
        adu *= weight;

        // add to local function
        wLocal.axpy( quadrature[ pt ], avu, adu );
        //! [Compute local contribution of operator]
      }
    }

    if( model().hasNeumanBoundary() && entity.hasBoundaryIntersections() )
    {
      const IntersectionIteratorType iitend = dfSpace.gridPart().iend( entity );
      for( IntersectionIteratorType iit = dfSpace.gridPart().ibegin( entity ); iit != iitend; ++iit )
      {
        const IntersectionType &intersection = *iit;
        if( !intersection.boundary() )
          continue;

        Dune::FieldVector<int,RangeRangeType::dimension> components(0);

        const bool hasDirichletComponent = model().isDirichletIntersection( intersection, components );

        const auto &intersectionGeometry = intersection.geometry();
        FaceQuadratureType quadInside( dfSpace.gridPart(), intersection, quadOrder, FaceQuadratureType::INSIDE );
        const std::size_t numQuadraturePoints = quadInside.nop();
        for( std::size_t pt = 0; pt < numQuadraturePoints; ++pt )
        {
          const auto &x = quadInside.localPoint( pt );
          double weight = quadInside.weight( pt ) * intersectionGeometry.integrationElement( x );
          DomainRangeType vu;
          uLocal.evaluate( quadInside[ pt ], vu );
          RangeRangeType alpha( 0 );
          model().alpha( quadInside[ pt ], vu, alpha );
          alpha *= weight;
          for( int k = 0; k < RangeRangeType::dimension; ++k )
            if( hasDirichletComponent && components[ k ] )
              alpha[ k ] = 0;
          wLocal.axpy( quadInside[ pt ], alpha );
        }
      }
    }
  }

  w.communicate();
  // apply constraints, e.g. Dirichlet contraints, to the result
  constraints()( u, w );
}



// Implementation of DifferentiableEllipticOperator
// ------------------------------------------------

template< class JacobianOperator, class Model, class Constraints >
template<class GF>
void DifferentiableEllipticOperator< JacobianOperator, Model, Constraints >
  ::apply ( const GF &u, JacobianOperator &jOp ) const
{
  typedef typename JacobianOperator::LocalMatrixType LocalMatrixType;
  typedef typename DomainDiscreteFunctionSpaceType::BasisFunctionSetType DomainBasisFunctionSetType;
  typedef typename RangeDiscreteFunctionSpaceType::BasisFunctionSetType  RangeBasisFunctionSetType;

  const DomainDiscreteFunctionSpaceType &domainSpace = jOp.domainSpace();
  const RangeDiscreteFunctionSpaceType  &rangeSpace = jOp.rangeSpace();

  Dune::Fem::DiagonalStencil<DomainDiscreteFunctionSpaceType,RangeDiscreteFunctionSpaceType> stencil( domainSpace, rangeSpace );
  jOp.reserve(stencil);
  jOp.clear();

  const int domainBlockSize = domainSpace.localBlockSize; // is equal to 1 for scalar functions
  std::vector< typename DomainLocalFunctionType::RangeType >         phi( domainSpace.blockMapper().maxNumDofs()*domainBlockSize );
  std::vector< typename DomainLocalFunctionType::JacobianRangeType > dphi( domainSpace.blockMapper().maxNumDofs()*domainBlockSize );
  const int rangeBlockSize = rangeSpace.localBlockSize; // is equal to 1 for scalar functions
  std::vector< typename RangeLocalFunctionType::RangeType >         rphi( rangeSpace.blockMapper().maxNumDofs()*rangeBlockSize );
  std::vector< typename RangeLocalFunctionType::JacobianRangeType > rdphi( rangeSpace.blockMapper().maxNumDofs()*rangeBlockSize );

  Dune::Fem::ConstLocalFunction< GF > uLocal( u );
  for( const EntityType &entity : rangeSpace )
  {
    if( !model().init( entity ) )
      continue;

    const GeometryType &geometry = entity.geometry();

    auto uGuard = Dune::Fem::bindGuard( uLocal, entity );
    LocalMatrixType jLocal = jOp.localMatrix( entity, entity );

    const DomainBasisFunctionSetType &domainBaseSet = jLocal.domainBasisFunctionSet();
    const RangeBasisFunctionSetType &rangeBaseSet  = jLocal.rangeBasisFunctionSet();
    const std::size_t domainNumBasisFunctions = domainBaseSet.size();

    QuadratureType quadrature( entity, domainSpace.order()+rangeSpace.order() );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      //! [Assembling the local matrix]
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
      const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

      // evaluate all basis functions at given quadrature point
      domainBaseSet.evaluateAll( quadrature[ pt ], phi );
      rangeBaseSet.evaluateAll( quadrature[ pt ], rphi );

      // evaluate jacobians of all basis functions at given quadrature point
      domainBaseSet.jacobianAll( quadrature[ pt ], dphi );
      rangeBaseSet.jacobianAll( quadrature[ pt ], rdphi );

      // get value for linearization
      DomainRangeType u0;
      DomainJacobianRangeType jacU0;
      uLocal.evaluate( quadrature[ pt ], u0 );
      uLocal.jacobian( quadrature[ pt ], jacU0 );

      RangeRangeType aphi( 0 );
      RangeJacobianRangeType adphi( 0 );
      for( std::size_t localCol = 0; localCol < domainNumBasisFunctions; ++localCol )
      {
        // if mass terms or right hand side is present
        model().linSource( u0, jacU0, quadrature[ pt ], phi[ localCol ], dphi[localCol], aphi );

        // if gradient term is present
        model().linDiffusiveFlux( u0, jacU0, quadrature[ pt ], phi[ localCol ], dphi[ localCol ], adphi );

        // get column object and call axpy method
        jLocal.column( localCol ).axpy( rphi, rdphi, aphi, adphi, weight );
      }
      //! [Assembling the local matrix]
    }

    if( model().hasNeumanBoundary() && entity.hasBoundaryIntersections() )
    {
      const IntersectionIteratorType iitend = rangeSpace.gridPart().iend( entity );
      for( IntersectionIteratorType iit = rangeSpace.gridPart().ibegin( entity ); iit != iitend; ++iit ) // looping over intersections
      {
        const IntersectionType &intersection = *iit;
        if( !intersection.boundary() )
          continue;

        Dune::FieldVector<int,RangeRangeType::dimension> components(0);
        bool hasDirichletComponent = model().isDirichletIntersection( intersection, components );

        const auto &intersectionGeometry = intersection.geometry();
        FaceQuadratureType quadInside( rangeSpace.gridPart(), intersection, domainSpace.order()+rangeSpace.order(), FaceQuadratureType::INSIDE );
        const std::size_t numQuadraturePoints = quadInside.nop();
        for( std::size_t pt = 0; pt < numQuadraturePoints; ++pt )
        {
          const auto &x = quadInside.localPoint( pt );
          const double weight = quadInside.weight( pt ) * intersectionGeometry.integrationElement( x );
          DomainRangeType u0;
          uLocal.evaluate( quadInside[ pt ], u0 );
          domainBaseSet.evaluateAll( quadInside[ pt ], phi );
          for( std::size_t localCol = 0; localCol < domainNumBasisFunctions; ++localCol )
          {
            RangeRangeType alpha( 0 );
            model().linAlpha( u0,quadInside[ pt ], phi[ localCol ], alpha );
            for( int k = 0; k < RangeRangeType::dimension; ++k )
              if( hasDirichletComponent && components[ k ] )
                alpha[ k ] = 0;
            jLocal.column( localCol ).axpy( phi, alpha, weight );
          }
        }
      }
    }
  } // end grid traversal

  // apply constraints to matrix operator
  constraints().applyToOperator( jOp );
  jOp.communicate();
}

#endif // #ifndef DUNE_FEM_SCHEMES_ELLIPTIC_HH
