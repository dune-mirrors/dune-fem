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
#ifndef DGRHS_HH
#define DGRHS_HH

#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem/schemes/rhs.hh>


// assembleRHS
// -----------

template< class Model, class Function, class Neuman, class DiscreteFunction >
// template< class Model, class Function, class Neuman, class Dirichlet, class DiscreteFunction >
void assembleDGRHS ( const Model &model, const Function &function, const Neuman &neuman, // const Dirichlet &dirichlet,
                     DiscreteFunction &rhs,
                     double penalty )
{
  assembleRHS( model, function, neuman, rhs );
  if ( ! model.hasDirichletBoundary() )
    return;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  static const int dimDomain = LocalFunctionType::dimDomain;
  static const int dimRange = LocalFunctionType::dimRange;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;

  const DiscreteFunctionSpaceType &dfSpace = rhs.space();
  const int quadOrder = 2*dfSpace.order()+1;

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    if ( !entity.hasBoundaryIntersections() )
      continue;

    const GeometryType &geometry = entity.geometry();
    double area = geometry.volume();

    LocalFunctionType rhsLocal = rhs.localFunction( entity );

    bool needsCalculation = model.init( entity );
    if (! needsCalculation )
      continue;

    const IntersectionIteratorType iitend = dfSpace.gridPart().iend( entity );
    for( IntersectionIteratorType iit = dfSpace.gridPart().ibegin( entity ); iit != iitend; ++iit ) // looping over intersections
    {
      const IntersectionType &intersection = *iit;
      if ( ! intersection.boundary() ) // i.e. if intersection is on boundary: nothing to be done for Neumann zero b.c.
        continue;                      // since [u] = 0  and grad u.n = 0
      Dune::FieldVector<int,dimRange> components(0);
      if ( ! model.isDirichletIntersection( intersection, components) )
        continue;
#if 0
      const typename Dirichlet::LocalFunctionType dirichletLocal = dirichlet.localFunction( entity);

      typedef typename IntersectionType::Geometry  IntersectionGeometryType;
      const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

      const double intersectionArea = intersectionGeometry.volume();
      const double beta = penalty * intersectionArea / area;

      FaceQuadratureType quadInside( dfSpace.gridPart(), intersection, quadOrder, FaceQuadratureType::INSIDE );
      const size_t numQuadraturePoints = quadInside.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
        const DomainType normal = intersection.integrationOuterNormal( x );

        const double weight = quadInside.weight( pt );

        RangeType value;
        JacobianRangeType dvalue,advalue;

        RangeType vuOut;
        dirichletLocal.evaluate( quadInside[pt], vuOut );
        for (int r=0;r<dimRange;++r)
          if (!components[r]) // do not use dirichlet constraints here
            vuOut[r] = 0;

        value = vuOut;
        value *= beta * intersectionGeometry.integrationElement( x );

        //  [ u ] * { grad phi_en } = -normal(u+ - u-) * 0.5 grad phi_en
        // here we need a diadic product of u x n
        for (int r=0;r<dimRange;++r)
          for (int d=0;d<dimDomain;++d)
            dvalue[r][d] = -0.5 * normal[d] * vuOut[r];

        model.diffusiveFlux( quadInside[ pt ], vuOut, dvalue, advalue );

        value *= weight;
        advalue *= weight;
        rhsLocal.axpy( quadInside[ pt ], value, advalue );
      }
#endif
    }
  }
  rhs.communicate();
}

#endif // #ifndef DGRHS_HH
