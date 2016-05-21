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
#ifndef RHS_HH
#define RHS_HH

#include <dune/fem/quadrature/cachingquadrature.hh>


// assembleRHS
// -----------

template< class Model, class Function, class Neuman, class DiscreteFunction >
void assembleRHS ( const Model &model, const Function &function, const Neuman &neuman, DiscreteFunction &rhs )
{
  const bool neumanBnd = model.hasNeumanBoundary();
  rhs.clear();
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  const DiscreteFunctionSpaceType &dfSpace = rhs.space();

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
  const int quadOrder = 2*dfSpace.order()+1;

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    const typename Function::LocalFunctionType localFunction =
             function.localFunction( entity);
    LocalFunctionType rhsLocal = rhs.localFunction( entity );
    typedef typename Function::RangeType RangeType;

    QuadratureType quadrature( entity, quadOrder );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      // obtain quadrature point
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );

      // evaluate f
      RangeType f;
      localFunction.evaluate( quadrature[ pt ], f );

      // multiply by quadrature weight
      f *= quadrature.weight( pt ) * geometry.integrationElement( x );

      // add f * phi_i to rhsLocal[ i ]
      rhsLocal.axpy( quadrature[ pt ], f );
    }
    if (neumanBnd)
    {
      if ( !entity.hasBoundaryIntersections() )
        continue;

      const IntersectionIteratorType iitend = dfSpace.gridPart().iend( entity );
      for( IntersectionIteratorType iit = dfSpace.gridPart().ibegin( entity ); iit != iitend; ++iit ) // looping over intersections
      {
        const IntersectionType &intersection = *iit;
        if ( ! intersection.boundary() )
          continue;
        Dune::FieldVector<bool,RangeType::dimension> components(true);
        // if ( model.isDirichletIntersection( intersection, components) )
        //   continue;
        bool hasDirichletComponent = model.isDirichletIntersection( intersection, components);

        const typename Neuman::LocalFunctionType neumanLocal = neuman.localFunction( entity);

        const typename IntersectionType::Geometry &intersectionGeometry = intersection.geometry();
        FaceQuadratureType quadInside( dfSpace.gridPart(), intersection, quadOrder, FaceQuadratureType::INSIDE );
        const size_t numQuadraturePoints = quadInside.nop();
        for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
        {
          const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
          RangeType nval;
          neumanLocal.evaluate(quadInside[pt], nval);
          nval *= quadInside.weight( pt ) * intersectionGeometry.integrationElement( x );
          for(int k = 0; k < RangeType::dimension; ++k)
            if ( hasDirichletComponent && components[k] )
              nval[k] = 0;
          rhsLocal.axpy( quadInside[ pt ], nval );
        }
      }
    }
  }
  rhs.communicate();
}

#endif // #ifndef RHS_HH
