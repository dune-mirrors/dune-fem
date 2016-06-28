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
#ifndef SPACE_LOCALKEYMAP_HH
#define SPACE_LOCALKEYMAP_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/space/mapper/code.hh>
#include <dune/fem/space/mapper/localkey.hh>

namespace Dune
{

  namespace Fem
  {

    template< int dim >
    struct BubbleElementLocalKeyMap
    {
      //! [Constructor of LocalKey tripple]
      BubbleElementLocalKeyMap ( int vertices )
      {
        for( int i = 0; i < vertices; ++i )
          map_.emplace_back( i, dim, 0 );
        map_.emplace_back( 0, 0, 0 );
      }
      //! [Constructor of LocalKey tripple]

      std::size_t size() const { return map_.size(); }

      LocalKey& localKey ( std::size_t i ) { return map_[ i ]; }
      const LocalKey& localKey ( std::size_t i ) const { return map_[ i ]; }

    private:
      std::vector< LocalKey > map_;
    };

    struct BubbleElementDofMapperCodeFactory
    {
      // return the shape functions for a given reference element. If this
      // is not possible an empty DofMapperCode is returned.
      template< class Field, int dim >
      DofMapperCode operator() ( const ReferenceElement< Field, dim > &refElement ) const
      {
        if( refElement.type().isSimplex() )
          return compile( refElement, BubbleElementLocalKeyMap< dim >(dim+1) );
        if( refElement.type().isCube() && refElement.type().dim() == 2)
          return compile( refElement, BubbleElementLocalKeyMap< dim >(pow(2,dim)) );
        else
          return DofMapperCode();
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef SPACE_LOCALKEYMAP_HH
