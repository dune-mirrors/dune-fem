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
#ifndef SPACE_LOCALINTERPOLATION_HH
#define SPACE_LOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{

  namespace Fem
  {

    template< class FunctionSpace >
    struct LocalBubbleElementInterpolation
    {
      typedef typename FunctionSpace::DomainType DomainType;
      typedef typename FunctionSpace::RangeType RangeType;
      static const int dimDomain = FunctionSpace::dimDomain;
      static const int dimRange = FunctionSpace::dimRange;

      LocalBubbleElementInterpolation ()
        : points_( dimDomain + 2, DomainType( 0.0 ) )
      {
        for( int i = 0; i < dimDomain; ++i )
          points_[ i + 1 ][ i ] = 1.0;

        points_[ dimDomain +1 ] = DomainType( 1.0 / ( dimDomain + 1.0 ) );
      }

      LocalBubbleElementInterpolation ( const LocalBubbleElementInterpolation & ) = default;
      LocalBubbleElementInterpolation ( LocalBubbleElementInterpolation && ) = default;

      //! [Evaluation of local interpolations]
      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &lf, LocalDofVector &ldv ) const
      //! [Evaluation of local interpolations]
      {
        int k = 0;
        for( const DomainType &x : points_ )
        {
          RangeType phi;
          lf.evaluate( x, phi );
          for( int i = 0; i < dimRange; ++i )
            ldv[ k++ ] = phi[ i ];
        }
      }

    private:
      std::vector< DomainType > points_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef SPACE_LOCALINTERPOLATION_HH
