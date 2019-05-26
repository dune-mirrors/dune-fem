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
#ifndef DUNE_FEM_EDGESPACE_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_EDGESPACE_SHAPEFUNCTIONSET_HH

#include <dune/common/exceptions.hh>

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/simple.hh>

namespace Dune
{

  namespace Fem
  {

    // SimplexBubbleElementShapeFunctionSet
    // ----------------------------------

    template< class FunctionSpace >
    class SimplexBubbleElementShapeFunctionSet
    {
      typedef SimplexBubbleElementShapeFunctionSet< FunctionSpace > ThisType;
    public:
      static const int dimDomain =  FunctionSpace :: dimDomain;
      static const int polynomialOrder = dimDomain + 1;
      static const int numShapeFunctions = dimDomain + 2;

      typedef FunctionSpace FunctionSpaceType;
      typedef typename FunctionSpace :: DomainType DomainType;
      typedef typename FunctionSpace :: RangeType RangeType;
      static_assert( RangeType::dimension == 1, "This is a scalar shapefunction set" );
      typedef typename FunctionSpace :: JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpace :: HessianRangeType HessianRangeType;

      SimplexBubbleElementShapeFunctionSet () {}

      template<class GeometryType >
      SimplexBubbleElementShapeFunctionSet ( const GeometryType& gt )
      {
        if( !gt.isSimplex() )
          DUNE_THROW( NotImplemented, "Wrong geometry type for Cube2DBubblelementShapeFunctionSet." );
      }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        DomainType xRef = coordinate( x );
        RangeType phi(1), phi0(1);
        for( int i=0; i< dimDomain; ++i )
        {
          functor( i+1, RangeType( xRef[ i ] ) );
          phi0[ 0 ] -= xRef[ i ];
          phi[ 0 ] *= xRef[ i ] ;
        }

        phi[ 0 ] *= phi0[ 0 ] / std::pow( ( dimDomain + 1.0 ), dimDomain + 1.0 );
        functor( 0, phi0 );
        functor( dimDomain +1, phi );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        DomainType xRef = coordinate( x );

        JacobianRangeType jac(0), jac0( -1 );
        RangeType phi0( 1 );

        functor( 0, jac0 );

        for( int i=0; i< dimDomain; ++i )
        {
          phi0[ 0 ] -= xRef[ i ];

          for( int j=1; j < dimDomain; ++j )
            jac0[ 0 ][ (i+j)%dimDomain ] *= xRef[ i ];

          jac[ 0 ][ i ] = 1;
          functor( i+1, jac );
          jac[ 0 ][ i ] = 0;
        }

        for( int i=0; i< dimDomain; ++i )
          jac0[ 0 ][ i ] *= -(phi0[ 0 ] - xRef[ i ]);
        jac0[ 0 ] *= 1.0 / std::pow( dimDomain + 1.0, dimDomain + 1.0 );
        functor( dimDomain +1, jac0 );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        DUNE_THROW( NotImplemented, "NotImplemented" );
        DomainType xRef = coordinate( x );
        HessianRangeType hes;
        functor( 0, hes );
        functor( 1, hes );
        functor( 2, hes );
        functor( 3, hes );
      }

      int order () const { return dimDomain + 1; }

      std::size_t size () const { return dimDomain +2; }
    };

    // 2d cube element taken from
    // http://arxiv.org/pdf/1507.04417v1.pdf
    template< class FunctionSpace >
    class Cube2DBubbleElementShapeFunctionSet
    {
      typedef Cube2DBubbleElementShapeFunctionSet< FunctionSpace > ThisType;
    public:
      static const int dimDomain =  FunctionSpace :: dimDomain;
      static const int polynomialOrder = dimDomain + 1;
      static const int numShapeFunctions = dimDomain + 2;

      typedef FunctionSpace FunctionSpaceType;
      typedef typename FunctionSpace :: DomainType DomainType;
      typedef typename FunctionSpace :: RangeType RangeType;
      static_assert( RangeType::dimension == 1, "This is a scalar shapefunction set" );
      typedef typename FunctionSpace :: JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpace :: HessianRangeType HessianRangeType;

      //! [Main methods for shape functions]
      Cube2DBubbleElementShapeFunctionSet () {}

      template<class GeometryType >
      Cube2DBubbleElementShapeFunctionSet ( const GeometryType& gt )
      {
        if( !gt.isCube() )
          DUNE_THROW( NotImplemented, "Wrong geometry type for Cube2DBubblelementShapeFunctionSet." );
        if( gt.dim() != 2 )
          DUNE_THROW( NotImplemented, "2DCubeBubbleElementShapeFunctionSet only implemented for dimension = 2." );
      }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      //! [Main methods for shape functions]
      {
        DomainType xRef = coordinate( x );
        RangeType phi(1), phi0(1);
        // (1-x)(1-y)  -> grad = ( (y-1),(x-1) )
        phi0[0] = (1.-xRef[0])*(1.-xRef[1]);
        functor( 0, phi0 );
        // x(1-y)      -> grad = ( (1-y),-x )
        phi[0] = xRef[0]*(1.-xRef[1]);
        functor( 1, phi );
        // (1-x)y      -> grad = ( (-y,(1-x) )
        phi[0] = (1.-xRef[0])*xRef[1];
        functor( 2, phi );
        // xy          -> grad = ( (y,x) )
        phi[0] = xRef[0]*xRef[1];
        functor( 3, phi );
        // 64 xy phi_0(x,y)^2       -> grad = 128 phi_0 xy grad_0 +
        //                                     64 phi_0^2 (y,x)
        phi[0] = 64.*phi0[0]*phi0[0]*xRef[0]*xRef[1];
        functor( 4, phi );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        DomainType xRef = coordinate( x );

        JacobianRangeType jac, jac0;
        RangeType phi0;
        phi0[0] = (1.-xRef[0])*(1.-xRef[1]);
        jac0[0] = {xRef[1]-1,xRef[0]-1};

        functor( 0, jac0 );
        jac[0] = {1-xRef[1],xRef[0]};
        functor( 1, jac );
        jac[0] = {xRef[1],1-xRef[0]};
        functor( 2, jac );
        jac[0] = {xRef[1],xRef[0]};
        functor( 3, jac );
        // 64 xy phi_0(x,y)^2       -> grad = 128 phi_0 xy grad_0 +
        //                                     64 phi_0^2 (y,x)
        jac0[0] *= 128.*phi0*xRef[0]*xRef[1];
        jac[0] = {xRef[1],xRef[0]};
        jac[0] *= 64.*phi0*phi0;
        jac[0] += jac0[0];
        functor( 4, jac );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        DUNE_THROW( NotImplemented, "NotImplemented" );
        DomainType xRef = coordinate( x );
        HessianRangeType hes;
        functor( 0, hes );
        functor( 1, hes );
        functor( 2, hes );
        functor( 3, hes );
      }

      int order () const { return 4; }

      std::size_t size () const { return 5; }
    };
  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_EDGESPACE_SHAPEFUNCTIONSET_HH
