#include <iostream>
#include <config.h>
#include <string>
#include <sstream>

static const int dimw = Dune::GridSelector::dimworld;

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/combinedspace.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>

#include <dune/fem/misc/double.hh>


using namespace Dune;
using namespace Fem;

// polynom approximation order of quadratures,
// at least polynom order of basis functions
const int polOrd = POLORDER;

//***********************************************************************
/*! L2 Projection of a function f:

  This is an example how to solve the equation on
  \f[\Omega = (0,1)^2 \f]

  \f[ \int_{\Omega} u \phi = \int_{\Omega} f \phi  \ \ \ in \Omega \f]
  \f[ f(x,y) = x ( 1 - x) y ( 1 - y ) \f]

  Here u is the L_2 projection of f.

  The Projection should converge to the given function f.
  with the finite element method using lagrangian elements of polynom order +1.
*/
//***********************************************************************

//! the index set we are using
typedef GridSelector::GridType MyGridType;
typedef AdaptiveLeafGridPart< MyGridType > GridPartType;

// see dune/common/functionspace.hh
typedef FunctionSpace < MyGridType::ctype, double, dimw, 2 > FuncSpace1;
typedef FunctionSpace < MyGridType::ctype, double, dimw, 1 > FuncSpace2;

typedef FuncSpace1::RangeType RangeType1;
typedef FuncSpace2::RangeType RangeType2;
typedef FuncSpace1::JacobianRangeType JacobianRangeType1;
typedef FuncSpace2::JacobianRangeType JacobianRangeType2;
typedef FuncSpace1::HessianRangeType HessianRangeType1;
typedef FuncSpace2::HessianRangeType HessianRangeType2;

//! define the function space our unkown belong to
//! see dune/fem/lagrangebase.hh
typedef LagrangeDiscontinuousGalerkinSpace<FuncSpace1, GridPartType,
	polOrd+1,CachingStorage> DiscreteFunctionSpaceType1;

typedef LagrangeDiscontinuousGalerkinSpace<FuncSpace2, GridPartType,
	polOrd,CachingStorage> DiscreteFunctionSpaceType2;

typedef TupleDiscreteFunctionSpace< DiscreteFunctionSpaceType1, DiscreteFunctionSpaceType2 > DiscreteFunctionSpaceType;

typedef DiscreteFunctionSpaceType::IteratorType::Entity EntityType;

typedef DiscreteFunctionSpaceType::FunctionSpaceType FuncSpace;
typedef FuncSpace::RangeType RangeType;
typedef FuncSpace::JacobianRangeType JacobianRangeType;
typedef FuncSpace::HessianRangeType HessianRangeType;

typedef DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
typedef DiscreteFunctionSpaceType1::BasisFunctionSetType BasisFunctionSetType1;
typedef DiscreteFunctionSpaceType2::BasisFunctionSetType BasisFunctionSetType2;

typedef CachingQuadrature<GridPartType,0> QuadratureType;

const int dimRange = FuncSpace::dimRange;
const int dimRange1 = FuncSpace1::dimRange;
const int dimRange2 = FuncSpace2::dimRange;

const int dimDomain = FuncSpace::dimDomain;



int checkBasisSet ( const EntityType &entity, const DiscreteFunctionSpaceType &space )
{
  int checkSum = 0;
  BasisFunctionSetType bSet = space.basisFunctionSet( entity );

  BasisFunctionSetType1 bSet1 = space.subDiscreteFunctionSpace<0>().basisFunctionSet( entity );
  BasisFunctionSetType2 bSet2 = space.subDiscreteFunctionSpace<1>().basisFunctionSet( entity );

  const int nrBasis = bSet.size();
  const int nrBasis1 = bSet1.size();
  const int nrBasis2 = bSet2.size();

  QuadratureType quad( entity, space.order() + 1 );

  std::vector< RangeType > ranges( nrBasis, RangeType(0.) );
  std::vector< RangeType1 > ranges1( nrBasis1, RangeType1(0.) );
  std::vector< RangeType2 > ranges2( nrBasis2, RangeType2(0.) );

  std::vector< JacobianRangeType > jacs( nrBasis, JacobianRangeType(0.) );
  std::vector< JacobianRangeType1 > jacs1( nrBasis1, JacobianRangeType1(0.) );
  std::vector< JacobianRangeType2 > jacs2( nrBasis2, JacobianRangeType2(0.) );

  HessianRangeType null;
  HessianRangeType1 null1;
  HessianRangeType2 null2;

  for( int i = 0; i< dimDomain; ++i )
    for( int j = 0; j< dimDomain; ++j )
    {
      for( int r = 0; r< dimRange; ++r )
        null[r][i][j] =.0;
      for( int r = 0; r< dimRange1; ++r )
        null1[r][i][j] =.0;
      for( int r = 0; r< dimRange2; ++r )
        null2[r][i][j] =.0;
      }

  std::vector< HessianRangeType > hessians( nrBasis, null );
  std::vector< HessianRangeType1 > hessians1( nrBasis1, null1 );
  std::vector< HessianRangeType2 > hessians2( nrBasis2, null2 );

  for( unsigned int qp =0; qp < quad.nop(); ++qp )
  {

    bSet.evaluateAll( quad[qp], ranges );
    bSet1.evaluateAll( quad[qp], ranges1 );
    bSet2.evaluateAll( quad[qp], ranges2 );

    bSet.jacobianAll( quad[qp], jacs );
    bSet1.jacobianAll( quad[qp], jacs1 );
    bSet2.jacobianAll( quad[qp], jacs2 );

    bSet.hessianAll( quad[qp], hessians );
    bSet1.hessianAll( quad[qp], hessians1 );
    bSet2.hessianAll( quad[qp], hessians2 );

    RangeType value(0);
    JacobianRangeType jac(0);
    HessianRangeType hessian( null );

    for( int i=0; i< nrBasis; ++i )
    {
      for( int r = 0; r< dimRange; ++r )
      {
        if( r < dimRange1 )
        {
          for( int j=0; j < dimDomain; ++j )
          {
            jac[ r ][j] = ( i < nrBasis1 )? jacs1[i][r][j] : 0. ;

            for( int k=0; k < dimDomain; ++k )
            {
              hessian[ r ][ j ][ k ] = ( i<nrBasis1 )? hessians1[i][r][j][k] : 0. ;
            }
          }

         value[ r ] = ( i<nrBasis1 )? ranges1[ i ][ r ] : 0.;
        }
        else
        {
          value[ r ] = ( i < nrBasis1 )? 0. : ranges2[ i - nrBasis1 ][ r - dimRange1 ];

          for( int j=0; j < dimDomain; ++j )
          {
            jac[ r ][j] = ( i<nrBasis1 )? 0. : jacs2[ i-nrBasis1 ][ r - dimRange1 ][j];

            for( int k=0; k < dimDomain; ++k )
            {
              hessian[ r ][ j ][ k ] = ( i<nrBasis1 )? 0. : hessians2[ i - nrBasis1 ][ r - dimRange1 ][j][k];
            }
          }
        }
      }

      // check function value
      value -= ranges[i];
      if( value.two_norm() > 1e-8)
      {
        ++checkSum;
        std::cout<<" Basisfunction evaluation missmatch at: "<< i << " on quad point "<< qp << "  and difference: "<<  value.two_norm() <<std::endl;

        std::cout<< ranges[i]<<std::endl;
        if( i < nrBasis1 )
          std::cout<< ranges1[i]<<std::endl;
        else
          std::cout<< ranges2[i-nrBasis1]<<std::endl;
      }

      JacobianRangeType jacHelp( jac );
      // check jac
      jacHelp -= jacs[i];
      if( jacHelp.frobenius_norm() > 1e-8)
      {
        ++checkSum;
        std::cout<<" Basisfunction jacobian missmatch at: "<< i << " on quad point "<< qp << "  and difference: "<<  jacHelp.frobenius_norm() <<std::endl;

        std::cout<< jacs[i]<<std::endl;
        std::cout<< jac <<std::endl;
        if( i < nrBasis1 )
          std::cout<< jacs1[i]<<std::endl;
        else
          std::cout<< jacs2[i-nrBasis1]<<std::endl;
      }

      for( int r = 0; r < dimRange; ++r )
        hessian[ r ] -= hessians[i][ r ];

      double val = 0;
      for( int r = 0; r< dimRange; ++r )
        val += hessian[r].frobenius_norm2();

      if( std::sqrt ( val ) > 1e-8 )
      {
        ++checkSum;
        std::cout<<" Basisfunction hessian missmatch at: "<< i << " on quad point "<< qp <<std::endl;
      }
    }
  }
  return checkSum;
}

int checkSpace( const DiscreteFunctionSpaceType & space )
{
  int checksum = 0;
  for( const auto& entity : space )
    checksum += checkBasisSet( entity, space );
  return checksum;
}







//**************************************************
//
//  main programm, run algorithm twice to calc EOC
//
//**************************************************
int main (int argc, char **argv)
{
  MPIManager :: initialize( argc, argv );
  try
  {
    std::stringstream tmp;
    tmp << dimw;
    std::string macroGridName (tmp.str());
    macroGridName += "dgrid.dgf";

    GridPtr< MyGridType > gridptr( macroGridName );
    MyGridType &grid = *gridptr;

    GlobalRefine::apply(grid,2);

    GridPartType part ( grid );
    DiscreteFunctionSpaceType space ( part );

    return checkSpace( space );
  }
  catch( const Exception &exception )
  {
    std::cerr << exception << std::endl;
    return 1;
  }
}

