#include <iostream>
#include <config.h>
#include <string>
#include <sstream>

static const int dimw = Dune::GridSelector::dimworld;

#ifndef USE_BASEFUNCTIONSET_CODEGEN
#define BASEFUNCTIONSET_CODEGEN_GENERATE 1
#else
#define BASEFUNCTIONSET_CODEGEN_GENERATE 0
#endif
#ifdef USE_BASEFUNCTIONSET_CODEGEN
// include overloaded default basis function set
#include <dune/fem/space/basisfunctionset/default_codegen.hh>
#endif
#include <dune/fem/space/basisfunctionset/codegen.hh>


#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/lagrange.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/adaptationmanager.hh>

#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/misc/l2norm.hh>

#include <dune/fem/function/common/localcontribution.hh>

#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/misc/double.hh>

// to use grape, set to WANT_GRAPE to 1
#ifndef WANT_GRAPE
#define WANT_GRAPE 0
#endif

#if HAVE_GRAPE
  #define USE_GRAPE WANT_GRAPE
#else
  #define USE_GRAPE 0
  #if WANT_GRAPE
    #warning "Grape was not found by configure."
  #endif
#endif

#if USE_GRAPE
  #include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/fem/test/testgrid.hh>

using namespace Dune;
using namespace Fem;

// polynom approximation order of quadratures,
// at least polynom order of basis functions
const int polOrd = 4;

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
//typedef HierarchicGridPart< MyGridType > GridPartType;
//typedef DGAdaptiveLeafGridPart< MyGridType > GridPartType;
typedef AdaptiveLeafGridPart< MyGridType > GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
//typedef MatrixFunctionSpace < double , double, dimw , 3,5 > FuncSpace;
typedef FunctionSpace < MyGridType::ctype, double, dimw, 4 > FuncSpace;

//! define the function space our unkown belong to
//! see dune/fem/lagrangebase.hh
//typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType, polOrd, CachingStorage> DiscreteFunctionSpaceType;
typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType, polOrd > DiscreteFunctionSpaceType;

//! define the type of discrete function we are using , see
//! dune/fem/discfuncarray.hh
typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;

//! Get the Dofmanager type
typedef DofManager< MyGridType > DofManagerType;


// the exact solution to the problem for EOC calculation
struct ExactSolution
: public Fem::Function< FuncSpace, ExactSolution >
{
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::RangeFieldType RangeFieldType;
  typedef FuncSpace::DomainType DomainType;

  //! f(x,y) = x*(1-x)*y*(1-y)
  void evaluate (const DomainType & x , RangeType & ret)  const
  {
    ret = 2.; // maximum of function is 2
    for (int j=0;j<RangeType::dimension; j++)
      for(int i=0; i<DomainType::dimension; i++)
        ret[j] *= pow(x[i]*(1.0 -x[i])*4.,double(j+1));
  }

  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
  {
    evaluate ( x , ret );
  }
};

template <class Array, int slot>
struct AssignSimdFunctor
{
  explicit AssignSimdFunctor ( Array &array )
    : array_( array )
  {}

  template< class Value >
  void operator() ( const std::size_t local, const Value &value ) const
  {
    array_[ local ].insert( slot, value );
  }
private:
  Array &array_;
};


template <int i>
struct GetDofs
{
  template <class DF, class EntityVec, class Array>
  static void apply( const DF& f, const EntityVec& entities, Array& array )
  {
    AssignSimdFunctor< Array, i > af( array );
    //f.getLocalDofsFunctor( entities[ i ], af );
  }
};

template <class Array, int slot>
struct LeftAddSimdFunctor
{
  explicit LeftAddSimdFunctor ( const Array &array )
    : array_( array )
  {}

  template< class Value >
  void operator() ( const std::size_t local, Value&& value ) const
  {
    value += array_[ local ][ slot ];
  }
private:
  const Array &array_;
};

template <int i>
struct SetDofs
{
  template <class DF, class EntityVec, class Array>
  static void apply( DF& f, const EntityVec& entities, const Array& array )
  {
    // LeftAddSimdFunctor< Array, i > laf( array );
    //f.setLocalDofsFunctor( entities[ i ], laf );
  }
};

#if 0
template <int simd>
void customStandardL2Projection( const DiscreteFunctionType& f, DiscreteFunctionType& dest )
{
  dest.clear();

  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: RangeType   RangeType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename GridPartType :: template Codim< 0 > :: EntityType  EntityType;

  typedef Dune::Fem::ConstLocalFunction< DiscreteFunctionType > LocalFunctionType;
  typedef Dune::Fem::AddLocalContribution< DiscreteFunctionType > AddLocalContributionType;

  LocalFunctionType lf( f );
  AddLocalContributionType ulocal( dest );

  std::vector< typename DiscreteFunctionType::DofType > tempDof;

  typedef CachingQuadrature< GridPartType, 0 > Quadrature;

  const int order = dest.space().order() * 2;

  //typedef Vec4d                      SimdVecType;
  typedef std::array< double, simd > SimdVecType;

  typedef Dune::FieldVector< SimdVecType, RangeType :: dimension > VecRangeType;

  std::array< std::vector< RangeType >, simd > values;
  std::vector< RangeType > val;
  std::vector< VecRangeType > weights;

  std::array< EntityType, simd > entities;

  std::vector< std::vector< RangeType > > basisFcts;
  std::vector< double > dofs;

  std::vector< double > lfSimd;


  const auto endit = dest.space().end();
  for( auto it = dest.space().begin(); it != endit; )
  {
    int usedSimd = 0;

    for( int i = 0; (i<simd) && (it != endit); ++i, ++it, ++usedSimd )
    {
      entities[ i ] = *it;
      Quadrature quad( entities[ i ], order );
      const int nop = quad.nop() ;
      values[ i ].resize( nop );
      if( i == 0 )
      {
        weights.resize( nop );
        val.resize( nop );
      }

      const auto geometry = entities[ i ].geometry();
      for( int qp = 0; qp < nop; ++ qp )
        for( int d=0; d<RangeType::dimension; ++d )
          weights[ qp ][ d ][ i ] = quad.weight( qp );

      if( basisFcts.empty() )
      {
        lf.bind( entities[ 0 ] );
        Quadrature quad( entities[ 0 ], order );
        const int qnop = quad.nop();

        const int nBaseFcts = lf.basisFunctionSet().size();

        basisFcts.resize( qnop );
        std::vector< RangeType > phi( nBaseFcts );

        dofs.resize( simd * nBaseFcts );

        for( int qp=0; qp<qnop; ++qp )
        {
          lf.basisFunctionSet().evaluateAll( quad[ qp ], phi );
          auto& baseFct = basisFcts[ qp ];
          baseFct.resize( nBaseFcts );

          for( int b=0; b<nBaseFcts; ++b )
          {
            for( int d=0; d<RangeType::dimension; ++d )
            {
              //baseFct[ b ][ d ].insert( i, phi[ b ][ d ] );
              baseFct[ b ][ d ] = phi[ b ][ d ];
            }
          }
        }
        lf.unbind();
      }

    }

    const int baseSize = basisFcts[ 0 ].size();
    Quadrature quad( entities[ 0 ], order );

    lfSimd.resize( baseSize * simd );
    // evaluate localfunctions of f
    for( int i = 0; i<usedSimd; ++i )
    {
      double* lfVal = lfSimd.data() + i*baseSize;
      f.getLocalDofs( entities[ i ], lfVal );
    }

    double* lf0 = lfSimd.data();
    double* lf1 = lf0 + baseSize;
    double* lf2 = lf1 + baseSize;
    double* lf3 = lf2 + baseSize;

    for( size_t qp = 0; qp < quad.nop(); ++qp )
    {
      auto& baseFct = basisFcts[ qp ];
      for( int i=0; i<simd; ++i )
        values[ i ][ qp ] = 0;

      for( int b=0; b<baseSize; ++b )
      {
        for( int d=0; d<RangeType::dimension; ++d )
        {
          values[ 0 ][ qp ][ d ] += lf0[ b ] * baseFct[ b ][ d ];
          values[ 1 ][ qp ][ d ] += lf1[ b ] * baseFct[ b ][ d ];
          values[ 2 ][ qp ][ d ] += lf2[ b ] * baseFct[ b ][ d ];
          values[ 3 ][ qp ][ d ] += lf3[ b ] * baseFct[ b ][ d ];
        }
      }
      //for( int d=0; d<RangeType::dimension; ++d )
      //  values[ qp ][ d ].insert( i, val[ qp ][ d ] );
        //values[ qp ][ d ][ i ] = val[ qp ][ d ];
    }

    // values * weights
    const int quadNop = values.size();
    for( int qp = 0; qp < quadNop; ++qp )
      for( int d=0; d<RangeType::dimension; ++d )
        for( int i=0; i<usedSimd; ++i )
          values[ i ][ qp ][ d ] *= weights[ qp ][ d ][ i ];


    // axpy operation
    {
      const int quadNop  = values.size();
      const int baseSize = basisFcts[ 0 ].size();

      for( int b=0; b<simd * baseSize; ++b )
      {
        dofs[ b ] = 0.0;
      }

      double* dof0 = dofs.data();
      double* dof1 = dof0 + baseSize;
      double* dof2 = dof1 + baseSize;
      double* dof3 = dof2 + baseSize;

      for( int qp = 0; qp< quadNop; ++qp )
      {
        auto& baseFct = basisFcts[ qp ];
        const auto& val0 = values[ 0 ][ qp ];
        const auto& val1 = values[ 1 ][ qp ];
        const auto& val2 = values[ 2 ][ qp ];
        const auto& val3 = values[ 3 ][ qp ];

        for( int b=0; b<baseSize; ++b )
        {
          for( int d=0; d<RangeType::dimension; ++d )
          {
            //dofs[ b ][ d ] += values[ qp ][ d ] * baseFct[ b ][ d ];
            //for( int i=0; i<usedSimd; ++i )
            dof0[ b ] += val0[ d ] * baseFct[ b ][ d ];
            dof1[ b ] += val1[ d ] * baseFct[ b ][ d ];
            dof2[ b ] += val2[ d ] * baseFct[ b ][ d ];
            dof3[ b ] += val3[ d ] * baseFct[ b ][ d ];
            //  dofs[ b ][ d ][ i ] += val[ d ][ i ] * baseFct[ b ][ d ];
          }
        }
      }
    }

    for( int i = 0; i<usedSimd; ++ i )
    {
      ulocal.bind( entities[ i ] );
      int dof = 0;
      double* dofVec = dofs.data() + i*baseSize;
      for( int b = 0; b < baseSize; ++b )
      {
        for( int d=0; d<RangeType::dimension; ++d, ++dof )
        {
          ulocal[ dof ] = dofVec[ b ];
        }
      }
      ulocal.unbind();
    }

    if( it == endit ) break ;
  }
}

void customL2Projection( const DiscreteFunctionType& f, DiscreteFunctionType& dest )
{
  dest.clear();

  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: RangeType   RangeType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename GridPartType :: template Codim< 0 > :: EntityType  EntityType;

  static const int simd = 4;

  typedef Dune::Fem::ConstLocalFunction< DiscreteFunctionType > LocalFunctionType;
  typedef Dune::Fem::AddLocalContribution< DiscreteFunctionType > AddLocalContributionType;

  LocalFunctionType lf( f );
  AddLocalContributionType ulocal( dest );

  std::vector< typename DiscreteFunctionType::DofType > tempDof;

  typedef CachingQuadrature< GridPartType, 0 > Quadrature;

  const int order = dest.space().order() * 2;

  typedef Vec4d                      SimdVecType;
  //typedef std::array< double, simd > SimdVecType;

  typedef Dune::FieldVector< SimdVecType, RangeType :: dimension > VecRangeType;

  std::vector< VecRangeType > values;
  std::vector< RangeType > val;
  std::vector< SimdVecType > weights;

  std::array< EntityType, simd > entities;

  std::vector< std::vector< VecRangeType > > basisFcts;
  std::vector< SimdVecType > dofs;

  typedef std::vector< SimdVecType > SimdArrayType;
  SimdArrayType lfSimd;

  const auto endit = dest.space().end();
  for( auto it = dest.space().begin(); it != endit; )
  {
    int usedSimd = 0;

    for( int i = 0; (i<simd) && (it != endit); ++i, ++it, ++usedSimd )
    {
      entities[ i ] = *it;
      Quadrature quad( entities[ i ], order );
      const int nop = quad.nop() ;
      if( i == 0 )
      {
        values.resize( nop );
        weights.resize( nop );
        val.resize( nop );
      }

      const auto geometry = entities[ i ].geometry();
      for( int qp = 0; qp < nop; ++ qp )
        weights[ qp ].insert( i, quad.weight( qp ) );

      if( basisFcts.empty() )
      {
        lf.bind( entities[ 0 ] );
        Quadrature quad( entities[ 0 ], order );
        const int qnop = quad.nop();

        const int nBaseFcts = lf.basisFunctionSet().size();

        basisFcts.resize( qnop );
        std::vector< RangeType > phi( nBaseFcts );

        dofs.resize( nBaseFcts );
        lfSimd.resize( nBaseFcts );

        for( int qp=0; qp<qnop; ++qp )
        {
          lf.basisFunctionSet().evaluateAll( quad[ qp ], phi );
          auto& baseFct = basisFcts[ qp ];
          baseFct.resize( nBaseFcts );

          for( int b=0; b<nBaseFcts; ++b )
          {
            for( int d=0; d<RangeType::dimension; ++d )
            {
              for( int i = 0; i<simd; ++i )
                baseFct[ b ][ d ].insert( i, phi[ b ][ d ] );
                //baseFct[ b ][ d ][ i ] = phi[ b ][ d ];
            }
          }
        }
        lf.unbind();
      }

    }

    const int baseSize = basisFcts[ 0 ].size();
    Quadrature quad( entities[ 0 ], order );

    Dune::Fem::ForLoop< GetDofs, 0, simd-1 >::apply( f, entities, lfSimd );

    //lf.evaluateQuadrature( quad, val );
    const int singleBase = baseSize / RangeType :: dimension;
    for( size_t qp = 0; qp < quad.nop(); ++qp )
    {
      auto& baseFct = basisFcts[ qp ];
      values[ qp ] = SimdVecType( 0 );
      int dof = 0;
      for( int b=0; b<singleBase; ++b )
      {
        for( int d=0; d<RangeType::dimension; ++d, ++dof  )
          values[ qp ][ d ] += lfSimd[ dof ] *  baseFct[ dof ][ d ];
        //values[ qp ].axpy( lfSimd[ b ], baseFct[ b ] );
        //+=
        //for( int d=0; d<RangeType::dimension; ++d )
       // {
        //  values[ qp ][ d ] += lfSimd[ b ] * baseFct[ b ][ d ];
        //}
      }
      //for( int d=0; d<RangeType::dimension; ++d )
      //  values[ qp ][ d ].insert( i, val[ qp ][ d ] );
        //values[ qp ][ d ][ i ] = val[ qp ][ d ];
    }

    // values * weights
    const int quadNop = values.size();
    for( int qp = 0; qp < quadNop; ++qp )
        //values[ qp ][ d ][ i ] *= weights[ qp ][ d ][ i ];
        values[ qp ] *= weights[ qp ];


    // axpy operation
    {
      const int quadNop  = values.size();
      const int baseSize = basisFcts[ 0 ].size();

      for( int b=0; b<baseSize; ++b )
      {
        dofs[ b ] = 0.0;
      }

      for( int qp = 0; qp< quadNop; ++qp )
      {
        auto& baseFct = basisFcts[ qp ];
        const auto& val = values[ qp ];
        int dof = 0;
        for( int b=0; b<singleBase; ++b )
        {
          for( int d=0; d<RangeType::dimension; ++d, ++dof )
            dofs[ dof ] += val[ d ] * baseFct[ dof ][ d ];
        }
      }
    }

    Dune::Fem::ForLoop< SetDofs, 0, simd-1 >::apply( dest, entities, dofs );
    /*
    for( int i = 0; i<usedSimd; ++ i )
    {
      // dest.setLocalDofs(
      ulocal.bind( entities[ i ] );
      int dof = 0;
      for( int b = 0; b < baseSize; ++b )
      {
        for( int d=0; d<RangeType::dimension; ++d, ++dof )
        {
          ulocal[ dof ] = dofs[ b ][ d ][ i ];
        }
      }
      ulocal.unbind();
    }
    */

    if( it == endit ) break ;
  }
}
#endif

void customInterpolate( const DiscreteFunctionType& f, DiscreteFunctionType& dest )
{
  ConstLocalFunction< DiscreteFunctionType > uLocal( f );
  AddLocalContribution< DiscreteFunctionType > vLocal( dest );

  for( const auto &entity : dest.space() )
  {
    auto uGuard = bindGuard( uLocal, entity );
    auto vGuard = bindGuard( vLocal, entity );

    // interpolate u
    dest.space().interpolation( entity )( uLocal, vLocal );
  }

}

// ********************************************************************
double algorithm ( MyGridType &grid, DiscreteFunctionType &solution, bool display )
{
   ExactSolution f;
   // create exact solution for error evaluation
   typedef GridFunctionAdapter< ExactSolution, GridPartType >  GridFunctionType;
   GridFunctionType exactSolution( "exact solution", f, solution.gridPart(), solution.space().order() );

   // L2 error class
   Dune :: Fem :: L2Norm< GridPartType > l2norm( solution.gridPart() );

   DiscreteFunctionType copy( solution );

   //! perform l2-projection
   interpolate(exactSolution, solution);

   Dune::Timer  timer;
   //! perform l2-projection
   interpolate(solution, copy);

   std::cout << "Standard Interpolation: " << timer.elapsed() << std::endl;

   copy.clear();

   timer.reset();
   //! perform l2-projection
   customInterpolate(solution, copy);

   std::cout << "Custom  Interpolation: " << timer.elapsed() << std::endl;

   /*
   timer.reset();
   // do l2 projection for one df to another
   customStandardL2Projection< 8 >( solution, copy );
   std::cout << "Fake simd 8 interpolation: " << timer.elapsed() << std::endl;
   */

   /*
   timer.reset();
   // do l2 projection for one df to another
   customStandardL2Projection< 4 >( solution, copy );
   std::cout << "Fake simd 4 interpolation: " << timer.elapsed() << std::endl;
   */

   /*
   timer.reset();
   // do l2 projection for one df to another
   customStandardL2Projection< 1 >( solution, copy );
   std::cout << "Fake simd 1 interpolation: " << timer.elapsed() << std::endl;
   */

   /*
   timer.reset();
   // do l2 projection for one df to another
   customL2Projection( solution, copy );
   std::cout << "Simd interpolation: " << timer.elapsed() << std::endl;
   */

   //interpolate( solution, copy );

   // calculation L2 error
   // pol ord for calculation the error should be higher than
   // pol for evaluation the basefunctions
   double error = l2norm.distance( exactSolution, solution );
   std::cout << "\nL2 Error: " << error << "\n\n";
   double error2 = l2norm.distance( exactSolution, copy );
   std::cout << "\nL2 Error(copy): " << error2 << "\n\n";

   Dune::Fem::VTKIO< GridPartType > output( solution.gridPart(), Dune::VTK::nonconforming );
   {
     output.addCellData( solution );
     output.write( "l2projection" );
   }

#if USE_GRAPE
   // if Grape was found, then display last solution
   if( display )
   {
     GrapeDataDisplay< MyGridType > grape( solution.space().gridPart() );
     grape.dataDisplay( solution );
   }
#endif

   return error;
}

template< class GridPartType >
void generateCode ( GridPartType &gridPart )
{
  DiscreteFunctionSpaceType space( gridPart );
  std::vector< int > elemQuadOrds = {{ 0, 1, 2, 3 ,4, space.order() + 2, space.order(), space.order()*2, space.order() * 2 + 2, space.order()+1, space.order() * 2 +1  }};
  std::vector< int > faceQuadOrds = {{ space.order(), space.order()*2, space.order() * 2 + 2, space.order()+1 }};
  const std::string path = ".";

  Dune::Fem::generateCode( space, elemQuadOrds, faceQuadOrds, path );
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

  int ml = 1;
  if( argc > 1 )
    ml = atoi( argv[1] );

  /*
  if( ml == -1 )
  {
    MyGridType &grid = Dune::Fem::TestGrid::grid();
    GridPartType part ( grid );
    generateCode( part );
    return 0;
  }
  */

  std::vector< double> error(ml);

  MyGridType &gridP = Dune::Fem::TestGrid::grid();
  //MyGridType grid( gridP.dualGrid() );
  MyGridType& grid = gridP;
  const int step = Dune::Fem::TestGrid::refineStepsForHalf();

  GridPartType part ( grid );
  DiscreteFunctionSpaceType linFuncSpace ( part );
  DiscreteFunctionType solution ( "sol", linFuncSpace );
  solution.clear();

  for(int i=0; i<ml; i+=step)
  {
    std::cout << "Global refine" << std::endl;
    GlobalRefine::apply(grid,step);
    std::cout << "Algorithm" << std::endl;
    error[i] = algorithm ( grid , solution , i==ml-1);
    if (i>0)
    {
      double eoc = log( error[i-step]/error[i]) / M_LN2;
      std::cout << "EOC = " << eoc << " \n";
    }
  }
  return 0;
  }
  catch( const Exception &exception )
  {
    std::cerr << exception << std::endl;
    return 1;
  }
}

