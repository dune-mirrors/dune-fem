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

#include <dune/fem/misc/double.hh>

#if HAVE_DUNE_VECTORCLASS
#include <dune/vectorclass/vectorclass.hh>
#endif
#include <dune/common/simd/vc.hh>

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
const int polOrd = 2;

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

#ifdef COUNT_FLOPS
typedef Dune::Fem::Double Field;
#else
typedef double Field;
#endif

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
//typedef MatrixFunctionSpace < double , double, dimw , 3,5 > FuncSpace;
typedef FunctionSpace < MyGridType::ctype, Field, dimw, 1 > FuncSpace;

//! define the function space our unkown belong to
//! see dune/fem/lagrangebase.hh
//typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType, polOrd, CachingStorage> DiscreteFunctionSpaceType;
typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType, polOrd, CachingStorage> DiscreteFunctionSpaceType;
//typedef HierarchicLegendreDiscontinuousGalerkinSpace<FuncSpace, GridPartType, polOrd, CachingStorage> DiscreteFunctionSpaceType;

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

/*
namespace Simd
{
  template <int d>
  double& lane( const int slot, FieldVector< double, d>& vec)
  {
    return vec[ slot ];
  }

  template <int d>
  const double& lane( const int slot, const FieldVector< double, d>& vec)
  {
    return vec[ slot ];
  }

};
*/

template <class Array, int slot>
struct AssignSimdFunctor
{
  explicit AssignSimdFunctor ( Array &array )
    : array_( array )
  {}

  template< class Value >
  void operator() ( const std::size_t local, const Value &value ) const
  {
    Simd::lane(slot, array_[ local ]) = value;
  }
private:
  Array &array_;
};


template <int i>
struct GetDofs
{
  template <class DF, class EntityVec, class Array>
  static void apply( const DF& f, const EntityVec& entityPtr, Array& array )
  {
    if( entityPtr[ i ] )
    {
      AssignSimdFunctor< Array, i > af( array );
      f.getLocalDofsFunctor( *entityPtr[ i ], af );
    }
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
    value += Simd::lane(slot, array_[ local ]);
  }
private:
  const Array &array_;
};

template <int i>
struct SetDofs
{
  template <class DF, class EntityVec, class Array>
  static void apply( DF& f, const EntityVec& entityPtr, const Array& array )
  {
    if( entityPtr[ i ] )
    {
      LeftAddSimdFunctor< Array, i > laf( array );
      f.setLocalDofsFunctor( *entityPtr[ i ], laf );
    }
  }
};

template<class SimdVecType>
void customL2Projection( const DiscreteFunctionType& f, DiscreteFunctionType& dest )
{
  dest.clear();

  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: RangeType       RangeType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType  RangeFieldType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename GridPartType :: template Codim< 0 > :: EntityType  EntityType;


  typedef Dune::Fem::ConstLocalFunction< DiscreteFunctionType > LocalFunctionType;
  typedef Dune::Fem::AddLocalContribution< DiscreteFunctionType > AddLocalContributionType;

  LocalFunctionType lf( f );
  AddLocalContributionType ulocal( dest );

  std::vector< typename DiscreteFunctionType::DofType > tempDof;

  typedef CachingQuadrature< GridPartType, 0 > Quadrature;

  const int order = dest.space().order() * 2;

  constexpr int lanes= Simd::lanes<SimdVecType>();

  typedef Dune::FieldVector< Simd::Scalar<SimdVecType>, RangeType :: dimension > VecRangeType;
  typedef Dune::FieldVector< SimdVecType, RangeType :: dimension > SimdVecRangeType;

  std::vector< SimdVecRangeType > values;
  std::vector< RangeType > val;
  std::vector< SimdVecType > weights;

  std::array< EntityType, lanes > entities;
  std::array< const EntityType*, lanes > entityPtr;

  std::vector< std::vector< RangeFieldType > > basisFcts;
  typedef std::vector< SimdVecType > SimdArrayType;
  SimdArrayType dofs;

  SimdArrayType lfSimd;

  Quadrature quad( *dest.space().begin(), order );
  const int nop = quad.nop() ;
  std::cout << "lanes = " << lanes << " qp = " << nop << std::endl;
  lf.bind( *dest.space().begin() );
  const int nBaseFcts = lf.basisFunctionSet().size();
  const int singleBase = nBaseFcts / RangeType :: dimension;
  basisFcts.resize( nop );
  std::vector< RangeType > phi( nBaseFcts );
  dofs.resize( nBaseFcts );
  lfSimd.resize( nBaseFcts );
  for( int qp=0; qp<nop; ++qp )
  {
    lf.basisFunctionSet().evaluateAll( quad[ qp ], phi );
    auto& baseFct = basisFcts[ qp ];
    baseFct.resize( nBaseFcts );

    for( int b=0; b<singleBase; ++b )
    {
      baseFct[ b ] = phi[ b*RangeType :: dimension ][ 0 ];
    }
  }
  lf.unbind();


  const auto endit = dest.space().end();
  /*
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  std::array< IteratorType, lanes > iterator;
  std::array< IteratorType, lanes > endIterator;
  iterator[ 0 ] = dest.space().begin();

  std::cout << "size = " << dest.space().gridPart().indexSet().size( 0 ) << std::endl;
  std::cout << "Lanes = " << lanes << std::endl;
  int count = 0;
  int quadCount = dest.space().gridPart().indexSet().size( 0 ) / lanes;
  int lane = 1;
  for( auto it = dest.space().begin(); it != endit; ++it, ++count )
  {
    if( count == quadCount )
    {
      std::cout << "Save lane = " << lane << " at count = " << count << std::endl;
      iterator[ lane ] = it;
      endIterator[ lane-1 ] = it;
      count = 0 ;
      ++lane;
    }
  }
  endIterator[ lanes-1 ] = endit;
  */

  std::cout << "lfSimd " << lfSimd.size() * sizeof(SimdVecType) << std::endl;
  std::cout << "dofs   " << dofs.size() * sizeof(SimdVecType) << std::endl;
  std::cout << "base   " << singleBase * nop *sizeof(double) << std::endl;

  std::cout << "Cache = " << (lfSimd.size() * sizeof(SimdVecType) + dofs.size() *
      sizeof(SimdVecType)+ singleBase * nop *sizeof(double) ) << std::endl;

  //while( iterator[ 0 ] != endIterator[ 0 ] )
  for( auto it = dest.space().begin(); it != endit; )
  {
    for( int i = 0; i<lanes; ++i ) entityPtr[ i ] = nullptr;;
    for( int i = 0; (i<lanes) && (it != endit); ++i, ++it)
    //for( int i = 0; i<lanes; ++i )
    {
      //entityPtr[ i ] = nullptr;
      //if( iterator[ i ] != endIterator[ i ])
      {
        //entities[ i ] = *(iterator[ i ]);
        entities[ i ] = *(it);

        entityPtr[ i ] = &entities[ i ];
        Quadrature quad( entities[ i ], order );
        const int nop = quad.nop() ;
        if( i == 0 )
        {
          values.resize( nop );
        }
      }
    }

    Quadrature quad( entities[ 0 ], order );

    Dune::Fem::ForLoop< GetDofs, 0, lanes-1 >::apply( f, entityPtr, lfSimd );

    //lf.evaluateQuadrature( quad, val );
    for( size_t qp = 0; qp < quad.nop(); ++qp )
    {
      auto& baseFct = basisFcts[ qp ];
      values[ qp ] = SimdVecType( 0 );
      int dof = 0;
      for( int b=0; b<singleBase; ++b )
      {
        for( int d=0; d<RangeType::dimension; ++d, ++dof  )
        {
          // TODO: use mul_add
          values[ qp ][ d ] += lfSimd[ dof ] * baseFct[ b ];
          //values[ qp ][ d ] = mul_add( lfSimd[ dof ], baseFct[ dof ][ d ], values[ qp ][ d ]);
        }
      }
    }

    // values * weights
    const int quadNop = quad.nop();
    for( int qp = 0; qp < quadNop; ++qp )
    {
      //values[ qp ][ d ][ i ] *= weights[ qp ][ d ][ i ];
      values[ qp ] *= quad.weight( qp );
    }

    // axpy operation
    {
      const int baseSize = singleBase * RangeType :: dimension;

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
          {
            dofs[ dof ] += val[ d ] * baseFct[ b ];
            //dofs[ dof ] = mul_add( val[ d ], baseFct[ dof ][ d ], dofs[ dof ]);
          }
        }
      }
    }

    Dune::Fem::ForLoop< SetDofs, 0, lanes-1 >::apply( dest, entityPtr, dofs );
    //SetDofsSimd::apply( dest, entityPtr, dofs );

    if( it == endit ) break ;
    /*
    for( int i = 0; i<lanes; ++i )
    {
      if( iterator[ i ] != endIterator[ i ])
        ++iterator[ i ];
    }
    */
  }
}

template<int lanes>
void customL2ProjectionNormal( const DiscreteFunctionType& f, DiscreteFunctionType& dest )
{
  dest.clear();

  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: RangeType       RangeType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType  RangeFieldType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef typename GridPartType :: template Codim< 0 > :: EntityType  EntityType;


  typedef Dune::Fem::ConstLocalFunction< DiscreteFunctionType > LocalFunctionType;
  typedef Dune::Fem::AddLocalContribution< DiscreteFunctionType > AddLocalContributionType;

  LocalFunctionType lf( f );
  AddLocalContributionType ulocal( dest );

  std::vector< typename DiscreteFunctionType::DofType > tempDof;

  typedef CachingQuadrature< GridPartType, 0 > Quadrature;

  const int order = dest.space().order() * 2;

  typedef Dune::FieldVector< RangeFieldType, lanes > SimdVecType;
  typedef Dune::FieldVector< SimdVecType, RangeType :: dimension > SimdVecRangeType;

  //std::vector< SimdVecRangeType > values;
  std::vector< RangeType > values;

  std::array< EntityType, lanes > entities;
  std::array< const EntityType*, lanes > entityPtr;

  std::vector< std::vector< RangeFieldType > > basisFcts;
  std::vector< RangeFieldType > dofs;

  typedef std::vector< RangeFieldType > SimdArrayType;
  SimdArrayType lfSimd;

  Quadrature quad( *dest.space().begin(), order );
  const int nop = quad.nop() ;
  lf.bind( *dest.space().begin() );
  const int nBaseFcts = lf.basisFunctionSet().size();
  const int singleBase = nBaseFcts / RangeType :: dimension;
  basisFcts.resize( nop );
  std::vector< RangeType > phi( nBaseFcts );
  dofs.resize( nBaseFcts * lanes );
  lfSimd.resize( nBaseFcts * lanes );
  for( int qp=0; qp<nop; ++qp )
  {
    lf.basisFunctionSet().evaluateAll( quad[ qp ], phi );
    auto& baseFct = basisFcts[ qp ];
    baseFct.resize( nBaseFcts );

    for( int b=0; b<singleBase; ++b )
    {
      baseFct[ b ] = phi[ b*RangeType :: dimension ][ 0 ];
    }
  }
  lf.unbind();

  const auto endit = dest.space().end();
  /*
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  std::array< IteratorType, lanes > iterator;
  std::array< IteratorType, lanes > endIterator;
  iterator[ 0 ] = dest.space().begin();

  std::cout << "size = " << dest.space().gridPart().indexSet().size( 0 ) << std::endl;
  std::cout << "Lanes = " << lanes << std::endl;
  int count = 0;
  int quadCount = dest.space().gridPart().indexSet().size( 0 ) / lanes;
  int lane = 1;
  for( auto it = dest.space().begin(); it != endit; ++it, ++count )
  {
    if( count == quadCount )
    {
      std::cout << "Save lane = " << lane << " at count = " << count << std::endl;
      iterator[ lane ] = it;
      endIterator[ lane-1 ] = it;
      count = 0 ;
      ++lane;
    }
  }
  endIterator[ lanes-1 ] = endit;
  */

  //while( iterator[ 0 ] != endIterator[ 0 ] )
  for( auto it = dest.space().begin(); it != endit; )
  {
    for( int i = 0; i<lanes; ++i ) entityPtr[ i ] = nullptr;;
    for( int i = 0; (i<lanes) && (it != endit); ++i, ++it)
    //for( int i = 0; i<lanes; ++i )
    {
      //entityPtr[ i ] = nullptr;
      //if( iterator[ i ] != endIterator[ i ])
      {
        //entities[ i ] = *(iterator[ i ]);
        entities[ i ] = *(it);

        entityPtr[ i ] = &entities[ i ];
        Quadrature quad( entities[ i ], order );
        const int nop = quad.nop() ;
        if( i == 0 )
        {
          values.resize( nop * lanes);
        }
      }
    }

    Quadrature quad( entities[ 0 ], order );

    //Dune::Fem::ForLoop< GetDofs, 0, lanes-1 >::apply( f, entityPtr, lfSimd );
    for( int i=0; i<lanes; ++i )
    {
      if( entityPtr[ i ] )
      {
        Dune::Fem::StaticArray< double > wrapper( nBaseFcts, lfSimd.data() + nBaseFcts * i);
        f.getLocalDofs( *(entityPtr[ i ]), wrapper );
      }
    }

    std::fill( values.begin(), values.end(), RangeType(0));

    //lf.evaluateQuadrature( quad, val );
    for( size_t qp = 0; qp < quad.nop(); ++qp )
    {
      auto& baseFct = basisFcts[ qp ];
      int dof = 0;
      for( int b=0; b<singleBase; ++b )
      {
        for( int d=0; d<RangeType::dimension; ++d, ++dof  )
        {
          for( int l=0; l<lanes; ++l )
          {
            // TODO: use mul_add
            values[ l*nop + qp ][ d ] += lfSimd[ l*nBaseFcts + dof ] * baseFct[ b ];
          }
        }
      }
    }

    // values * weights
    const int quadNop = quad.nop();
    {
      for( int qp = 0; qp < quadNop*lanes; ++qp )
      {
        //values[ qp ][ d ][ i ] *= weights[ qp ][ d ][ i ];
        values[ qp ] *= quad.weight( qp );
      }
    }

    // axpy operation
    {
      const int baseSize = singleBase * RangeType :: dimension;

      for( int b=0; b<baseSize; ++b )
      {
        dofs[ b ] = 0.0;
      }

      for( int qp = 0; qp< quadNop; ++qp )
      {
        auto& baseFct = basisFcts[ qp ];
        //const auto& val = values[ qp ];
        int dof = 0;
        for( int b=0; b<singleBase; ++b )
        {
          for( int d=0; d<RangeType::dimension; ++d, ++dof )
          {
            for( int l=0; l<lanes; ++l )
            {
              dofs[ l*nBaseFcts + dof ] += values[ l*nop + qp ][ d ] * baseFct[ b ];
              //dofs[ dof ] = mul_add( val[ d ], baseFct[ dof ][ d ], dofs[ dof ]);
            }
          }
        }
      }
    }

    //Dune::Fem::ForLoop< GetDofs, 0, lanes-1 >::apply( f, entityPtr, lfSimd );
    for( int i=0; i<lanes; ++i )
    {
      if( entityPtr[ i ] )
      {
        Dune::Fem::StaticArray< double > wrapper( nBaseFcts, lfSimd.data() + nBaseFcts * i);
        dest.setLocalDofs( *(entityPtr[ i ]), wrapper );
      }
    }
    //Dune::Fem::ForLoop< SetDofs, 0, lanes-1 >::apply( dest, entityPtr, dofs );
    //SetDofsSimd::apply( dest, entityPtr, dofs );

    if( it == endit ) break ;
    /*
    for( int i = 0; i<lanes; ++i )
    {
      if( iterator[ i ] != endIterator[ i ])
        ++iterator[ i ];
    }
    */
  }
}

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

double gFlops( const unsigned long flop, const double time )
{
  // theoretical 4 * 2.7 * 16 = 172.8 GFLOPs (43.2 per core)

  //double socket = 1;
  //double coresPerSocket = 4;
  //double cyclesPerSec = 2.7;

  return double(flop) * 1e-9 / time;
}

// ********************************************************************
double algorithm ( MyGridType &grid, DiscreteFunctionType &solution, bool display )
{
   ExactSolution f;
   // create exact solution for error evaluation
   typedef GridFunctionAdapter< ExactSolution, GridPartType >  GridFunctionType;
   GridFunctionType exactSolution( "exact solution", f, solution.gridPart(), 2 );

   // L2 error class
   Dune :: Fem :: L2Norm< GridPartType > l2norm( solution.gridPart() );

   DiscreteFunctionType copy( solution );

   Dune::Fem::Double::reset();
   //! perform l2-projection
   interpolate(exactSolution, solution);

   Dune::Timer  timer;
   //! perform l2-projection
   interpolate(solution, copy);

   double time = timer.elapsed();
   unsigned long flop = Dune::Fem::Double::reset();
   std::cout << "Standard Interpolation: " << timer.elapsed() << " flops = " << flop << std::endl;
   std::cout << "GFLOPs: " << gFlops( flop, time ) << std::endl;

   // calculation L2 error
   // pol ord for calculation the error should be higher than
   // pol for evaluation the basefunctions
   double error = l2norm.distance( exactSolution, solution );
   std::cout << "L2 Error: " << error << "\n";
   double error2 = l2norm.distance( exactSolution, copy );
   std::cout << "L2 Error(copy): " << error2 << "\n\n";


   /*
   timer.reset();
   // do l2 projection for one df to another
   customL2ProjectionNormal< 16 >( solution, copy );
   std::cout << "Fake simd 4 interpolation: " << timer.elapsed() << std::endl;
   error2 = l2norm.distance( exactSolution, copy );
   std::cout << "L2 Error: " << error2 << "\n\n";
   */


#if HAVE_VC && !defined(COUNT_FLOPS)
   timer.reset();
   customL2Projection<Vc::SimdArray<double, 24>>( solution, copy );
   std::cout << "Simd interpolation (Vc::SimdArray<double, 24>): " << timer.elapsed() << std::endl;
   error2 = l2norm.distance( exactSolution, copy );
   std::cout << "L2 Error: " << error2 << "\n\n";

   timer.reset();
   customL2Projection<Vc::SimdArray<double, 20>>( solution, copy );
   std::cout << "Simd interpolation (Vc::SimdArray<double, 20>): " << timer.elapsed() << std::endl;
   error2 = l2norm.distance( exactSolution, copy );
   std::cout << "L2 Error: " << error2 << "\n\n";

   timer.reset();
   customL2Projection<Vc::SimdArray<double, 16>>( solution, copy );
   std::cout << "Simd interpolation (Vc::SimdArray<double, 16>): " << timer.elapsed() << std::endl;
   error2 = l2norm.distance( exactSolution, copy );
   std::cout << "L2 Error: " << error2 << "\n\n";

   timer.reset();
   customL2Projection<Vc::SimdArray<double, 12>>( solution, copy );
   std::cout << "Simd interpolation (Vc::SimdArray<double, 12>): " << timer.elapsed() << std::endl;
   error2 = l2norm.distance( exactSolution, copy );
   std::cout << "L2 Error: " << error2 << "\n\n";

   timer.reset();
   customL2Projection<Vc::SimdArray<double, 8>>( solution, copy );
   std::cout << "Simd interpolation (Vc::SimdArray<double, 8>): " << timer.elapsed() << std::endl;
   error2 = l2norm.distance( exactSolution, copy );
   std::cout << "L2 Error: " << error2 << "\n\n";

   timer.reset();
   customL2Projection<Vc::SimdArray<double, 4>>( solution, copy );
   std::cout << "Simd interpolation (Vc::SimdArray<double, 4>): " << timer.elapsed() << std::endl;
   error2 = l2norm.distance( exactSolution, copy );
   std::cout << "L2 Error: " << error2 << "\n\n";
#endif

#if HAVE_DUNE_VECTORCLASS && !defined(COUNT_FLOPS)
   timer.reset();
   customL2Projection<Vec4d>( solution, copy );
   std::cout << "Simd interpolation (Vec4d): " << timer.elapsed() << std::endl;
   error2 = l2norm.distance( exactSolution, copy );
   std::cout << "L2 Error: " << error2 << "\n\n";

   timer.reset();
   customL2Projection<Vec8d>( solution, copy );
   std::cout << "Simd interpolation (Vec8d): " << timer.elapsed() << std::endl;
   error2 = l2norm.distance( exactSolution, copy );
   std::cout << "L2 Error: " << error2 << "\n\n";
#endif

   //interpolate( solution, copy );

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

  int ml = 2;
  if( argc > 1 )
    ml = atoi( argv[1] );
  if( ml == -1 )
  {
    MyGridType &grid = Dune::Fem::TestGrid::grid();
    GridPartType part ( grid );
    generateCode( part );
    return 0;
  }

  std::vector< double> error(ml);

  MyGridType &grid = Dune::Fem::TestGrid::grid();
  const int step = Dune::Fem::TestGrid::refineStepsForHalf();

  GridPartType part ( grid );
  DiscreteFunctionSpaceType linFuncSpace ( part );
  DiscreteFunctionType solution ( "sol", linFuncSpace );
  solution.clear();

  for(int i=0; i<ml; i+=step)
  {
    GlobalRefine::apply(grid,step);
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

