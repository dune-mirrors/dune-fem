#ifndef DUNE_FEM_SAPCE_EVALUATECALLERDEFAULTIMPL_HH
#define DUNE_FEM_SAPCE_EVALUATECALLERDEFAULTIMPL_HH

#include <iostream>

#include <dune/fem/space/basisfunctionset/evaluatecallerdeclaration.hh>

namespace Dune {
namespace Fem {
namespace Codegen {

/////////////////////////////////////////////////////////////////////////
//
//  evaluate and store results in a vector
//
/////////////////////////////////////////////////////////////////////////
template <class BaseFunctionSet, class Geometry, int dimR, int numRows, int numCols>
struct EvaluateRanges
{
  template< class QuadratureType,
            class RangeVectorType,
            class RangeFactorType,
            class LocalDofVectorType>
  static void eval( const QuadratureType& quad,
                    const RangeVectorType& rangeStorageTransposed,
                    const LocalDofVectorType& dofs,
                    RangeFactorType &rangeVector)
  {
    std::cout << "Default implementation of EvaluateRanges" << std::endl;
    std::abort();
  }
};

template <class BaseFunctionSet, int dimR, int numRows, int numCols>
struct EvaluateRanges<BaseFunctionSet, EmptyGeometry, dimR, numRows, numCols >
{
  static const int dimRange = dimR ;      // dimRange
  static const int quadNop  = numRows ;   // number of quadrature points
  static const int scalarBasis = numCols; // number of scalar basis functions
  static const int lanes = 4;             // simdWidth
  static const int basisLoop = (scalarBasis / lanes) * lanes;
  static const int nResult = dimRange * quadNop;
  static const int dofSkip = dimRange * lanes;

  template< class QuadratureType,
            class RangeVectorType,
            class RangeFactorType,
            class LocalDofVectorType>
  static void eval( const QuadratureType& quad,
                    const RangeVectorType& rangeStorageTransposed,
                    const LocalDofVectorType& dofs,
                    RangeFactorType &rangeVector)
  {
    typedef double Field;
    static thread_local std::vector< Field > memory( nResult );
    Field* resultTmp = memory.data();

    for( int i=0; i < nResult; ++i ) resultTmp[ i ] = 0;

    Field* result[ dimRange ];
    for( int i=0; i<dimRange; ++i )
      result[ i ] = resultTmp + i*quadNop;

    const Field* baseData = rangeStorageTransposed.data();
    if constexpr ( basisLoop > 0 )
    {
      for( int col = 0, dof = 0 ; col < basisLoop; col += lanes, dof += dofSkip )
      {
        const Field* base[ lanes ];
        for( int l=0; l<lanes; ++l )
          base[ l ] = baseData + ((col + l )*quadNop);

        evalRangesLoop( &dofs[ dof ], base, result );
      }
    }

    if constexpr ( scalarBasis > basisLoop )
    {
      // remainder iteration
      for( int col = basisLoop, dof = basisLoop*dimRange ; col < scalarBasis ; ++col )
      {
        const Field* base0 = baseData + (col * quadNop);
        // Loop over all quadrature points
        for(int qp = 0; qp < quadNop ; ++qp )
        {
          const double phi0 = base0[ qp ];
          for( int i=0; i<dimRange; ++i, ++dof )
            result[ i ][ qp ] += phi0 * dofs[ dof ] ;
        }
      }
    }

    // store result
    for(int qp = 0; qp < quadNop ; ++qp )
    {
      auto& res = rangeVector[ qp ];
      for( int i=0; i<dimRange; ++i )
        res[ i ] = result[ i ][ qp ];
    }
  }

private:
  static void evalRangesLoop(const double* dofs,
                             const double** __restrict__ basis,
                             double** __restrict__ result)
  {
    // Loop over all quadrature points
    for(int qp = 0; qp < quadNop ; ++qp )
    {
      for( int l=0, dof = 0; l<lanes; ++l )
      {
        const double base = basis[ l ][ qp ];
        for( int i=0; i<dimRange; ++i, ++dof )
          result[ i ][ qp ] += base * dofs[ dof ];
      }
    }
  }
};

////////////////////////////////////////////////////////////////////
//
//  --evaluateJacobians
//
////////////////////////////////////////////////////////////////////
template <class BaseFunctionSet, class Geometry,
          int dimRange, int numRows, int numCols>
struct EvaluateJacobians
{
  template< class QuadratureType,
            class JacobianRangeVectorType,
            class JacobianRangeFactorType,
            class LocalDofVectorType>
  static void eval( const QuadratureType& quad,
                    const Geometry& geometry,
                    const JacobianRangeVectorType& jacobianStorage,
                    const LocalDofVectorType& dofs,
                    JacobianRangeFactorType &jacFactors)
  {
    std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::evaluateJacobians< "
              << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
    std::abort();
  }
};

template <class BaseFunctionSet,
          int dimRange, int numRows, int numCols>
struct EvaluateJacobians< BaseFunctionSet, EmptyGeometry, dimRange, numRows, numCols >
{
  template< class QuadratureType,
            class JacobianRangeVectorType,
            class JacobianRangeFactorType,
            class LocalDofVectorType>
  static void eval( const QuadratureType&,
                    const EmptyGeometry&,
                    const JacobianRangeVectorType&,
                    const LocalDofVectorType&,
                    const JacobianRangeFactorType& )
  {
    std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::evaluateJacobians< "
              << "EmptyGeo, " << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
    std::abort();
  }
};

/////////////////////////////////////////////////////////////
//
//  --axpyRanges -- add a vector of ranges to the dof vector
//
/////////////////////////////////////////////////////////////
template <class BaseFunctionSet, class Geometry,
          int dimRange, int numRows, int numCols>
struct AxpyRanges
{
  template< class QuadratureType,
            class RangeVectorType,
            class RangeFactorType,
            class LocalDofVectorType>
  static void axpy( const QuadratureType& quad,
                    const RangeVectorType& rangeStorage,
                    const RangeFactorType &rangeFactors,
                    LocalDofVectorType& dofs)
  {
    std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::axpyRanges <"
              << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
    std::abort();
  }
};

template <class BaseFunctionSet,
          int dimRange, int numRows, int numCols>
struct AxpyRanges<BaseFunctionSet, EmptyGeometry, dimRange, numRows, numCols>
{
  template< class QuadratureType,
            class RangeVectorType,
            class RangeFactorType,
            class LocalDofVectorType>
  static void axpy( const QuadratureType& quad,
                    const RangeVectorType& rangeStorage,
                    const RangeFactorType &rangeFactors,
                    LocalDofVectorType& dofs)
  {
    std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::axpyRanges <"
              << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
    std::abort();
  }
};


///////////////////////////////////////////////////////////
//  applyAxpy Jacobian
///////////////////////////////////////////////////////////
template <class BaseFunctionSet, class Geometry,
          int dimRange, int numRows, int numCols>
struct AxpyJacobians
{
  template< class QuadratureType,
            class JacobianRangeVectorType,
            class JacobianRangeFactorType,
            class LocalDofVectorType>
  static void axpy( const QuadratureType& quad,
                    const Geometry& geometry,
                    const JacobianRangeVectorType& jacobianStorage,
                    const JacobianRangeFactorType& jacFactors,
                    LocalDofVectorType& dofs)
  {
    std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::axpyJacobian <"
              << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
    std::abort();
  }
};

template <class BaseFunctionSet,
          int dimRange, int numRows, int numCols>
struct AxpyJacobians< BaseFunctionSet, EmptyGeometry, dimRange, numRows, numCols >
{
  template< class QuadratureType,
            class JacobianRangeVectorType,
            class JacobianRangeFactorType,
            class LocalDofVectorType>
  static void axpy( const QuadratureType&,
                    const EmptyGeometry&,
                    const JacobianRangeVectorType&,
                    const JacobianRangeFactorType &,
                    LocalDofVectorType&)
  {
    std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::axpyJacobians" << std::endl;
    std::abort();
  }
};

}}}

#endif
