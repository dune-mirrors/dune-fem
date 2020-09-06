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
template <class BaseFunctionSet, class Geometry, int dimRange, int numRows, int numCols>
struct EvaluateRanges
{
  template< class QuadratureType,
            class RangeVectorType,
            class LocalDofVectorType,
            class RangeFactorType>
  static void eval( const QuadratureType& quad,
                    const RangeVectorType& rangeStorage,
                    const LocalDofVectorType& dofs,
                    RangeFactorType &rangeFactors )
  {
    std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::evaluateRanges< "
              << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
    std::abort();
  }
};

template <class BaseFunctionSet, int dimRange, int numRows, int numCols>
struct EvaluateRanges<BaseFunctionSet, EmptyGeometry, dimRange, numRows, numCols >
{
  template< class QuadratureType,
            class RangeVectorType,
            class LocalDofVectorType,
            class RangeFactorType>
  static void eval( const QuadratureType& quad,
                    const RangeVectorType& rangeStorage,
                    const LocalDofVectorType& dofs,
                    RangeFactorType &rangeFactors)
  {
    std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::evaluateRanges< "
              << "EmptyGeo, " << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
    std::abort();
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
