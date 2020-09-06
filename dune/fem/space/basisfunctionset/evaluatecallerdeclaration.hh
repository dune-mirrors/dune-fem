#ifndef DUNE_FEM_EVALUATECALLER_DECLARATION_HH
#define DUNE_FEM_EVALUATECALLER_DECLARATION_HH

#include <cstdlib>
#include <iostream>

namespace Dune
{
  namespace Fem {

    namespace Codegen {

      // empty class for specialization of evaluation classes in basefunctionsets.hh
      class EmptyGeometry {};

      template <class BaseFunctionSet, class Geometry, int dimRange, int numRows, int numCols>
      struct EvaluateRanges;

      template <class BaseFunctionSet, class Geometry,
                int dimRange, int numRows, int numCols>
      struct EvaluateJacobians;

      template <class BaseFunctionSet, class Geometry,
                int dimRange, int numRows, int numCols>
      struct AxpyRanges;

      template <class BaseFunctionSet, class Geometry,
                int dimRange, int numRows, int numCols>
      struct AxpyJacobians;
    }
  }
}
#endif // DUNE_FEM_EVALUATECALLER_DECLARATION_HH
