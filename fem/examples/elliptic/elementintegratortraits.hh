/**************************************************************************
**       Title: Example of a ElementIntegratorTraits class implementation
**    $RCSfile$
**   $Revision$$Name$
**       $Date$
**   Copyright: GPL $Author$
** Description: the file contains a default implementation of a
**              ElementIntegratorTraits class to be used with an 
**              appropriate Model in an FEOp or RhsAssembler for 
**              solving a general elliptic problem.
**
**************************************************************************/

#ifndef DUNE_ELEMENTINTEGRATORTRAITS_HH
#define DUNE_ELEMENTINTEGRATORTRAITS_HH

#include <config.h>
#include <dune/grid/io/file/dgfparser/gridtype.hh>
#include <dune/fem/quadrature/elementquadrature.hh>
#include <dune/fem/space/common/filteredgrid.hh>
#include <dune/fem/operator/matrixadapter.hh>

namespace Dune
{

/*======================================================================*/
/*!
 *  \class DefaultElementIntegratorTraits
 *  \brief The DefaultElementIntegratorTraits provides Type-Information 
 *         for the ElementMatrices and FEOp operator.
 *
 *  default implementation of a ElementIntegratorTraits class to be used 
 *  with an appropriate Model in an FEOp for solving a general elliptic 
 *  problem.
 *
 *  It is only 
 *  considered to yield information by its types, no member variables or 
 *  methods are provided, neither is it instantiated at any time.
 * 
 *  Currently scalar functions and Lagrange-Basis of degree 1 are used, 
 *  elementquadratures are chosen for any quadrature in the FEOp and 
 *  ElementIntegrators.
 *
 *  All types are derived from GridType and dimworld obtained by 
 *  inclusion of gridtype.hh
 *
 *  Essential Datatypes without explicit interface:
 *
 *  required for ElementQuadratureTypes:
 *     constructor with arguments (entity, order)
 *     nop, weight, point methods
 *
 *  required for IntersectionQuadratureTypes:
 *     enum INSIDE
 *     constructor with arguments (entity, order, INSIDE)
 *     nop, weight, point methods
 *     localpoint, geometry methods
 *
 *  required for LocalMatrixType (e.g. satisfied by 
 *  FieldMatrixAdapter<Fieldmatrix<...>>):
 *     constructor without arguments
 *     rows(), cols() methods
 *     add(rown, coln, value) allows writable additive access 
 *     to ij-th component.
 *
 *  The class can be taken as an example for own implementations.
 */
/*======================================================================*/

template< int order >
class EllipticElementIntegratorTraits
{
public:
  //! specific choices for types, based on GridType from inclusion of 
  //! gridtype.hh

  typedef  LeafGridPart<GridType> GridPartType;

  typedef FunctionSpace < double, double, dimworld, 1 > FunctionSpaceType;
  typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType,
                                         order > DiscreteFunctionSpaceType;
  typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > 
          DiscreteFunctionType;

  enum { elementMatrixSize = 1000 };  
  typedef FieldMatrixAdapter< FieldMatrix< double, elementMatrixSize,
                                           elementMatrixSize> >
    ElementMatrixType; 

  typedef ElementQuadrature< GridPartType, 0 > ElementQuadratureType; 
  typedef ElementQuadrature< GridPartType, 1 > IntersectionQuadratureType; 
  enum { quadDegree = 5 }; //<! degree of quadrature
  
  //! derived types
  typedef typename GridPartType :: IntersectionIteratorType
    IntersectionIteratorType;
  typedef typename GridType :: Codim< 0 > :: Entity EntityType;
  typedef typename GridType :: Codim< 0 > :: EntityPointer EntityPointerType;
  typedef EntityType::ctype CoordType; 
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
    BaseFunctionSetType;
  typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
    JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;

  //! default definitions, not to be changed: 
  enum { dim = GridType :: dimension };
  // do not remove any of the following bnd-type, as selection according to 
  // this is happening in the operator
  enum BoundaryType { Dirichlet, Neumann, Robin };
};
 
}; // end namespace Dune

#endif



