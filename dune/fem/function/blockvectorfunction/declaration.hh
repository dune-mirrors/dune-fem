#ifndef DUNE_FEM_BLOCKVECTORFUNCTION_DECLARATION_HH
#define DUNE_FEM_BLOCKVECTORFUNCTION_DECLARATION_HH

#include <dune/common/fvector.hh>

namespace Dune
{
  namespace Fem
  {
    //**********************************************************************
    //! @ingroup BlockVectorDFunction
    //  --ISTLBlockVectorDiscreteFunction
    //
    //! this is one special implementation of a discrete function using an
    //! array for storing the dofs.
    //! \note The class Block needs to fullfil the Dune::DenseVector interface.
    //!
    //**********************************************************************
    template <class DiscreteFunctionSpace,
              class Block = Dune::FieldVector< typename DiscreteFunctionSpace::RangeFieldType, DiscreteFunctionSpace::localBlockSize > >
    class ISTLBlockVectorDiscreteFunction;
  }
}
#endif
