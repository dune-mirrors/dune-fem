#ifndef DUNE_FEM_DERITYPE_HH
#define DUNE_FEM_DERITYPE_HH

#include <dune/common/fvector.hh>

namespace Dune
{

  //! type of derivative component chooser 
  typedef int deriType;

  //! type of derivative specializer 
  template< int dim >
  struct DiffVariable
  {
    typedef FieldVector< deriType, dim > Type;
  };

}

#endif // #ifndef DUNE_FEM_DERITYPE_HH
