#ifndef DUNE_FEM_DERITYPE_HH
#define DUNE_FEM_DERITYPE_HH

#warning "Deprecated header, use int instead of deriType"

#include <dune/common/fvector.hh>

namespace Dune
{

  //! type of derivative component chooser 
  typedef int deriType;

  //! type of derivative specializer 
  template< int dim >
  struct DiffVariable
  {
    typedef FieldVector< int, dim > Type;
  };

}

#endif // #ifndef DUNE_FEM_DERITYPE_HH
