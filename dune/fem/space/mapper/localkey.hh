#ifndef DUNE_FEM_SPACE_DOFMAPPER_LOCALKEY_HH
#define DUNE_FEM_SPACE_DOFMAPPER_LOCALKEY_HH

#if HAVE_DUNE_LOCALFUNCTIONS
#include <dune/localfunctions/common/localkey.hh>
#endif // #if HAVE_DUNE_LOCALFUNCTIONS

namespace Dune
{

  namespace Fem
  {

#if HAVE_DUNE_LOCALFUNCTIONS

    using Dune::LocalKey;

#else // #if HAVE_DUNE_LOCALFUNCTIONS

    struct LocalKey
    {
      LocalKey ( unsigned int subEntity, unsigned int codim, unsigned int index )
      : subEntity_( subEntity ), codim_( codim ), index_( index )
      {}

      unsigned int subEntity () const { return subEntity_; }
      unsigned int codim () const { return codim_; }
      unsigned int index () const { return index_; }

    private:
      unsigned int subEntity_, codim_, index_;
    };

#endif // #else // #if HAVE_DUNE_LOCALFUNCTIONS

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DOFMAPPER_LOCALKEY_HH
