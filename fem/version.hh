#ifndef DUNE_FEM_VERSION_HH
#define DUNE_FEM_VERSION_HH

// These two lines define the dune-fem version.
#define DUNE_FEM_VERSION_MAJOR 0
#define DUNE_FEM_VERSION_MINOR 10
#define DUNE_FEM_VERSION_REVISION 0

#include <sstream>

namespace Dune
{
  
  class DuneFEM
  {
  private:
    typedef DuneFEM ThisType;

  public:
    enum { MajorVersion = DUNE_FEM_VERSION_MAJOR };
    enum { MinorVersion = DUNE_FEM_VERSION_MINOR };
    enum { Revision = DUNE_FEM_VERSION_REVISION };

  private:
    DuneFEM ();
    DuneFEM ( const ThisType &other );

  public:
    static inline std :: string version ()
    {
      std :: ostringstream s;
      s << "dune-fem " << MajorVersion << "." << MinorVersion
        << "." << Revision;
      return s.str();
    }
    
  };
  
};

#endif
