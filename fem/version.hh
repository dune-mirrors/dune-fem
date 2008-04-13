#ifndef DUNE_FEM_VERSION_HH
#define DUNE_FEM_VERSION_HH

// These two lines define the dune-fem version.
#define DUNE_FEM_VERSION_MAJOR 0
#define DUNE_FEM_VERSION_MINOR 9
#define DUNE_FEM_VERSION_REVISION 2

#include <sstream>

namespace Dune
{
  
  class DuneFEM
  {
  private:
    typedef DuneFEM ThisType;

  public:
    static const unsigned int MajorVersion = DUNE_FEM_VERSION_MAJOR;
    static const unsigned int MinorVersion = DUNE_FEM_VERSION_MINOR;
    static const unsigned int Revision = DUNE_FEM_VERSION_REVISION;

  private:
    DuneFEM ();
    DuneFEM ( const ThisType &other );

  public:
    static inline std :: string version ()
    {
      std :: ostringstream s;
      s << "dune-fem " << MajorVersion
        << "." << MinorVersion
        << "." << Revision;
      return s.str();
    }

    static inline std :: string version ( unsigned int versionId )
    {
      std :: ostringstream s;
      s << "dune-fem " << (versionId >> 24)
        << "." << ((versionId >> 16) & 0xff)
        << "." << (versionId & 0xffff);
      return s.str();
    }

    static inline unsigned int versionId ( unsigned int majorVersion,
                                           unsigned int minorVersion,
                                           unsigned int revision )
    {
      return (majorVersion << 24) + (minorVersion << 16) + revision;
    }

    static inline unsigned int versionId ()
    {
      return versionId( MajorVersion, MinorVersion, Revision );
    }
  };
  
};

#endif
