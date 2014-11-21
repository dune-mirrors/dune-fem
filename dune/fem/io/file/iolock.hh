// (c) Robert Kloefkorn
#ifndef DUNE_FEM_INPUTOUPUTLOCK_HH
#define DUNE_FEM_INPUTOUPUTLOCK_HH

//- system includes
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <string>

namespace Dune
{

  namespace Fem
  {

    //! creates and removes lock file during backup process
    class FileIOLock
    {
      std::string filename_;
      FileIOLock ( const FileIOLock & );
    public :
      //! creates lock file
      FileIOLock(const std::string& fn);
      //! removes lock file
     ~FileIOLock() ;

      //! suffix that is appended to lock files
      static const char * suffix() { return "lock"; }
    };

    //! check if lock file exists and aborts if true
    class FileIOCheckError
    {
      FileIOCheckError( const FileIOCheckError & );
    public :
      //! creates lock file
      FileIOCheckError(const std::string& fn)
      {
        std::string lockfile(fn);
        lockfile += ".";
        lockfile += FileIOLock::suffix();

        std::ifstream file ( lockfile.c_str () );
        if( file.is_open() )
        {
          std::cerr << "ERROR: data set `"<<fn<<"' not complete, lock file exists! " << std::endl;
          abort();
        }
      }
    };

    ////////////////////////////////////////////////////////////////
    //
    //  INLINE
    //
    ////////////////////////////////////////////////////////////////

    // create lock file
    inline FileIOLock :: FileIOLock (const std::string& fn)
      : filename_(fn)
    {
      if( filename_ == "" )
      {
        filename_ = suffix();
      }
      else   {
        filename_ += ".";
        filename_ += suffix();
      }

      std::ofstream file ( filename_.c_str() );
      if( !file )
      {
        std::cerr << "WARNING: Couldn't open lock file `"<<filename_<<"' in: ";
        std::cerr << __FILE__<< " line: "<< __LINE__ << std::endl;
      }
      else
      {
        file.close();
      }
      return ;
    }

    // remove lock file
    inline FileIOLock :: ~FileIOLock ()
    {
      if (filename_ != "")
      {
        int test = remove (filename_.c_str()) ;
        if (test != 0)
        {
          std::cerr << "WARNING: Couldn't remove lock file `"<<filename_<<"' in: ";
          std::cerr <<__FILE__<<" line: " <<__LINE__ << std::endl ;
        }
      }
      return ;
    }

  } // namespace Fem

} // namespace Dune
#endif // #ifndef DUNE_FEM_INPUTOUPUTLOCK_HH
