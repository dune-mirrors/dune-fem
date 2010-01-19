#include <config.h>

#include <iostream>

#include <sys/stat.h>  
#include <dirent.h>

#include <dune/fem/io/io.hh>

namespace Dune
{

  bool createDirectory ( const std::string &name )
  {
    // try to open directory (returns null pointer, if path does not exist)
    DIR *dir = opendir( name.c_str() );
    if( dir != 0 )
    {
      if( closedir( dir ) < 0 )
        std::cerr << "Error: Could not close directory." << std::endl;
      return true;
    }

    // try to create the father directory
    size_t pos = name.rfind( '/' );
    if( pos != std::string::npos )
    {
      const std::string father = name.substr( 0, pos );
      if( !createDirectory( father ) )
        return false;
    }

    // try to create the new child directory
    mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
    return (mkdir( name.c_str(), mode ) >= 0);
  }

}
