#include <config.h>

#include <fstream>
#include <iostream>
#include <cstdio>

#include <sys/stat.h>
#include <dirent.h>

#include <dune/fem/io/io.hh>
#include <dune/common/exceptions.hh>

namespace Dune
{

  namespace Fem
  {

    bool createDirectory ( const std::string &inName )
    {
      std::string name = inName;

      // strip of last character if it is a '/'
      if( name[ name.size() - 1 ] == '/' )
        name = name.substr( 0, name.size() -1 );

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


    bool fileExists ( const std::string &name )
    {
      std::ifstream file( name.c_str() );
      return file.is_open();
    }


    bool directoryExists ( const std::string &name )
    {
      // if directory does not exist return false
      DIR* directory = opendir( name.c_str() );
      const bool directoryExists = (directory != 0);
      // close directory again
      if ( directoryExists )
        closedir( directory );
      return directoryExists;
    }


    std::string executeCommand ( const std::string &command )
    {
      std::string returnString;
      FILE *pipe = popen( command.c_str(), "r" );

      if( !pipe )
        DUNE_THROW( IOError, "Unable to execute '" << command << "'." );

      std::size_t size;
      do
      {
        char buffer[ 4096 ];
        size = fread( buffer, sizeof( char ), sizeof( buffer ) / sizeof( char ), pipe );
        returnString.append( buffer, size );
      }
      while( size > std::size_t( 0 ) );

      const int status = pclose( pipe );
      if( status != 0 )
        DUNE_THROW( IOError, "Command '" << command << "' returned unsuccessfully (" << status << ")." );

      return returnString;
    }

  } // namespace Fem

} // namespace Dune
