// (c) Robert Kloefkorn
#ifndef DUNE_INPUTOUPUTLOCK_HH
#define DUNE_INPUTOUPUTLOCK_HH

//- system includes 
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

namespace Dune {

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
};

// create lock file 
inline FileIOLock :: FileIOLock (const std::string& fn) 
  : filename_(fn)  
{
  if( filename_ == "" )  
  {  
    filename_ = "lock";  
  }
  else   {
    filename_ += ".lock";
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

} // end namespace Dune 
#endif
