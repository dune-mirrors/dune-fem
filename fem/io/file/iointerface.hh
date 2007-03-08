#ifndef DUNE_IOINTERFACE_HH
#define DUNE_IOINTERFACE_HH

//- system includes 
#include <sys/stat.h>  
#include <sys/types.h>
#include <dirent.h>

//- Dune includes 
// defines function readParameter 
#include <dune/fem/io/file/asciiparser.hh>

// grape data io 
#include <dune/fem/io/file/grapedataio.hh>
// input and output of tuples 
#include <dune/fem/io/file/grapetuple.hh>

// if grape was configured then include headers 
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

namespace Dune {
  
/** 
 \brief IOInterface to write data to hard disk 
*/ 
class IOInterface {

protected: 
  //! default constructor 
  IOInterface() {}
  
public: 
  //! destructor 
  virtual ~IOInterface () {}

  //! write data to disc
  //! \param time actual time of computation 
  //! \param timestep current number of time step 
  //! \param balancecounter current number in load balance cycle 
  virtual void write(double time, int timestep, int balancecounter) const = 0; 

  //! display data if HAVE_GRAPE is 1  
  virtual void display() const = 0; 
  
  //! standard path reading and creation method 
  //! rank is added to output path 
  static std::string readPath(const std::string& paramfile, int rank)
  {
    std::string path;
    // read output path from parameter file 
    if( readParameter(paramfile,"OutputPath",path) ) 
    {
      // add proc number to path 
      {
        path += "_";
        std::stringstream rankDummy;
        rankDummy << rank; 
        path += rankDummy.str();
      }
      
      {
        // try to open directory 
        DIR * dir = opendir(path.c_str());
        // if path does not exists, null pointer is returned
        if( !dir )
        {
          // see <sys/stat.h> for definition 
          // this stands for mode 755
          mode_t mode = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
          // try to create dicretory
          // mkdir returns int < 0 if creation failed 
          if( mkdir( path.c_str(), mode ) < 0) 
          {
            std::cerr << "Failed to create path `" << path << "' !" << std::endl;
          }
        }
        else 
        {
          // close dir, if fails number < 0 is returned 
          if( closedir(dir) < 0 )
          {
            std::cerr << "Couldn't close output path! " << std::endl;
          }
        }
      }
      return path;
    }
    else 
    {
      std::cerr << "Couldn't read output path, exiting... " << std::endl;
    }
    // returning empty path 
    return path;
  }

}; // end class IOInterface 

} // end namespace Dune  
#endif
