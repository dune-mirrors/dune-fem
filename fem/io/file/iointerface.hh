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

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

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

  //! \brief save structured macro grid 
  //! \param macroFileName is the macro which should be saved in DGF format  
  virtual void saveMacroGrid(const std::string macroFileName) const = 0;
  
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

  //! if grid is structured grid, write macro file 
  template <class GridImp>
  static void writeMacroGrid(const GridImp& grid, 
                             const std::string& macroname,
                             const std::string& path, 
                             const std::string& prefix) 
  {
    // do nothing for unstructured grids 
    if( Capabilities::IsUnstructured<GridType>::v ) return;
    
    // create file descriptor 
    std::ifstream gridin(macroname.c_str());
    if( !gridin) 
    {
      std::cerr << "Couldn't open file `" << macroname << "' ! \n";
      return ;
    } 
        
    // read interval information of structured grid 
    IntervalBlock interval(gridin);
    if(!interval.isactive()) 
    {
      std::cerr<<"Did not find IntervalBlock in macro grid file `" << macroname << "' ! \n";
      return;
    }
    
    std::string filename(path);
    filename += "/g";
    filename += prefix;
    filename += ".macro";

    enum { dimworld = GridImp :: dimensionworld };

    saveMacroGridImp<dimworld> (interval,filename); 
    return;
  }

protected:
  template <int dimworld> 
  static void saveMacroGridImp(IntervalBlock& interval, std::string filename) 
  {
    FieldVector<double,dimworld> lang;
    FieldVector<int,dimworld>    anz;
    FieldVector<double,dimworld>  h;
    
    // set values 
    for (int i=0;i<dimworld;i++) 
    {
      lang[i] = interval.length(i);
      anz[i]  = interval.segments(i);
      h[i] = lang[i]/anz[i];
    }

#if HAVE_MPI 
    // write sub grid 
    {
      typedef FieldVector<int,dimworld> iTupel;

      // origin is zero 
      iTupel o(0);

      iTupel o_interior;
      iTupel s_interior;

      enum { tag = MultiYGrid<dimworld,double> ::tag };
      Torus<dimworld> torus(MPI_COMM_WORLD,tag,anz);
      torus.partition( torus.rank() , o,anz,
                       o_interior,s_interior);

      FieldVector<double,dimworld> origin(0.0);
      FieldVector<double,dimworld> sublang(0.0);
      for(int i=0; i<dimworld; ++i)
      {
        origin[i] = o_interior[i] * h[i];
        sublang[i] = origin[i] + (s_interior[i] * h[i]);
      }

      writeStructuredGrid(filename,origin,sublang,s_interior);
    }
#endif
    {
      // write global file for recovery 
      filename += ".global";
      FieldVector<double,dimworld> zero(0.0);
      writeStructuredGrid(filename,zero,lang,anz);
    }
  }

  //! write structured grid as DGF file 
  template <int dimworld>
  static void writeStructuredGrid(const std::string& filename,
                           const FieldVector<double,dimworld>& origin,
                           const FieldVector<double,dimworld>& lang,
                           const FieldVector<int,dimworld>& anz)
  {
    std::ofstream file (filename.c_str());
    if( file.is_open())
    {
      file << "DGF" << std::endl;
      file << "Interval" << std::endl;
      // write first point 
      for(int i=0;i<dimworld; ++i)
      {
        file << origin[i] << " ";
      }
      file << std::endl;
      // write second point 
      for(int i=0;i<dimworld; ++i)
      {
        file << lang[i] << " ";
      }
      file << std::endl;
      // write number of intervals in each direction 
      for(int i=0;i<dimworld; ++i)
      {
        file << anz[i] << " ";
      }
      file << std::endl;
      file << "#" << std::endl;

      file << "BoundaryDomain" << std::endl;
      file << "default 1" << std::endl;
      file << "#" << std::endl;
    }
    else
    {
      std::cerr << "Couldn't open file `" << filename << "' !\n";
    }
  }
}; // end class IOInterface 

} // end namespace Dune  
#endif
