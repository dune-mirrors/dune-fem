#ifndef DUNE_GRAPEDATAIO_HH
#define DUNE_GRAPEDATAIO_HH

//- system includes 

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/utility/grapedataioformattypes.hh>

#include <dune/fem/space/common/dofmanager.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh> 

//- Local includes 
#include "asciiparser.hh"

namespace Dune {

template <int dim, int dimworld, class GridImp, bool hasBackupRestore> 
class GrapeDataIOImp 
{
public:  
  typedef GridImp GridType;
   /** Write Grid with GridType file filename and time 
   *
   * This method uses the Grid Interface Method writeGrid 
   * to actually write the grid, within this method the real file name is
   * generated out of filename and timestep 
   */
  inline static bool writeGrid (const GridType & grid, 
    const GrapeIOFileFormatType ftype, const GrapeIOStringType & fnprefix 
      , double time=0.0, int timestep=0, int precision = 6);

  //! get Grid from file with time and timestep , return true if ok 
  inline static bool readGrid (GridType & grid, 
      const GrapeIOStringType & fnprefix , double & time , int timestep);
  
  //! get Grid from file with time and timestep , return grid pointer if ok 
  inline static GridType * restoreGrid (
      const GrapeIOStringType & fnprefix , double & time , int timestep)
  {
    GridType * grid = new GridType (); 
    assert( grid );  
    readGrid(*grid,fnprefix,time,timestep);
    return grid;
  }
};

template <int dim, int dimworld, class GridImp> 
class GrapeDataIOImp<dim,dimworld,GridImp,false>
{
protected:  
  struct HackIt : public MacroGrid
  {
    typedef GridImp GridType;
    typedef FieldVector<typename
      GridType::ctype,GridType::dimensionworld> DomainType;
    public: 
    // grid auto pointer
    mutable std::auto_ptr<GridType> gridptr_;
    std::vector<double> emptyParam;
    // element and vertex parameters
    std::vector<std::pair<DomainType, std::vector<double> > > elParam,vtxParam;
    int nofElParam_,nofVtxParam_;
  };

  //! derive from GridPtr to implement release method 
  //! which is needed below
  class ImprovedGridPtr : public GridPtr<GridImp> 
  { 
    // type of base class 
    typedef GridPtr<GridImp> BaseType;
    // cast this class to HackIt to have access to internal data
    HackIt& asHack() { return *((HackIt *) this); }
  public:
    //! constructor creating macro grid 
    ImprovedGridPtr(const std::string filename) 
      : BaseType(filename)
    {
    }
    
    //! release pointer from auto pointer 
    GridImp* release () { return asHack().gridptr_.release(); } 
  };
  
  typedef GridImp GridType;
public:  
   /** Write structurd grid to file filename 
    NOTE: the macro grid file of this structured grid 
          has to stored in a file named "filename.macro" and
          the format must be Dune DGF. 
   */
  inline static bool writeGrid (const GridType & grid, 
    const GrapeIOFileFormatType ftype, const GrapeIOStringType & fnprefix 
      , double time=0.0, int timestep=0, int precision = 6)
  {
    // write dof manager, that corresponds to grid 
    bool hasDm = false;
    {
      typedef DofManager<GridImp> DofManagerType; 
      typedef DofManagerFactory<DofManagerType> DMFactoryType; 

      std::string dmname(fnprefix);
      dmname += "_dm";
      hasDm = DMFactoryType::writeDofManager(grid,dmname,timestep);
    }
   
    // write Grid itself 
    {
      std::ofstream file (fnprefix.c_str());
      if( file.is_open() )
      {
        file << "Grid: "   << grid.name() << std::endl;
        file << "Format: " << ftype <<  std::endl;
        file << "Precision: " << precision << std::endl;
        int writeDm = (hasDm)? 1 : 0;
        file << "DofManager: " << writeDm << std::endl; 
        file.close();
      }
      else {
        std::cerr << "Couldn't open file `" << fnprefix << "' ! \n";
        return false;
      }
    }
      
    // write max level and current time to grid specific file 
    { 
      const char * path = "";
      GrapeIOStringType fnstr = genFilename(path,fnprefix,timestep,precision);
      std::ofstream gridfile (fnstr.c_str());
      
      if( gridfile.is_open() )
      {
        gridfile << "MaxLevel: " << grid.maxLevel() << std::endl;
        gridfile << "Time: " << std::scientific << time << std::endl;
        gridfile.close();
      }
      else 
      {
        std::cerr << "Couldn't open file `" << fnstr << "' ! \n";
        return false;
      }
    }
    return true;
  }

  //! get Grid from file with time and timestep , return true if ok 
  inline static bool readGrid (GridType & grid, 
      const GrapeIOStringType & fnprefix , double & time , int timestep)
  {
    std::string gridname;

    bool readGridName = readParameter(fnprefix,"Grid",gridname);
    if(! readGridName ) 
    {
      std::cerr << "ERROR: Couldn't open file '"<<fnprefix<<"' !" << std::endl;
      abort();
    }

    if(grid.name() != gridname)
    {
      std::cerr << "\nERROR: '" << grid.name() << "' tries to read '" << gridname << "' file. \n";
      abort();
    }

    int precision = 6;
    readParameter(fnprefix,"Precision",precision);

    int hasDm = 0;
    readParameter(fnprefix,"DofManager",hasDm);

    {
      const char * path = "";
      GrapeIOStringType fnstr = genFilename(path,fnprefix,timestep,precision);
      
      {
        int maxLevel = 0;
        readParameter(fnstr,"MaxLevel",maxLevel);
        // refine grid 
        grid.globalRefine(maxLevel);
      }
      readParameter(fnstr,"Time",time);
    }
    return true;
  }

  //! get Grid from file with time and timestep , return true if ok 
  inline static GridType * restoreGrid (
      const GrapeIOStringType & fnprefix , double & time , int timestep)
  {
    std::string macroName (fnprefix);
    macroName += ".macro";

    GridType * grid = 0 ;
    {
      // create macro grid 
      ImprovedGridPtr gridptr(macroName);
      // release pointer 
      grid = gridptr.release();
    }
    assert( grid );
    readGrid(*grid,fnprefix,time,timestep);
    return grid;
  }
};

template <class GridImp>
class GrapeDataIO 
{
public:
   typedef GridImp GridType;

   GrapeDataIO () {};

   /** Write Grid with GridType file filename and time 
   *
   * This method uses the Grid Interface Method writeGrid 
   * to actually write the grid, within this method the real file name is
   * generated out of filename and timestep 
   */
  inline bool writeGrid (const GridType & grid, 
    const GrapeIOFileFormatType ftype, const GrapeIOStringType fnprefix 
      , double time=0.0, int timestep=0, int precision = 6) const
  {
    const bool hasBackupRestore = Capabilities::hasBackupRestoreFacilities<GridType>::v;
    return GrapeDataIOImp<GridType::dimension,GridType::dimensionworld,GridType,hasBackupRestore>::
      writeGrid(grid,ftype,fnprefix,time,timestep,precision);
  }

  //! get Grid from file with time and timestep , return true if ok 
  inline bool readGrid (GridType & grid, 
      const GrapeIOStringType fnprefix , double & time , int timestep)
  {
    const bool hasBackupRestore = Capabilities::hasBackupRestoreFacilities<GridType>::v;
    return GrapeDataIOImp<GridType::dimension,GridType::dimensionworld,GridType,hasBackupRestore>::
      readGrid(grid,fnprefix,time,timestep);
  }

  //! get Grid from file with time and timestep , return true if ok 
  inline GridType * restoreGrid(const GrapeIOStringType fnprefix , double & time , int timestep)
  {
    const bool hasBackupRestore = Capabilities::hasBackupRestoreFacilities<GridType>::v;
    return GrapeDataIOImp<GridType::dimension,GridType::dimensionworld,GridType,hasBackupRestore>::
          restoreGrid(fnprefix,time,timestep);
  }

  /**
    Write DiscreteFunctions  
  */ 
  //! write disc func information file and write dofs to file+timestep
  //! this method use the write method of the implementation of the
  //! discrete function
  template <class DiscreteFunctionType>
  inline bool writeData(const DiscreteFunctionType & df,
     const GrapeIOFileFormatType ftype, const GrapeIOStringType filename, 
      int timestep, int precision = 6);

  //! same as write only read
  template <class DiscreteFunctionType>
  inline bool readData(DiscreteFunctionType & df,
        const GrapeIOStringType filename, int timestep);
};

template <int dim, int dimworld, class GridImp, bool hasBackupRestore>
inline bool GrapeDataIOImp<dim,dimworld,GridImp,hasBackupRestore> :: writeGrid 
(const GridImp & grid,
  const GrapeIOFileFormatType ftype, const GrapeIOStringType & fnprefix , 
  double time, int timestep, int precision )
{
  bool hasDm = false;
  // write dof manager, that corresponds to grid 
  {
    typedef DofManager<GridImp> DofManagerType; 
    typedef DofManagerFactory<DofManagerType> DMFactoryType; 

    std::string dmname(fnprefix);
    dmname += "_dm";
    hasDm = DMFactoryType::writeDofManager(grid,dmname,timestep);
  }
 
  // write Grid itself 
  {
    const char *path = "";
    std::fstream file (fnprefix.c_str(),std::ios::out);
    if( file.is_open() )
    {
      file << "Grid: "   << grid.name() << std::endl;
      file << "Format: " << ftype <<  std::endl;
      file << "Precision: " << precision << std::endl;
      int writeDm = (hasDm)? 1 : 0;
      file << "DofManager: " << writeDm << std::endl; 

      GrapeIOStringType fnstr = genFilename(path,fnprefix,timestep,precision);
      
      file.close();
      switch (ftype)
      {
        case xdr  :   return grid.template writeGrid<xdr>  (fnstr,time);
        case ascii:   return grid.template writeGrid<ascii>(fnstr,time);
        default:
            {
              std::cerr << ftype << " GrapeIOFileFormatType not supported at the moment! " << __FILE__ << __LINE__ << "\n";
              assert(false);
              abort();
              return false;
            }
      }
    }
    else 
    {
      std::cerr << "Couldn't open file `" << fnprefix << "' ! \n";
      return false;
    }
  }
}

template <int dim, int dimworld, class GridImp, bool hasBackupRestore>
inline bool GrapeDataIOImp<dim,dimworld,GridImp,hasBackupRestore> :: readGrid 
(GridImp & grid, const GrapeIOStringType & fnprefix , double & time , int timestep)
{
  int helpType;

  std::string gridname;

  bool readGridName = readParameter(fnprefix,"Grid",gridname);
  if(! readGridName ) 
  {
    std::cerr << "ERROR: Couldn't open file '"<<fnprefix<<"' !" << std::endl;
    abort();
  }

  if(grid.name() != gridname)
  {
    std::cerr << "\nERROR: '" << grid.name() << "' tries to read '" << gridname << "' file. \n";
    abort();
  }

  readParameter(fnprefix,"Format",helpType);
  GrapeIOFileFormatType ftype = (GrapeIOFileFormatType) helpType;

  int precision = 6;
  readParameter(fnprefix,"Precision",precision);

  int hasDm = 0;
  readParameter(fnprefix,"DofManager",hasDm);

  const char *path = "";
  GrapeIOStringType fn = genFilename(path,fnprefix,timestep,precision);
  std::cout << "Read file: fnprefix = `" << fn << "' \n";

  bool succeded = false;
  switch (ftype)
  {
    case xdr  :   succeded = grid.template readGrid<xdr>  (fn,time); break;
    case ascii:   succeded = grid.template readGrid<ascii>(fn,time); break;
    default:
        {
          std::cerr << ftype << " GrapeIOFileFormatType not supported at the moment! \n";
          assert(false);
          abort();
          return false;
        }
  }
 
  // write dof manager, that corresponds to grid 
  if(hasDm)
  {
    typedef DofManager<GridImp> DofManagerType; 
    typedef DofManagerFactory<DofManagerType> DMFactoryType; 
    
    std::string dmname(fn);
    dmname += "_dm";
    //std::cout << "Read DofManager from file " << dmname << "\n";
    // this call creates DofManager if not already existing 
    DMFactoryType::getDofManager(grid);
    succeded = DMFactoryType::writeDofManager(grid,dmname,timestep);
  }
  return succeded;
}

template <class GridType>
template <class DiscreteFunctionType> 
inline bool GrapeDataIO<GridType> :: writeData(const DiscreteFunctionType & df, 
const GrapeIOFileFormatType ftype, const GrapeIOStringType filename, int timestep, int  precision )
{
  {
    typedef typename DiscreteFunctionType::FunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

    enum { n = DiscreteFunctionSpaceType::DimDomain };
    enum { m = DiscreteFunctionSpaceType::DimRange };

    std::ofstream file( filename.c_str() );

    if( file.is_open())
    {
      GrapeIOStringType d = typeIdentifier<DomainFieldType>();
      GrapeIOStringType r = typeIdentifier<RangeFieldType>();

      file << "DomainField: " << d << std::endl;
      file << "RangeField: " << r << std::endl;
      file << "Dim_Domain: " << n << std::endl;
      file << "Dim_Range: " << m << std::endl;
      file << "Space: " << df.space().type() << std::endl;
      file << "Format: " << ftype << std::endl;
      file << "Precision: " << precision << std::endl;
      file << "Polynom_order: " << df.space().order() << std::endl;
      file << "DataBase: " << df.name() << std::endl;
      file.close();
    }
    else 
    {
      std::cerr << "Couldn't open file `" << filename << "' ! \n";
      return false;
    }
  }

  const char * path = "";
  GrapeIOStringType fn = genFilename(path,filename,timestep,precision);

  if(ftype == xdr)
    return df.write_xdr(fn);
  if(ftype == ascii)
    return df.write_ascii(fn);
  if(ftype == pgm)
    return df.write_pgm(fn);

  return false;
}

template <class GridType>
template <class DiscreteFunctionType> 
inline bool GrapeDataIO<GridType> :: 
readData(DiscreteFunctionType & df, const GrapeIOStringType filename, int timestep)
{
    typedef typename DiscreteFunctionType::FunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

    enum { tn = DiscreteFunctionSpaceType::DimDomain };
    enum { tm = DiscreteFunctionSpaceType::DimRange };

    int n,m;
    GrapeIOStringType r,d;
    GrapeIOStringType tr (typeIdentifier<RangeFieldType>());
    GrapeIOStringType td (typeIdentifier<DomainFieldType>());

    readParameter(filename,"DomainField",d,false);
    readParameter(filename,"RangeField",r,false);
    readParameter(filename,"Dim_Domain",n,false);
    readParameter(filename,"Dim_Range",m,false);
    int space;
    readParameter(filename,"Space",space,false);
    int filetype;
    readParameter(filename,"Format",filetype,false);
    GrapeIOFileFormatType ftype = static_cast<GrapeIOFileFormatType> (filetype);
    int precision;
    readParameter(filename,"Precision",precision,false);

    if((d != td) || (r != tr) || (n != tn) || (m != tm) )
    {
      std::cerr << d << " | " << td << " DomainField in read!\n";
      std::cerr << r << " | " << tr << " RangeField  in read!\n";
      std::cerr << n << " | " << tn << " in read!\n";
      std::cerr << m << " | " << tm << " in read!\n";
      std::cerr << "Can not initialize DiscreteFunction with wrong FunctionSpace! \n";
      abort();
    }

  const char * path = "";
  GrapeIOStringType fn = genFilename(path,filename,timestep,precision);

  if(ftype == xdr)
    return df.read_xdr(fn);
  if(ftype == ascii)
    return df.read_ascii(fn);
  if(ftype == pgm)
    return df.read_pgm(fn);

  std::cerr << ftype << " GrapeIOFileFormatType not supported at the moment! in file " << __FILE__ << " line " << __LINE__ << "\n"; 
  abort();

  return false;
}


} // end namespace

#endif
