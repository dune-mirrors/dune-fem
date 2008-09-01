#ifndef DUNE_GRAPEDATAIO_HH
#define DUNE_GRAPEDATAIO_HH

//- system includes 

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/defaultindexsets.hh>
#include <dune/grid/utility/grapedataioformattypes.hh>

#include <dune/fem/space/common/dofmanager.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh> 

//- Local includes 
#include "asciiparser.hh"

namespace Dune {

inline std::string generateFilename(const std::string& fn,
                               int ntime,
                               int precision = 6)
{
  const char * fakePath = "";
  return genFilename(fakePath,fn,ntime,precision);
}

///////////////////////////////////////////////////
//
//  IndexSet Names
//
///////////////////////////////////////////////////
template <class IndexSetImp>
std::string indexSetToName(const IndexSetImp& set)
{
  return set.name();
}

template <class GridImp>
std::string indexSetToName(const WrappedLevelIndexSet<GridImp>& set)
{
  return "LevelIndexSet";
}

template <class GridImp>
std::string indexSetToName(const WrappedLeafIndexSet<GridImp>& set)
{
  return "LeafIndexSet";
}

template <class GridImp>
std::string indexSetToName(const WrappedHierarchicIndexSet<GridImp>& set)
{
  return "HierarchicIndexSet";
}


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
      const GrapeIOStringType & fnprefix , double & time , int timestep, bool verbose = true );
  
  //! get Grid from file with time and timestep , return grid pointer if ok 
  inline static GridType * restoreGrid (
      const GrapeIOStringType & fnprefix , double & time , int timestep, bool verbose = false )
  {
    // todo MPI_Comm pass to grid type 
    GridType * grid = new GridType (); 
    assert( grid );  
    readGrid(*grid,fnprefix,time,timestep, verbose);
    return grid;
  }
};

template <int dim, int dimworld, class GridImp> 
class GrapeDataIOImp<dim,dimworld,GridImp,false>
{
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
      GrapeIOStringType fnstr = generateFilename(fnprefix,timestep,precision);
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
      const GrapeIOStringType & fnprefix , double & time , int timestep, bool verbose = true )
  {
    std::string gridname;

    bool readGridName = readParameter(fnprefix,"Grid",gridname);
    if(! readGridName ) 
    {
      std::cerr << "P["<< grid.comm().rank() << "] ERROR: Couldn't open file '"<<fnprefix<<"' !" << std::endl;
      return false;
    }

    if( grid.name() != gridname )
    {
      if( (grid.name() == "SGrid") && (gridname == "YaspGrid") ) 
      {
        std::cerr << "WARNING: YaspGrid is read as SGrid! \n";
      }
      else 
      {
        std::cerr << "\nERROR: '" << grid.name() << "' tries to read '" << gridname << "' file. \n";
        abort();
      }
    }

    int precision = 6;
    readParameter(fnprefix,"Precision",precision);

    int hasDm = 0;
    readParameter(fnprefix,"DofManager",hasDm);

    {
      GrapeIOStringType fnstr = generateFilename(fnprefix,timestep,precision);
      
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
      const GrapeIOStringType & fnprefix , double & time , int timestep, bool verbose = false )
  {
    std::string macroName (fnprefix);
#if HAVE_MPI 
    // read global file, only differs for YaspGrid
    macroName += ".macro.global";
#else 
    // read sub grid file 
    macroName += ".macro";
#endif
    
    GridType * grid = 0;
    {
      // create macro grid 
      GridPtr<GridType> gridptr(macroName);
      std::cout << "Created Structured Macro Grid `" << macroName << "' !\n";
      // release pointer 
      grid = gridptr.release();
    }
    assert( grid );
    readGrid(*grid,fnprefix,time,timestep,verbose);
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
      const GrapeIOStringType fnprefix , double & time , int timestep, bool verbose = true )
  {
    const bool hasBackupRestore = Capabilities::hasBackupRestoreFacilities<GridType>::v;
    return GrapeDataIOImp<GridType::dimension,GridType::dimensionworld,GridType,hasBackupRestore>::
      readGrid(grid,fnprefix,time,timestep);
  }

  //! get Grid from file with time and timestep , return true if ok 
  inline GridType * restoreGrid(const GrapeIOStringType fnprefix , double & time , int timestep, bool verbose = false )
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
    std::fstream file (fnprefix.c_str(),std::ios::out);
    if( file.is_open() )
    {
      file << "Grid: "   << grid.name() << std::endl;
      file << "Format: " << ftype <<  std::endl;
      file << "Precision: " << precision << std::endl;
      int writeDm = (hasDm)? 1 : 0;
      file << "DofManager: " << writeDm << std::endl; 

      GrapeIOStringType fnstr = generateFilename(fnprefix,timestep,precision);
      
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
(GridImp & grid, const GrapeIOStringType & fnprefix , double & time , int timestep, bool verbose )
{
  int helpType = (int) xdr;
  std::string gridname;

  bool readGridName = readParameter(fnprefix,"Grid",gridname, verbose);
  if(! readGridName ) 
  {
    if(grid.comm().rank() == 0)
    {
      std::cerr << "P["<< grid.comm().rank() << "] ERROR: Couldn't open file '"<<fnprefix<<"' !" << std::endl;
      abort();
    }
    else if( verbose )
    {
      // on all other procs on print warning
      std::cerr << "P["<< grid.comm().rank() << "] WARNING: Couldn't open file '"<<fnprefix<<"' !" << std::endl;
    }
  }
  else 
  {
    if(grid.name() != gridname)
    {
      std::cerr << "\nERROR: '" << grid.name() << "' tries to read '" << gridname << "' file. \n";
      abort();
    }
  }

  int precision = 6;
  int hasDm = 0;
  readParameter(fnprefix,"Format",    helpType, verbose);
  readParameter(fnprefix,"Precision", precision,verbose);
  readParameter(fnprefix,"DofManager",hasDm,    verbose);

  GrapeIOFileFormatType ftype = (GrapeIOFileFormatType) helpType;

  GrapeIOStringType fn = generateFilename(fnprefix,timestep,precision);
  if( verbose ) 
  {
    std::cout << "Read file: fnprefix = `" << fn << "' \n";
  }

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
  /*
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
  */
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

    enum { n = DiscreteFunctionSpaceType::dimDomain };
    enum { m = DiscreteFunctionSpaceType::dimRange };

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
      file << "IndexSet: " << indexSetToName(df.space().indexSet()) << std::endl;
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

  GrapeIOStringType fn = generateFilename(filename,timestep,precision);

  if(ftype == xdr)
    return df.write_xdr(fn);
  if(ftype == ascii)
    return df.write_ascii(fn);
#if DUNE_FEM_COMPATIBILITY
  if(ftype == pgm)
    return df.write_pgm(fn);
#endif

  DUNE_THROW( NotImplemented, "GrapeIOFileFormatType " << ftype
                               << " currently not supported." );
  return false;
}

template <class GridType>
template <class DiscreteFunctionType> 
inline bool GrapeDataIO<GridType> :: 
readData(DiscreteFunctionType & df, const GrapeIOStringType filename, int timestep)
{
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

    enum { tn = DiscreteFunctionSpaceType::dimDomain };
    enum { tm = DiscreteFunctionSpaceType::dimRange };

    int n,m;
    GrapeIOStringType r,d;
    GrapeIOStringType tr (typeIdentifier<RangeFieldType>());
    GrapeIOStringType td (typeIdentifier<DomainFieldType>());

    std::ifstream check ( filename.c_str() );
    if( !check )
    {
      std::cerr << "WARNING: Couldn't open file `"<< filename << "'! \n";
      return false;
    }

    readParameter(filename,"DomainField",d,false);
    readParameter(filename,"RangeField",r,false);
    readParameter(filename,"Dim_Domain",n,false);
    readParameter(filename,"Dim_Range",m,false);
    int space;
    readParameter(filename,"Space",space,false);
    int order; 
    readParameter(filename,"Polynom_order",order,false);
    if( space != (int) df.space().type() || (order != df.space().order()) )
    {
      derr << "GrapeDataIO::readData: Wrong SpaceType, read (space = " << space << ", order = " << order << ")";
      derr << " but program has (space = "<< df.space().type() << ", order = ";
      derr << df.space().order() << ") "<< std::endl; 
      abort();
    }
    
    std::string indexSetName;
    if( ! readParameter(filename,"IndexSet",indexSetName,false) ) 
    {
      // if parameter not available skip test 
      indexSetName = indexSetToName(df.space().indexSet());
    }
    if( indexSetName != indexSetToName(df.space().indexSet()) )
    {
      derr << "GrapeDataIO::readData: Wrong IndexSet, stored data type is `" << indexSetName;
      derr << "' but type to restore is (change this type) `"<< indexSetToName(df.space().indexSet()) << "'."<< std::endl; 
      abort();
    }
    
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

  GrapeIOStringType fn = generateFilename(filename,timestep,precision);

  if(ftype == xdr)
    return df.read_xdr(fn);
  if(ftype == ascii)
    return df.read_ascii(fn);
#if DUNE_FEM_COMPATIBILITY
  if(ftype == pgm)
    return df.read_pgm(fn);
#endif

  DUNE_THROW( NotImplemented, "GrapeIOFileFormatType " << ftype
                               << " currently not supported." );
  return false;
}


} // end namespace

#endif
