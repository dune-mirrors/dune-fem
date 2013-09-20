#ifndef DUNE_FEM_BINARYDATAIO_HH
#define DUNE_FEM_BINARYDATAIO_HH

#warning "BinaryDataIO is deprecated and will be removed after dune-fem 1.4."

//- Dune includes 
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/utility/grapedataioformattypes.hh>

#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/gridpart/defaultindexsets.hh>
#include <dune/fem/misc/gridname.hh> 

//- Local includes 
#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{
  namespace Fem 
  {

    inline std::string generateFilename(const std::string& fn,
                                        int ntime,
                                        int precision = 6)
    {
      std::ostringstream name;
      name << fn << std::setw(precision) << std::setfill('0') << ntime;

      // Return the string corresponding to the stringstream
      return name.str();
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


    template< int dim, int dimworld, class GridImp, bool hasBackupRestore >
    struct BinaryDataIOImp;


    template< int dim, int dimworld, class GridImp >
    struct BinaryDataIOImp< dim, dimworld, GridImp, true >
    {
      typedef GridImp GridType;


       /** Write Grid with GridType file filename and time 
       *
       * This method uses the Grid Interface Method writeGrid 
       * to actually write the grid, within this method the real file name is
       * generated out of filename and timestep 
       */
      inline static bool writeGrid (const GridType & grid, 
        const GrapeIOFileFormatType ftype, const std::string & fnprefix 
          , double time=0.0, int timestep=0, int precision = 6);

      //! get Grid from file with time and timestep , return true if ok 
      inline static bool readGrid (GridType & grid, 
          const std::string & fnprefix , double & time , int timestep, bool verbose = true );
      
      //! get Grid from file with time and timestep , return grid pointer if ok 
      inline static GridType * restoreGrid (
          const std::string & fnprefix , double & time , int timestep, bool verbose = false )
      {
        // todo MPI_Comm pass to grid type 
        GridType * grid = new GridType (); 
        assert( grid );  
        readGrid(*grid,fnprefix,time,timestep, verbose);
        return grid;
      }
    };

    template< int dim, int dimworld, class GridImp >
    struct BinaryDataIOImp< dim, dimworld, GridImp, false >
    {
      typedef GridImp GridType;

       /** Write structurd grid to file filename 
        NOTE: the macro grid file of this structured grid 
              has to stored in a file named "filename.macro" and
              the format must be Dune DGF. 
       */
      inline static bool writeGrid (const GridType & grid, 
        const GrapeIOFileFormatType ftype, const std::string & fnprefix 
          , double time=0.0, int timestep=0, int precision = 6);

      //! get Grid from file with time and timestep , return true if ok 
      inline static bool readGrid (GridType & grid, 
          const std::string & fnprefix , double & time , int timestep, bool verbose = true );

      //! get Grid from file with time and timestep , return true if ok 
      inline static GridType * restoreGrid (
          const std::string & fnprefix , double & time , int timestep, bool verbose = false );
    };

    template< class GridImp >
    struct BinaryDataIO 
    {
      typedef GridImp GridType;

    private:
      static const char endAsciiHeader = 0x00 ;
      static const int dimGrid = GridType::dimension;
      static const int dimWorld = GridType::dimensionworld;

      static const bool useBackupRestore
        = Dune::Capabilities::hasBackupRestoreFacilities< GridType >::v ;

      typedef BinaryDataIOImp< dimGrid, dimWorld, GridType, useBackupRestore > Impl;

    public:
      BinaryDataIO () {}

      /** Write Grid with GridType file filename and time 
       *
       *  This method uses the Grid Interface Method writeGrid 
       *  to actually write the grid, within this method the real file name is
       *  generated out of filename and timestep 
       */
      inline bool writeGrid (const GridType & grid, 
        const GrapeIOFileFormatType ftype, const std::string fnprefix 
          , double time=0.0, int timestep=0, int precision = 6) const
      {
        return Impl::writeGrid( grid, ftype, fnprefix, time, timestep, precision );
      }

      //! get Grid from file with time and timestep , return true if ok 
      inline bool readGrid (GridType & grid, 
          const std::string fnprefix , double & time , int timestep, bool verbose = true )
      {
        return Impl::readGrid( grid, fnprefix, time, timestep );
      }

      //! get Grid from file with time and timestep , return true if ok 
      inline GridType * restoreGrid(const std::string fnprefix , double & time , int timestep, bool verbose = false )
      {
        return Impl::restoreGrid( fnprefix, time, timestep );
      }

      /**
        Write DiscreteFunctions  
      */ 
      //! write disc func information file and write dofs to file+timestep
      //! this method use the write method of the implementation of the
      //! discrete function
      template <class DiscreteFunctionType>
      inline bool writeData(const DiscreteFunctionType & df,
         const GrapeIOFileFormatType ftype, const std::string filename, 
          int timestep, int precision = 6);

      /**
        Write data file header   
      */ 
      //! write disc func information file and write dofs to file+timestep
      //! this method use the write method of the implementation of the
      //! discrete function
      template <class StreamType, class DiscreteFunctionType>
      inline void writeDataHeader(StreamType& out, 
                                  const DiscreteFunctionType & df,
                                  const GrapeIOFileFormatType ftype,
                                  const int multipleFiles,
                                  const int precision = 6 ) const;

      //! check header of data file 
      template <class DiscreteFunctionType>
      inline bool checkDataHeader(const DiscreteFunctionType & df,
                                  const std::string& filename, 
                                  GrapeIOFileFormatType& ftype,
                                  int& precision, int& multipleFiles ) const;

      //! same as write only read
      template <class DiscreteFunctionType>
      inline bool readData(DiscreteFunctionType & df,
                           const std::string filename, int timestep);
    };


    /** \cond */
    template< int dim, int dimworld, class GridImp >
    inline bool BinaryDataIOImp< dim, dimworld, GridImp, true >::writeGrid
    (const GridImp & grid,
      const GrapeIOFileFormatType ftype, const std::string & fnprefix , 
      double time, int timestep, int precision )
    {
      bool hasDm = false;
      // write dof manager, that corresponds to grid 
      {
        typedef DofManager<GridImp> DofManagerType; 

        std::string dmname(fnprefix);
        dmname += "_dm";
        hasDm = DofManagerType :: write(grid,dmname,timestep);
      }
     
      // write Grid info only on rank 0 to avoid to many files 
      if( grid.comm().rank() == 0 )
      {
        std::fstream file (fnprefix.c_str(),std::ios::out);
        if( file.is_open() )
        {
          file << "Grid: "   << Fem::gridName(grid) << std::endl;
          file << "Format: " << ftype <<  std::endl;
          file << "Precision: " << precision << std::endl;
          int writeDm = (hasDm)? 1 : 0;
          file << "DofManager: " << writeDm << std::endl; 
        }
        else 
        {
          std::cerr << "Couldn't open file `" << fnprefix << "' ! \n";
          return false;
        }

        file.close();
      }

      std::string fnstr = generateFilename(fnprefix,timestep,precision);
      
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

    template< int dim, int dimworld, class GridImp >
    inline bool BinaryDataIOImp< dim, dimworld, GridImp, true >::readGrid 
    (GridImp & grid, const std::string & fnprefix , double & time , int timestep, bool verbose )
    {
      int helpType = (int) xdr;
      int precision = 6;
      int hasDm = 1;

      if( grid.comm().rank() == 0 )
      {
        std::string gridname;
        bool readGridName = readParameter(fnprefix,"Grid",gridname, verbose);
        if( ! readGridName ) 
        {
          if( verbose )
          {
            // on all other procs on print warning
            std::cerr << "P["<< grid.comm().rank() << "] WARNING: Couldn't open file '"<<fnprefix<<"' !" << std::endl;
          }
        }
        else 
        {
          std::string grName ( Fem::gridName( grid ) );
          if( grName != gridname)
          {
            std::cerr << "\nWARNING: '" << grName << "' tries to read '" << gridname << "' file. \n";
          }
        }

        if( readGridName ) 
        {
          readParameter(fnprefix,"Format",    helpType, verbose);
          readParameter(fnprefix,"Precision", precision,verbose);
          readParameter(fnprefix,"DofManager",hasDm,    verbose);
        }
      }

      GrapeIOFileFormatType ftype = (GrapeIOFileFormatType) helpType;
      std::string fn = generateFilename(fnprefix,timestep,precision);
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
     
      return succeded;
    }


    template< int dim, int dimworld, class GridImp >
    inline bool BinaryDataIOImp< dim, dimworld, GridImp, false >
      ::writeGrid ( const GridType &grid,
                    const GrapeIOFileFormatType ftype, const std::string &fnprefix,
                    double time, int timestep, int precision )
    {
      // write dof manager, that corresponds to grid 
      bool hasDm = false;
      {
        typedef DofManager<GridImp> DofManagerType; 

        std::string dmname(fnprefix);
        dmname += "_dm";
        hasDm = DofManagerType :: write(grid,dmname,timestep);
      }
     
      // only write grid info on rank 0 
      if( grid.comm().rank() == 0 )
      {
        std::ofstream file (fnprefix.c_str());
        if( file.is_open() )
        {
          file << "Grid: "   << Fem::gridName(grid) << std::endl;
          file << "Format: " << ftype <<  std::endl;
          file << "Precision: " << precision << std::endl;
          int writeDm = (hasDm)? 1 : 0;
          file << "DofManager: " << writeDm << std::endl; 
          // also write time and maxLevel needed for restore 
          file << "MaxLevel: " << grid.maxLevel() << std::endl;
          file << "Time: " << std::scientific << time << std::endl;
          file.close();
        }
        else 
        {
          std::cerr << "Couldn't open file `" << fnprefix << "' ! \n";
          return false;
        }
      }
      return true;
    }


    template< int dim, int dimworld, class GridImp >
    inline bool BinaryDataIOImp< dim, dimworld, GridImp, false >
      ::readGrid ( GridType &grid,
                   const std::string &fnprefix,
                   double &time, int timestep, bool verbose )
    {
      std::string gridname;

      const int rank = grid.comm().rank() ;

      // read stored maxLevel 
      int maxLevel = 0;

      if( rank == 0 ) 
      {
        bool readGridName = readParameter(fnprefix,"Grid",gridname);
        if( ! readGridName )  
        {
          std::cerr << "P["<< rank << "] ERROR: Couldn't open file '"<<fnprefix<<"' !" << std::endl;
          abort();
        }

        std::string grName ( Fem::gridName( grid ) );
        if( grName != gridname )
        {
          if( (grName == "SGrid") && (gridname == "YaspGrid") ) 
          {
            std::cerr << "WARNING: YaspGrid is read as SGrid! \n";
          }
          else 
          {
            std::cerr << "\nWARNING: '" << grName << "' tries to read '" << gridname << "' file. \n";
          }
        }

        int precision = 6;
        readParameter(fnprefix,"Precision",precision);

        int hasDm = 0;
        readParameter(fnprefix,"DofManager",hasDm);

        const bool foundMaxLevel = readParameter(fnprefix,"MaxLevel",maxLevel);

        // also read time 
        readParameter(fnprefix,"Time",time);

        // if we have old format 
        if( ! foundMaxLevel ) 
        {
          std::string fnstr = generateFilename(fnprefix,timestep,precision);
          readParameter(fnstr,"Time",time);
          readParameter(fnstr,"MaxLevel",maxLevel);
        }

        // calculate level to achieve 
        maxLevel -= grid.maxLevel();
      }

      // get max level for each core 
      maxLevel = grid.comm().max( maxLevel );

      if( maxLevel < 0 ) 
      {
        DUNE_THROW(InvalidStateException,"maxLevel of grid is already to big!");
      }

      // refine grid 
      grid.globalRefine( maxLevel );
      return true;
    }


    template< int dim, int dimworld, class GridImp >
    inline typename BinaryDataIOImp< dim, dimworld, GridImp, false >::GridType *
    BinaryDataIOImp< dim, dimworld, GridImp, false >
      ::restoreGrid ( const std::string &fnprefix, double &time, int timestep, bool verbose )
    {
      std::string macroName (fnprefix);

      if( MPIManager :: size() > 1 )
      {
        // read global file, only differs for YaspGrid
        macroName += ".macro.global";
      }
      else 
      {
        // read sub grid file 
        macroName += ".macro";
      }
      
      GridType * grid = 0;
      {
        if( Parameter :: verbose () )
          std::cout << "Read grid from " << macroName << std::endl;

        {
          std::ifstream testfile ( macroName.c_str() );
          if( ! testfile.is_open() ) 
          {
            std::cerr << "Macro grid file `" << macroName << "' not found!" << std::endl;
          }
        }

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
    /** \endcond */

    template <class GridType>
    template <class StreamType, class DiscreteFunctionType> 
    inline void BinaryDataIO<GridType> :: 
    writeDataHeader( StreamType& file, 
                     const DiscreteFunctionType & df, 
                     const GrapeIOFileFormatType ftype,
                     const int multipleFiles,
                     const int precision ) const
    {
      typedef typename DiscreteFunctionType::FunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

      enum { n = DiscreteFunctionSpaceType::dimDomain };
      enum { m = DiscreteFunctionSpaceType::dimRange };

      std::string d = typeIdentifier<DomainFieldType>();
      std::string r = typeIdentifier<RangeFieldType>();

      file << "DUNE-FEM-Version-Id: " << DUNE_MODULE_VERSION_ID(DUNE_FEM) << std::endl;
      file << "DomainField: " << d << std::endl;
      file << "RangeField: " << r << std::endl;
      file << "Dim_Domain: " << n << std::endl;
      file << "Dim_Range: " << m << std::endl;
      file << "Space: " << df.space().type() << std::endl;
      file << "IndexSet: " << Fem::indexSetToName(df.space().indexSet()) << std::endl;
      file << "Format: " << ftype << std::endl;
      file << "Precision: " << precision << std::endl;
      file << "Polynom_order: " << df.space().order() << std::endl;
      file << "DataBase: " << df.name() << std::endl;
      file << "MultipleFiles: " << multipleFiles << std::endl;
      // write character code for end of ascii header 
      file << endAsciiHeader ; 
    }

    template <class GridType>
    template <class DiscreteFunctionType> 
    inline bool BinaryDataIO<GridType> :: writeData(const DiscreteFunctionType & df, 
    const GrapeIOFileFormatType ftype, const std::string filename, int timestep, int  precision )
    {
      const int multipleFiles = Parameter :: getValue< int > ("fem.io.multiplefiles", 0 );

      {
        std::ofstream file( filename.c_str() );
        if( file.is_open() )
        {
          writeDataHeader( file, df, ftype, multipleFiles, precision );
        }
        else 
        {
          std::cerr << "Couldn't open file `" << filename << "' ! \n";
          return false;
        }
      }

      // for single file use XDR streams
      if( multipleFiles == 0 ) 
      {
        // create xdr stream 
        const bool append = true ;
        XDRFileOutStream xdrOut( filename, append );

        // write discrete function 
        df.write( xdrOut );
      }
      else 
      {
        std::string fn = generateFilename(filename,timestep,precision);

        if(ftype == xdr)
          return df.write_xdr(fn);
        if(ftype == ascii)
          return df.write_ascii(fn);

        DUNE_THROW( NotImplemented, "GrapeIOFileFormatType " << ftype
                                     << " currently not supported." );
      }

      return false;
    }

    template <class GridType>
    template <class DiscreteFunctionType> 
    inline bool BinaryDataIO<GridType> :: 
    checkDataHeader(const DiscreteFunctionType & df, 
                    const std::string& filename, 
                    GrapeIOFileFormatType& ftype,
                    int& precision, int& multipleFiles ) const
    {
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

      enum { tn = DiscreteFunctionSpaceType::dimDomain };
      enum { tm = DiscreteFunctionSpaceType::dimRange };

      int n,m;
      std::string r,d;
      std::string tr (typeIdentifier<RangeFieldType>());
      std::string td (typeIdentifier<DomainFieldType>());

      unsigned int versionId = 0;
      readParameter(filename,"DUNE-FEM-Version-Id",versionId,false);

      // check version 
      unsigned int duneFemId = DUNE_MODULE_VERSION_ID(DUNE_FEM);
      if( versionId != 0 && duneFemId < versionId ) 
      {
        std::cerr << "WARNING: DUNE-FEM-Version-Id from data set is newer than the program ones!" << std::endl;
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
        derr << "BinaryDataIO::readData: Wrong SpaceType, read (space = " << space << ", order = " << order << ")";
        derr << " but program has (space = "<< df.space().type() << ", order = ";
        derr << df.space().order() << ") "<< std::endl; 
        abort();
      }
      
      std::string indexSetName;
      if( ! readParameter(filename,"IndexSet",indexSetName,false) ) 
      {
        // if parameter not available skip test 
        indexSetName = Fem::indexSetToName(df.space().indexSet());
      }
      if( indexSetName != Fem::indexSetToName(df.space().indexSet()) )
      {
        derr << "BinaryDataIO::readData: Wrong IndexSet, stored data type is `" << indexSetName;
        derr << "' but type to restore is (change this type) `"<< Fem::indexSetToName(df.space().indexSet()) << "'."<< std::endl; 
        abort();
      }
      
      int filetype;
      readParameter(filename,"Format",filetype,false);
      ftype = static_cast<GrapeIOFileFormatType> (filetype);
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

      readParameter(filename,"MultipleFiles",multipleFiles,false);

      return true ;
    }

    template <class GridType>
    template <class DiscreteFunctionType> 
    inline bool BinaryDataIO<GridType> :: 
    readData(DiscreteFunctionType & df, const std::string filename, int timestep)
    {
      std::ifstream check ( filename.c_str() );
      if( !check )
      {
        std::cerr << "WARNING: Couldn't open file `"<< filename << "'! \n";
        return false;
      }

      GrapeIOFileFormatType ftype = xdr ;
      int precision = 6 ;
      int multipleFiles = 0 ;

      checkDataHeader( df, filename, ftype, precision, multipleFiles );

      if( multipleFiles == 0 )
      {
        std::ifstream posSeek ( filename.c_str() );

        // search for starting position 
        char c = posSeek.get();
        while( posSeek && c != endAsciiHeader ) 
        {
          c = posSeek.get();
        }

        if( ! posSeek ) 
        {
          std::cerr << "BinaryDataIO::readData: reached unexpected end of file " << filename << std::endl; 
          abort();
        }

        // get current position 
        const size_t pos = posSeek.tellg();
        // close file
        posSeek.close();

        // create xdr streams 
        XDRFileInStream xdrIn( filename, pos );

        // read discrete function
        df.read( xdrIn );
      }
      else 
      {
        // old way of file reading 
        std::string fn = generateFilename(filename,timestep,precision);

        if(ftype == xdr)
          return df.read_xdr(fn);
        if(ftype == ascii)
          return df.read_ascii(fn);

        DUNE_THROW( NotImplemented, "GrapeIOFileFormatType " << ftype
                                     << " currently not supported." );
      }
      return false;
    }


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BINARYDATAIO_HH
