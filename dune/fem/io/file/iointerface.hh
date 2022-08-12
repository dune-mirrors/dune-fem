#ifndef DUNE_FEM_IOINTERFACE_HH
#define DUNE_FEM_IOINTERFACE_HH

//- system includes
#include <dirent.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <utility>
#include <sys/stat.h>
#include <sys/types.h>


//- Dune includes
#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>


// defines Parameter
#include <dune/fem/io/parameter.hh>

// binary data io
#include <dune/fem/io/io.hh>

#include <dune/fem/misc/capabilities.hh>

// if grape was configured then include headers
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#endif

namespace Dune
{

  namespace Fem
  {

    // generateFilename
    // ----------------

    inline std::string generateFilename ( const std::string &fn,
                                          int ntime,
                                          int precision = 6 )
    {
      std::ostringstream name;
      name << fn << std::setw( precision ) << std::setfill( '0' ) << ntime;
      return name.str();
    }



    class TimeProviderBase;

    /** @addtogroup DiscFuncIO
       The package dune-fem provides a number
       of possibilities to write data to disk either
       for visualization and debugging purposes or
       for checkpointing.
       Visualization output of discrete functions
       is provided at the moment for
       - GraPE
       - VTK, e.g., for ParaView, VisIt, or others
       - Matlab
       .
       In addition data can be visualized in GraPE online during
       a simulation by setting the parameter
       \b fem.io.grapedisplay to one.
       Except for Matlab output or Dune::Checkpointer (always binary)
       the format is choosen through the parameter \b fem.io.outputformat (see below).

       Utilities for Matlab output are availabe
       through the Dune::MatlabHelper.

       The \ref Dune::Fem::CheckPointer "checkpointing facility"
       writes both the grid state, a number of
       discrete function, and other parameters
       provided by the user. This information can
       then be read from disc to continue the simulation
       from the saved state. Using the \ref Visualization
       "datadisp" programm the checkpoint files can also
       be used for visualization purposes.

       \remark The interface class for general output
               of DiscreteFunctions to disc is given
               by IOInterface. A general purpose data writter
               is provided by the class Dune::DataWriter.

       Since the binary output format is lossless
       it is also used by the Dune::CheckPointer
       class which also writes data to files
       but alternates between to filenames -
       while the Dune::DataWriter should be used
       to store data for postprocessing, the
       Dune::CheckPointer facility should be used
       to be able to restart computations after
       a unexpected termination of a simulation.

       Data files are generated in the directory
       given by the \b fem.prefix parameter and
       with the file prefix chosen via the parameter
       \b fem.io.datafileprefix. The data to be
       saved to disk is given to the Dune::DataWriter instance
       through a reference to a std::tuple of
       discrete function pointer.

       For a time series, data can be either written for a fixed
       time step, or after a fixed number of iterations using the
       parameters
       \b fem.io.savestep or \b fem.io.savecount,
       respectivly.
       If a series of data is to be written without a real
       time variable available, e.g., a series of refined grids,
       startTime=0, endTime=100,
       \b fem.io.savestep=-1 and \b fem.io.savecount=1
       is a good choice to make; data is then written
       using \b datawriter.write(step,step).

       The following code snippet demonstrated the
       general usage of the Dune::DataWriter:
       \code
       typedef std::tuple< DestinationType > IOTupleType;
       IOTupleType dataTup ( &U );
       typedef DataWriter< GridType, IOTupleType > DataWriterType;
       DataWriterType dataWriter( grid,gridfilename,dataTup,startTime,endTime );
       for (counter=0; time<endTime; ++counter) {
         ...
         dataWriter.write(time,counter);
       }
       \endcode

        \femparam{fem.prefix, path used for all file output, .}
        \femparam{fem.io.datafileprefix, prefix used for all data files}
        \femparam{fem.io.outputformat, output format, vtk-vertex}
                  values are:
                - vtk-cell =  VTK cell data
                - vtk-vertex = VTK vertex data (might require interpolation into Lagrange space)
                - sub-vtk-cell =  VTK cell data including a specified level of subsampling
                - binary = write binary data (only available with DataWriter and CheckPointer)
                - gnuplot = write gunplot compatible data
                - none = no data output
                .
        \femparam{fem.io.grapedisplay, use grape for online visualization; default is false}
        \femparam{fem.io.savestep, interval for writting data files}
         use value <0 to deativate
        \femparam{fem.io.savecount, number of time steps between writting
                                    file}
        use value <0 to deactivate
    **/


    /** @ingroup DiscFuncIO
     \brief IOInterface to write data to hard disk
     \interfaceclass
    */
    class IOInterface {

    protected:
      //! default constructor
      IOInterface() {}

    public:
      //! destructor
      virtual ~IOInterface () {}

      /** \brief write data with a given sequence stamp
          \param sequenceStamp stamp for the data set
      */
      virtual void writeData ( double sequenceStamp ) const = 0;

      /** \brief write given data to disc, evaluates parameter savecount and savestep
          \param tp  time provider for time and step
      */
      virtual void write( const TimeProviderBase& tp ) const = 0;

      /** \brief write given data to disc, evaluates parameter savecount
      */
      virtual void write() const = 0;

      //! return FEM key for macro grid reading
      static std::string defaultGridKey ( int dimension, bool check = true )
      {
        return defaultGridKey( dimension, Parameter::container(), check );
      }

      static std::string defaultGridKey ( int dimension, const ParameterReader &parameter, bool check = true )
      {
        return defaultGridKey( "fem.io.macroGridFile", dimension, parameter, check );
      }

      static std::string defaultGridKey ( std::string base, int dimension, bool check = true )
      {
        return defaultGridKey( std::move( base ), dimension, Parameter::container(), check );
      }

      //! return FEM key for macro grid reading
      static std::string defaultGridKey ( std::string base, int dimension, const ParameterReader &parameter, bool check = true )
      {
        const std::string oldGridKey( base );

        std::ostringstream gridKeyStream;
        gridKeyStream << oldGridKey << "_" << dimension << "d";
        const std::string newGridKey( gridKeyStream.str() );

        // check for old parameter
        if( parameter.exists( oldGridKey ) )
        {
          if( parameter.exists( newGridKey ) )
          {
            std::cerr << "WARNING: ignoring `" << oldGridKey << "' because `"
                      << newGridKey << "' was also found in parameter file." << std::endl;
            return newGridKey;
          }
          else
          {
            std::cerr << "WARNING: change `" << oldGridKey << "' to `"  << newGridKey
                      << "' in parameter file." << std::endl;
            return oldGridKey;
          }
        }

        // check for parameter with dimension
        if( check && !parameter.exists( newGridKey ) )
        {
          std::cerr << "ERROR: Parameter `" << newGridKey << "' not found." << std::endl;
          DUNE_THROW( ParameterNotFound, "Parameter `" << newGridKey << "' not found." );
        }
        return newGridKey;
      }

      //! create given path in combination with rank
      static void createPath ( const std::string &path )
      {
        if( !createDirectory( path ) )
          std::cerr << "Failed to create path `" << path << "'." << std::endl;
      }

      //! create given path in combination with rank
      static std::string createPathName(const std::string& pathPref, int rank )
      {
        std::string path(pathPref);

        // add proc number to path
        {
          path += "_";
          std::stringstream rankDummy;
          rankDummy << rank;
          path += rankDummy.str();
        }
        return path;
      }

      //! standard path reading and creation method
      //! rank is added to output path
      static std::string readPath()
      {
        return Parameter::commonOutputPath();
      }

      /** \brief create global path for data output */
      template <class CommunicatorType>
      static void createGlobalPath(const CommunicatorType& comm,
              const std::string& path)
      {
        // only rank 0 creates global dir
        if( comm.rank() <= 0 )
        {
          // create directory
          if( !createDirectory( path ) )
            std::cerr << "Failed to create path `" << path << "'." << std::endl;
        }

        // wait for all procs to arrive here
        comm.barrier ();
      }

      // copy path to filename and add a slash if necessary
      static std::string copyPathToFilename( const std::string& path )
      {
        // first proc creates directory
        std::string filename( path );

        const char lastToken = filename.c_str()[ filename.size() - 1 ];
        const char* slash = "/";
        // add / if necessary
        if( lastToken != slash[0] )
          filename += "/";

        return filename;
      }

      // creates path and processor sub pathes
      template <class CommunicatorType>
      static std::string createPath(const CommunicatorType& comm,
              const std::string& pathPrefix,
              const std::string& dataPrefix,
              const int step,
              const bool alsoCreateRankPath = true )
      {
        // first proc creates directory
        std::string filename( copyPathToFilename( pathPrefix ));

        filename += dataPrefix;
        std::string path = generateFilename( filename, step );

        // create global path
        createGlobalPath( comm, path );

        // also create path for each rank
        if( alsoCreateRankPath )
        {
          // append path with p for proc
          path += "/p";

          // create path if not exists
          path = createPathName( path, comm.rank() );

          // create path if not exits
          if( !createDirectory( path ) )
            std::cerr << "Failed to create path `" << path << "'." << std::endl;
        }
        return path;
      }

      // creates path and processor sub pathes
      static std::string createRecoverPath(
              const std::string& pathPrefix,
              const int rank,
              const std::string& dataPrefix,
              const int step,
              const bool alsoUseRankPath = true )
      {
        // first proc creates directory
        std::string filename( copyPathToFilename( pathPrefix ));

        filename += dataPrefix;
        std::string path = generateFilename( filename, step );

        if( alsoUseRankPath )
        {
          // append path with p for proc
          path += "/p";

          // create proc dir
          return createPathName( path , rank );
        }
        else
          return path;
      }

      //! if grid is structured grid, write macro file
      template <class GridImp>
      static void writeMacroGrid(const GridImp& grid,
                                 std::ostream& out,
                                 const std::string& macroname,
                                 const std::string& path,
                                 const std::string& prefix,
                                 const bool writeSubFiles = false )
      {
        // do nothing if SaveParallelCartesianGrid is not specified
        if( ! SaveParallelCartesianGrid< GridImp > :: saveMacroGrid ) return ;

        // create file descriptor
        std::ifstream gridin(macroname.c_str());
        if( !gridin)
        {
          std::cerr << "Couldn't open file `" << macroname << "' ! \n";
          return ;
        }

        // read interval information of structured grid
        dgf::IntervalBlock interval(gridin);
        if(!interval.isactive())
        {
          std::cerr<<"Did not find IntervalBlock in macro grid file `" << macroname << "' ! \n";
          return;
        }

        std::string filename(path);
        filename += "/";
        filename += prefix;
        filename += "_grid";

        saveCartesianGrid( grid, out, interval, filename, writeSubFiles );
        return;
      }

      //! if grid is structured grid, write macro file
      template <class GridImp>
      static void copyMacroGrid(const GridImp& g,
                                const std::string& macroGrid,
                                const std::string& orgPath,
                                const std::string& destPath,
                                const std::string& prefix)
      {
        // do nothing if SaveParallelCartesianGrid is not specified
        if( ! SaveParallelCartesianGrid< GridImp > :: saveMacroGrid ) return ;

        if( macroGrid != "" )
        {
          std::string destFilename(destPath);
          destFilename += "/";
          destFilename += prefix;
          destFilename += "_grid.macro";

          std::ofstream file( destFilename.c_str() );
          if( file.is_open() )
            file << macroGrid;
          else
          {
            if( Parameter :: verbose () )
              std::cerr << "Unable to open: '" << destFilename << "'." << std::endl;
          }
        }
      }

    protected:
      //! create string containing rank
      static std::string strRank(const int rank)
      {
        std::stringstream tmp;
        tmp << "." << rank;
        return tmp.str();
      }

      template <class Grid>
      struct SaveParallelCartesianGrid
      {
        static const bool saveMacroGrid = false ;
        typedef FieldVector<int, Grid::dimension> iTupel;

        static void getCoordinates(const Grid&, const iTupel& , iTupel&, iTupel& ,iTupel& )
        {
          DUNE_THROW(NotImplemented,"SaveParallelCartesianGrid not implemented for choosen GridType");
        }
      };

      template < int dim, class CoordCont >
      struct SaveParallelCartesianGrid< YaspGrid< dim, CoordCont > >
      {
        static const bool saveMacroGrid = true  ;
        typedef YaspGrid< dim, CoordCont >  Grid;
        typedef FieldVector<int, dim> iTupel;

        static void getCoordinates(const Grid& grid, const iTupel& anz,
                                   iTupel& origin, iTupel& originInterior,
                                   iTupel& lengthInterior )
        {
          // Yasp only can do origin = 0
          origin = 0;

          enum { tag = Grid::tag };
          YLoadBalance< dim > loadBalancer;
          Torus< typename Grid::Communication, dim > torus( grid.comm(), tag, anz, &loadBalancer );
          torus.partition( torus.rank(), origin, anz, originInterior, lengthInterior );
        }
      };

#if HAVE_DUNE_SPGRID
      template< class ct, int dim, template< int > class Strategy, class Comm >
      struct SaveParallelCartesianGrid< SPGrid< ct, dim, Strategy, Comm > >
      {
        static const bool saveMacroGrid = true;

        typedef SPGrid< ct, dim, Strategy, Comm > Grid;
        typedef FieldVector< int, dim > iTupel;

        static void getCoordinates( const Grid &grid, const iTupel &anz,
                                    iTupel &origin, iTupel &originInterior,
                                    iTupel &lengthInterior )
        {
#if HAVE_MPI
          typedef typename Grid::MultiIndex MultiIndex;
          const MultiIndex begin = grid.gridLevel( 0 ).localCube().begin();
          const MultiIndex end = grid.gridLevel( 0 ).localCube().end();
          for( int i = 0; i < dim; ++i )
          {
            originInterior[ i ] = begin[ i ];
            lengthInterior[ i ] = end[ i ] - begin[ i ];
          }
#endif // #if HAVE_MPI
        }
      };
#endif // #if HAVE_DUNE_SPGRID

      //! write my partition as macro grid
      template <class GridImp>
      static void saveCartesianGrid (const GridImp& grid,
                                     std::ostream& out,
                                     dgf::IntervalBlock& intervalBlock,
                                     std::string filename,
                                     const bool writeSubFiles )
      {
        enum { dimension = GridImp :: dimension };
        const int rank = grid.comm().rank();

        FieldVector<double,dimension> lang;
        FieldVector<int,dimension>    anz;
        FieldVector<double,dimension>  h;
        FieldVector<int,dimension>    orig;

        if( intervalBlock.numIntervals() != 1 )
        {
          std::cerr << "Warning: Only 1 interval block is handled by "
                    << "IOInterface::saveMacroGridImp" << std::endl;
        }

        typedef typename dgf::IntervalBlock::Interval Interval;
        const Interval &interval = intervalBlock.get( 0 );
        for( int i = 0; i < dimension; ++i )
        {
          orig[ i ] = interval.p[ 0 ][ i ];
          lang[ i ] = interval.p[ 1 ][ i ] - interval.p[ 0 ][ i ];
          anz[ i ] = interval.n[ i ];
          h[ i ] = lang[ i ] / anz[ i ];
        }

        // write sub grid for this rank
        {
          std::string subfilename;
          if( writeSubFiles )
          {
            subfilename = filename;
            // add rank
            subfilename += strRank(rank);
          }

#if HAVE_MPI
          {
            typedef FieldVector<int,dimension> iTupel;

            // origin is zero
            iTupel o( orig );

            iTupel o_interior;
            iTupel s_interior;

            SaveParallelCartesianGrid< GridImp > ::
              getCoordinates( grid, anz, o, o_interior, s_interior );

            FieldVector<double,dimension> origin;
            for(int i=0; i<dimension; ++i)
              origin[ i ] = o[ i ];

            FieldVector<double,dimension> sublang(0.0);
            for(int i=0; i<dimension; ++i)
            {
              origin[i] = o_interior[i] * h[i];
              sublang[i] = origin[i] + (s_interior[i] * h[i]);
            }

            writeStructuredGrid(subfilename,out,origin,sublang,s_interior);
          }
#else
          {
            // in serial this should be zero
            assert( rank == 0 );
            FieldVector<double,dimension> zero(0.0);
            writeStructuredGrid(subfilename,out,zero,lang,anz);
          }
#endif
        }

        // write global grid on rank 0
        if (rank == 0 )
        {
          // write global file for recovery
          filename += ".global";
          FieldVector<double,dimension> zero(0.0);
          std::stringstream dummy;
          writeStructuredGrid(filename,dummy,zero,lang,anz);
        }
      }

      template <int dimension>
      static void writeToStream(std::ostream& file,
                                const FieldVector<double,dimension>& origin,
                                const FieldVector<double,dimension>& lang,
                                const FieldVector<int,dimension>& anz)
      {
        file << "DGF" << std::endl;
        file << "Interval" << std::endl;
        // write first point
        for(int i=0;i<dimension; ++i)
        {
          file << origin[i] << " ";
        }
        file << std::endl;
        // write second point
        for(int i=0;i<dimension; ++i)
        {
          file << lang[i] << " ";
        }
        file << std::endl;
        // write number of intervals in each direction
        for(int i=0;i<dimension; ++i)
        {
          file << anz[i] << " ";
        }
        file << std::endl;
        file << "#" << std::endl;

        file << "BoundaryDomain" << std::endl;
        file << "default 1" << std::endl;
        file << "#" << std::endl;
      }

      //! write structured grid as DGF file
      template <int dimension>
      static void writeStructuredGrid(const std::string& filename,
                                      std::ostream& out,
                                      const FieldVector<double,dimension>& origin,
                                      const FieldVector<double,dimension>& lang,
                                      const FieldVector<int,dimension>& anz)
      {
        writeToStream( out, origin, lang, anz);

        if( filename != "" )
        {
          std::ofstream file (filename.c_str());
          if( file.is_open())
          {
            writeToStream( file, origin, lang, anz);
          }
          else
          {
            std::cerr << "Couldn't open file `" << filename << "' !\n";
          }
        }
      }
    }; // end class IOInterface

  } // end namespace Fem

} // end namespace Dune
#endif // #ifndef DUNE_FEM_IOINTERFACE_HH
