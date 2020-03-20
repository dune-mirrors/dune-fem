#ifndef DUNE_FEM_INPUTOUTPUTTUPLES_HH
#define DUNE_FEM_INPUTOUTPUTTUPLES_HH

#include <sstream>
#include <string>
#include <type_traits>
#include <tuple>
#include <utility>

#include <dune/common/hybridutilities.hh>

#include <dune/grid/common/backuprestore.hh>

#include <dune/fem/io/file/persistencemanager.hh>
#include <dune/fem/io/file/iointerface.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune
{

  namespace Fem
  {

    // IOTupleBase
    // -----------

    struct IOTupleBase
    {
      // take binary stream types from PersistenceManager (overload with define)
      typedef PersistenceManager :: BackupStreamType   OutStreamType ;
      typedef PersistenceManager :: RestoreStreamType  InStreamType ;
      static const bool singleBackupRestoreFile = PersistenceManager :: singleBackupRestoreFile;

      // return path/ prefix name as one string
      static std::string pathAndName(const std::string& path, const std::string& name, const std::string& suffix)
      {
        std::string comboname;
        if (path != "")
          comboname += path;
        else
          comboname += ".";

        comboname += "/";
        comboname += name;
        comboname += suffix;
        return comboname;
      }

      // return grids name
      static std::string gridName(const std::string& path, const std::string& name)
      {
        return pathAndName(path,name,"_grid");
      }

      // return datas name
      static std::string dataName(const std::string& path, const std::string& name)
      {
        return pathAndName(path,name,"_data");
      }

      // return name with rank
      static std::string rankName(const std::string& path, const std::string& name,
                                  const int rank, const int size )
      {
        std::stringstream rankStr;
        const int number = ( singleBackupRestoreFile ) ? size : rank ;
        // add rank if output/input is with different files
        rankStr << "." << number ;
        return pathAndName(path,name,rankStr.str());
      }

      template < class GridType >
      static void
      restoreGrid(GridType *&grid,
                  InStreamType& inStream,
                  const std::string& filename )
      {
        try
        {
          // get standard istream
          std::istream& stream = inStream.stream();
          // restore grid from stream (grid implementation has to take care of byte order)
          grid = BackupRestoreFacility< GridType > :: restore( stream );
        }
        catch ( NotImplemented )
        {
          // write grid to file with filename + extension
          grid = BackupRestoreFacility< GridType > :: restore( filename );
        }

        if( grid == 0 )
          DUNE_THROW(IOError,"Unable to restore grid" );
      }

      template <class GridType>
      static void
      restoreDofManager ( const GridType& grid,
                          InStreamType& inStream )
      {
        // type of DofManager
        typedef DofManager<GridType> DofManagerType;

        // read DofManager's index sets
        DofManagerType& dm =  DofManagerType :: instance ( grid );

        // read data
        dm.read( inStream );

        // resize data
        dm.resizeForRestrict();
      }

      //! write grid and data to given directory
      template < class GridType >
      static
      void writeGrid( const GridType& grid,
                      OutStreamType& outStream,
                      const std::string& filename )
      {
        try
        {
          // get standard ostream
          std::ostream& stream = outStream.stream();
          // write grid to stream (grid implementation has to take care of byte order)
          BackupRestoreFacility< GridType > :: backup( grid, stream );
        }
        catch ( const NotImplemented& )
        {
          // write grid to file with filename + extension
          BackupRestoreFacility< GridType > :: backup( grid, filename );
        }

        // type of DofManager
        typedef DofManager< GridType > DofManagerType;

        // get DofManager reference
        DofManagerType& dm = DofManagerType :: instance ( grid );

        // write DofManager's index sets
        dm.write( outStream );
      }
    };



    // IOTuple
    // -------

    template< class Tuple >
    class IOTuple
    : public IOTupleBase
    {
      template< int N > struct CreateData;
      template< int N > struct RestoreStream;
      template< int N > struct OutputStream;
      template< int N > struct AddToDisplay;
      template< int N > struct AddToDisplayOrRemove;
      template< int N > struct RemoveData;

    public:
      static const int length = std::tuple_size< Tuple >::value;

      typedef Tuple ReturnType ;

      template < class GridType >
      static Tuple *input ( GridType *&grid,
                            double& time,
                            const int rank,
                            const int size,
                            const std::string &path,
                            const std::string &name )
      {
        // true if grid has to be read
        const bool newGrid = (grid == 0);
        assert( newGrid );

        // get filename
        std::string filename = rankName( path, name, rank, size );

        if( Parameter :: verbose () )
          std::cout << "IOTuple: Reading data from " << filename << std::endl;

        // create binary stream
        InStreamType& inStream = *(Fem :: StreamFactory< InStreamType > :: create( filename, rank ));

        if( newGrid )
        {
          // create and read grid
          IOTupleBase::restoreGrid( grid, inStream, filename );
        }

        // create all data
        Tuple *ret = new Tuple;
        Hybrid::forEach( std::make_index_sequence< length >{},
                         [&]( auto i ) { CreateData<i>::apply( *grid, *ret ); });

        if( newGrid )
        {
          // now read dofmanager and index sets
          IOTupleBase::restoreDofManager( *grid, inStream );
        }

        // get simulation time of data
        inStream >> time ;

        // now read all data
        Hybrid::forEach( std::make_index_sequence< length >{},
                         [&]( auto i ) { RestoreStream<i>::apply( inStream, *ret ); });

        if( Parameter :: verbose() )
          std::cout << "    FINISHED!" << std::endl;

        // delete stream which was created by the StreamFactory
        delete &inStream;

        return ret;
      }

      //! restore all data in tupel
      template< class GridType >
      static void restoreData ( Tuple &data,
                                const GridType &grid,
                                const std::string &path,
                                const std::string &name )
      {
        // get filename
        std::string filename =
          rankName( path, name, grid.comm().rank(), grid.comm().size() );

        // create binary stream
        InStreamType inStream( filename );

        // read all data now
        Hybrid::forEach( std::make_index_sequence< length >{},
                         [&]( auto i ) { RestoreStream<i>::apply( inStream, data ); });
      }

      //! write grid and data to given directory
      template< class GridType >
      static void output ( const GridType &grid,
                           const double time,
                           const std::string &path,
                           const std::string &name,
                           const Tuple &tuple )
      {
        // get filename
        std::string filename =
          rankName( path, name, grid.comm().rank(), grid.comm().size() );

        // create binary stream
        OutStreamType& outStream = *(Fem :: StreamFactory< OutStreamType > :: create( filename ));

        // write grid, either to binaryStream or with given filename
        writeGrid( grid, outStream, filename );

        // save simulation time of data
        outStream << time;

        // write data to stream
        Hybrid::forEach( std::make_index_sequence< length >{},
                         [&]( auto i ) { OutputStream<i>::apply( outStream, tuple ); });

        // delete stream created by StreamFactory
        delete &outStream;
      }

      template< class Disp, class DINFO >
      static void addToDisplay ( Disp &disp, const DINFO *dinf, double time, Tuple &tuple )
      {
        Hybrid::forEach( std::make_index_sequence< length >{},
                         [&]( auto i ) { AddToDisplay<i>::apply( disp, dinf, time, tuple ); });
      }

      template< class Disp, class DINFO >
      static void addToDisplayOrRemove ( Disp &disp, const DINFO *dinf, double time, Tuple &tuple )
      {
        Hybrid::forEach( std::make_index_sequence< length >{},
                         [&]( auto i ) { AddToDisplayOrRemove<i>::apply( disp, dinf, time, tuple ); });
      }

      template< class Disp >
      static void addToDisplay ( Disp &disp, Tuple &tuple )
      {
        Hybrid::forEach( std::make_index_sequence< length >{},
                         [&]( auto i ) { AddToDisplay<i>::apply( disp, tuple ); });
      }

      static void removeData ( Tuple &tuple )
      {
        Hybrid::forEach( std::make_index_sequence< length >{},
                         [&]( auto i ) { RemoveData<i>::apply( tuple ); });
      }
    };


    template< class Tuple >
    template< int N >
    struct IOTuple< Tuple >::CreateData
    {
      typedef typename std::remove_pointer< typename std::tuple_element< N, Tuple >::type >::type DiscreteFunction;
      typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
      typedef typename DiscreteFunctionSpace::GridPartType GridPart;

      template< class Grid >
      static void apply ( Grid &grid, Tuple &tuple )
      {
        GridPart *gridPart = new GridPart( grid );
        DiscreteFunctionSpace *space = new DiscreteFunctionSpace( *gridPart );
        std::get< N >( tuple ) = new DiscreteFunction( "", *space );
      }
    };

    template< class Tuple >
    template< int N >
    struct IOTuple< Tuple >::RestoreStream
    {
      typedef typename std::remove_pointer< typename std::tuple_element< N, Tuple >::type >::type DiscreteFunction;

      template< class StreamTraits >
      static void apply ( InStreamInterface< StreamTraits > &inStream, Tuple &tuple )
      {
        DiscreteFunction *df = std::get< N >( tuple );
        bool readDF = false ;
        inStream >> readDF ;
        // if discrete function was written, we can read it
        if( readDF )
        {
          if( df )
            df->read( inStream );
          else
            DUNE_THROW(InvalidStateException,"no discrete function on input");
        }
      }
    };

    template< class Tuple >
    template< int N >
    struct IOTuple< Tuple >::OutputStream
    {
      typedef typename std::remove_pointer< typename std::tuple_element< N, Tuple >::type >::type DiscreteFunction;

      //! apply function writing to stream
      template< class StreamTraits >
      static void apply ( OutStreamInterface< StreamTraits > &outStream, const Tuple &tuple )
      {
        const DiscreteFunction *df = std::get< N >( tuple );
        const bool writeDF = ( df != 0 );

        // write flag whether discrete function was written, for later restore
        outStream << writeDF ;

        // if pointer is valid: write function to stream
        if( writeDF )
        {
          df->write( outStream );
        }
      }
    };

    template< class Tuple >
    template< int N >
    struct IOTuple< Tuple >::AddToDisplay
    {
      // revert tuple order to reverse deletion to creation
      static const int pos = std::tuple_size< Tuple >::value - 1 - N;

      typedef typename std::remove_pointer< typename std::tuple_element< pos, Tuple >::type >::type DiscreteFunction;

      template< class Disp, class DINFO >
      static void apply ( Disp &disp, const DINFO *&dinf, const double &time, Tuple &tuple )
      {
        DiscreteFunction *df = std::get< pos >( tuple );
        if( df )
        {
          assert( dinf->comp );
          std::cout << "adding to display " << dinf->name << std::endl;
          disp.addData( *df, dinf, time );
        }
      }

      template< class Disp >
      static void apply ( Disp &disp, Tuple &tuple )
      {
        DiscreteFunction *df = std::get< pos >( tuple );
        if( df )
        {
          disp.addData( *df );
        }
      }
    };


    template< class Tuple >
    template< int N >
    struct IOTuple< Tuple >::AddToDisplayOrRemove
    {
      template< class Disp, class DINFO >
      static void apply ( Disp &disp, const DINFO *&dinf, const double &time, Tuple &tuple )
      {
        if( dinf->comp )
          AddToDisplay< N >::apply( disp, dinf, time, tuple );
        else
          // RemoveData reverts N itself
          RemoveData< N >::apply( tuple );
      }
    };


    template< class Tuple >
    template< int N >
    struct IOTuple< Tuple >::RemoveData
    {
      // revert tuple order to reverse deletion to creation
      static const int pos = std::tuple_size< Tuple >::value - 1 - N;
      typedef typename std::remove_pointer< typename std::tuple_element< pos, Tuple >::type >::type DiscreteFunction;
      typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
      typedef typename DiscreteFunctionSpace::GridPartType GridPart;

      static void apply ( Tuple &tuple )
      {
        DiscreteFunction *&df = std::get< pos >( tuple );
        if( df )
        {
          const DiscreteFunctionSpace *space = &(df->space());
          const GridPart *gridPart = &(space->gridPart());

          delete df;
          delete space;
          delete gridPart;
          df = 0;
        }
      }
    };



    // IOTuple for empty tuple
    // -----------------------

    template<>
    class IOTuple< std::tuple<> >
    : public IOTupleBase
    {
      typedef std::tuple<> Tuple;

    public:
      static const int length = 0;

      template <class DataIO,class GridType>
      static Tuple *input ( DataIO &dataio, GridType *&grid, double &t, int n,
                            const std::string &path, const std::string &name )
      {
        return new Tuple();
      }

      //! restore all data in tuple
      template< class GridType >
      static void restoreData ( Tuple &data,
                                const GridType &grid,
                                const std::string &path,
                                const std::string &name )
      {
        // nothing to do here
      }

      //! write grid and data to given directory
      template< class GridType >
      static void output ( const GridType &grid,
                           const double time,
                           const std::string &path,
                           const std::string &name,
                           const Tuple &tuple )
      {
        // get filename
        std::string filename =
          rankName( path, name, grid.comm().rank(), grid.comm().size() );

        // create binary stream
        OutStreamType outStream( filename );

        // write grid, either to binaryStream or with given filename
        writeGrid( grid, outStream, filename );

        // save simulation time of data
        outStream << time ;
      }

      template< class Disp, class DINFO >
      static void addToDisplay ( Disp &disp, const DINFO *dinf, double time, Tuple &tuple )
      {}

      template< class Disp, class DINFO >
      static void addToDisplayOrRemove ( Disp &disp, const DINFO *dinf, double time, Tuple &tuple )
      {}

      template< class Disp >
      static void addToDisplay ( Disp &disp, Tuple &tuple )
      {}

      static void removeData ( Tuple &tuple )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_INPUTOUTPUTTUPLES_HH
