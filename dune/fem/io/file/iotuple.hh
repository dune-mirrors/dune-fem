#ifndef DUNE_INPUTOUTPUTTUPLES_HH
#define DUNE_INPUTOUTPUTTUPLES_HH

//- system includes 
#include <string>

//- Dune includes 
#include <dune/common/forloop.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/tuples.hh>

//- Dune grid includes 
#include <dune/grid/common/backuprestore.hh>

//- Dune fem includes 
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/io/file/persistencemanager.hh>
#include <dune/fem/io/file/iointerface.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/grid/common/backuprestore.hh>

namespace Dune
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
      catch ( NotImplemented ) 
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
    static const int length = tuple_size< Tuple >::value;

  public:
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
      ForLoop< CreateData, 0, length-1 >::apply( *grid, *ret );
      
      if( newGrid ) 
      {
        // now read dofmanager and index sets 
        IOTupleBase::restoreDofManager( *grid, inStream );
      }

      // get simulation time of data 
      inStream >> time ;
      
      // now read all data 
      ForLoop< RestoreStream, 0, length-1 >::apply( inStream, *ret );
      
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
      ForLoop< RestoreStream, 0, length-1 >::apply( inStream, data );
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
      ForLoop< OutputStream, 0, length-1 >::apply( outStream, tuple );

      // delete stream create by StreamFactory
      delete &outStream; 
    }
    
    template< class Disp, class DINFO >
    static void addToDisplay ( Disp &disp, const DINFO *dinf, double time, Tuple &tuple )
    {
      ForLoop< AddToDisplay, 0, length-1 >::apply( disp, dinf, time, tuple );
    }

    template< class Disp, class DINFO >
    static void addToDisplayOrRemove ( Disp &disp, const DINFO *dinf, double time, Tuple &tuple )
    {
      ForLoop< AddToDisplayOrRemove, 0, length-1 >::apply( disp, dinf, time, tuple );
    }

    template< class Disp >
    static void addToDisplay ( Disp &disp, Tuple &tuple )
    {
      ForLoop< AddToDisplay, 0, length-1 >::apply( disp, tuple );
    }
    
    static void removeData ( Tuple &tuple ) 
    {
      ForLoop< RemoveData, 0, length-1 >::apply( tuple );
    }
  };


  template< class Tuple >
  template< int N >
  struct IOTuple< Tuple >::CreateData
  {
    typedef typename TypeTraits< typename tuple_element< N, Tuple >::type >::PointeeType DiscreteFunction;
    typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
    typedef typename DiscreteFunctionSpace::GridPartType GridPart;

    template< class Grid >
    static void apply ( Grid &grid, Tuple &tuple )
    {
      GridPart *gridPart = new GridPart( grid );
      DiscreteFunctionSpace *space = new DiscreteFunctionSpace( *gridPart );
      get< N >( tuple ) = new DiscreteFunction( "", *space );
    }
  };

  template< class Tuple >
  template< int N >
  struct IOTuple< Tuple >::RestoreStream
  {
    typedef typename TypeTraits< typename tuple_element< N, Tuple >::type >::PointeeType DiscreteFunction;

    template< class StreamTraits >
    static void apply ( InStreamInterface< StreamTraits > &inStream, Tuple &tuple ) 
    {
      DiscreteFunction *df = get< N >( tuple );
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
    typedef typename TypeTraits< typename tuple_element< N, Tuple >::type >::PointeeType DiscreteFunction;

    //! apply function writing to stream 
    template< class StreamTraits >
    static void apply ( OutStreamInterface< StreamTraits > &outStream, const Tuple &tuple )
    {
      const DiscreteFunction *df = get< N >( tuple );
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
    typedef typename TypeTraits< typename tuple_element< N, Tuple >::type >::PointeeType DiscreteFunction;

    template< class Disp, class DINFO >
    static void apply ( Disp &disp, const DINFO *&dinf, const double &time, Tuple &tuple )
    {
      DiscreteFunction *df = get< N >( tuple );
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
      DiscreteFunction *df = get< N >( tuple );
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
        RemoveData< N >::apply( tuple );
    }
  };

    
  template< class Tuple >
  template< int N >
  struct IOTuple< Tuple >::RemoveData
  {
    typedef typename TypeTraits< typename tuple_element< N, Tuple >::type >::PointeeType DiscreteFunction;
    typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;
    typedef typename DiscreteFunctionSpace::GridPartType GridPart;

    static void apply ( Tuple &tuple )
    {
      DiscreteFunction *&df = get< N >( tuple );
      const GridPart *gridPart = &(df->space().gridPart());
      const DiscreteFunctionSpace *space = &(df->space());

      delete df;
      delete space;
      delete gridPart;
      df = 0;
    }
  };



  // IOTuple for empty tuple
  // -----------------------

  template<>
  class IOTuple< tuple<> >
  : public IOTupleBase 
  {
    typedef tuple<> Tuple;

  public:
    static const int length = 0;

    template <class DataIO,class GridType>
    static Tuple *input ( DataIO &dataio, GridType *&grid, double &t, int n,
                          const std::string &path, const std::string &name )
    {
      return new Tuple();
    }

    //! restore all data in tupel 
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

} // end namespace Dune 

#endif // #ifndef DUNE_INPUTOUTPUTTUPLES_HH
