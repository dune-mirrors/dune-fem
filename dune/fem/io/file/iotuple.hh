#ifndef DUNE_INPUTOUTPUTTUPLES_HH
#define DUNE_INPUTOUTPUTTUPLES_HH

//- system includes 
#include <string>

//- Dune includes 
#include <dune/common/forloop.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/tuples.hh>

//- Dune fem includes 
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/io/file/iolock.hh>
#include <dune/fem/io/file/iointerface.hh>
#include <dune/fem/io/parameter.hh>

namespace Dune
{

  // IOTupleBase
  // -----------

  struct IOTupleBase 
  {
    // return path/ prefix name as one string   
    static std::string pathAndName(const std::string& path, const std::string& name, const std::string suffix) 
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

    template <class DataIO>
    static typename DataIO :: GridType*  
    restoreGrid(DataIO& dataio,
                double& t,int n,
                std::string path,
                std::string name) 
    {
      std::string gname ( gridName(path,name) );
      // check if lock file exists, and if exit 
      FileIOCheckError check( gname );
        
      // grid is created inside of restore grid method
      return dataio.restoreGrid( gname, t, n); 
    }

    template <class GridType>
    static void
    restoreDofManager ( const GridType& grid, int n,
                        const std::string &path, const std::string &name )
    {
      if( Parameter :: verbose() ) 
        std::cout << "Reading Dof Manager" << std::endl;
      
      typedef DofManager<GridType> DofManagerType;

      std::string dmname;
      dmname = gridName(path,name) + "_dm";
      DofManagerType& dm = DofManagerType :: instance(grid);

      if( Parameter :: verbose() ) 
        std::cout << "    from file " << dmname << std::endl;

      // read dofmanager, i.e. read all index sets 
      DofManagerType :: read(grid,dmname,n);
      
      // resize all data because size of index set might have changed  
      // NOTE: avoid resize of index sets by using resizeForRestict 
      dm.resizeForRestrict();
    }

    //! write grid and data to given directory 
    template <class DataIO,class GridType>
    static void 
    writeGrid(DataIO& dataio,GridType& grid,
              double t,int n,
              const std::string& path,
              const std::string& name) 
    {
      std::string gname( gridName( path, name ) );
      
      if( Parameter :: verbose() ) 
      {
        std::cout << "Writing grid to " << gname << std::endl;
      }

      // extra braces to destroy lock before data is written 
      {
        // create lock file which is removed after sucessful backup 
        FileIOLock lock( gname );
        
        // write grid 
        dataio.writeGrid(grid, xdr, gname, t, n);
      }
    }
  };



  // IOTuple
  // -------

  template< class Tuple >
  class IOTuple
  : public IOTupleBase 
  {
    template< int N > struct CreateData;
    template< int N > struct Restore;
    template< int N > struct Output;
    template< int N > struct AddToDisplay;
    template< int N > struct AddToDisplayOrRemove;
    template< int N > struct RemoveData;

    static const int length = tuple_size< Tuple >::value;

  public:
    template <class DataIO,class GridType>
    static Tuple *input ( DataIO &dataio, GridType *&grid, double &t, int n,
                          const std::string &path, const std::string &name )
    {
      // true if grid has to be read 
      const bool newGrid = (grid == 0);

      if( newGrid ) 
      {
        // create and read grid 
        grid = IOTupleBase::restoreGrid(dataio,t,n,path,name);
      }
      
      std::string dname( dataName(path,name) );
      std::cout << "Reading data from " << dname << std::endl;

      // create all data
      Tuple *ret = new Tuple;
      ForLoop< CreateData, 0, length-1 >::apply( dataio, dname, n, *grid, *ret );
      
      if( newGrid ) 
      {
        // now read dofmanager and index sets 
        IOTupleBase::restoreDofManager(*grid,n,path,name);
      }

      // now read all data 
      ForLoop< Restore, 0, length-1 >::apply( *ret, dataio, dname, n );
      
      if( newGrid ) 
      {
        typedef DofManager<GridType> DofManagerType;

        // get dof manager 
        DofManagerType& dm = DofManagerType :: instance(*grid);
     
        // compress all data 
        dm.compress();
      }
      
      std::cout << "    FINISHED!" << std::endl;
      return ret;
    }

    //! restore all data in tupel 
    template< class DataIO, class GridType >
    static void restoreData ( Tuple &data,
                              DataIO &dataio, const GridType &grid, int n,
                              const std::string &path,
                              const std::string &name ) 
    {
      std::string dname( dataName(path,name) );
      if( Parameter :: verbose() )
      {
        std::cout << "P["<< grid.comm().rank()<< "] Reading data from " << dname << std::endl;
      }

      // read dofmanager and index sets 
      IOTupleBase::restoreDofManager(grid,n,path,name);

      // read all data now 
      ForLoop< Restore, 0, length-1 >::apply( data, dataio, dname, n );

      typedef DofManager<GridType> DofManagerType;

      // get dof manager 
      DofManagerType& dm = DofManagerType :: instance(grid);
      
      // compress all data 
      dm.compress();

      if( Parameter :: verbose() ) 
      {
        std::cout << "P["<<grid.comm().rank()<< "]  FINISHED!" << std::endl;
      }
    }

    //! write grid and data to given directory 
    template< class DataIO, class GridType >
    static void output ( DataIO &dataio, GridType &grid,
                         double t, int n,
                         const std::string &path,
                         const std::string &name,
                         const Tuple &tuple )
    {
      // write grid first 
      writeGrid( dataio, grid, t, n, path, name );

      std::string dname( dataName(path, name ) );

      if( Parameter::verbose() ) 
        std::cout << "Writing data to " << dname << std::endl;
      
      // write data
      ForLoop< Output, 0, length-1 >::apply( dataio, dname, n, tuple );
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

    template< class DataIO, class Grid >
    static void apply ( DataIO &dataio, const std::string &name, const int &n, Grid &grid, Tuple &tuple )
    {
      GridPart *gridPart = new GridPart( grid );
      DiscreteFunctionSpace *space = new DiscreteFunctionSpace( *gridPart );
      get< N >( tuple ) = new DiscreteFunction( "", *space );
    }
  };


  template< class Tuple >
  template< int N >
  struct IOTuple< Tuple >::Restore
  {
    typedef typename TypeTraits< typename tuple_element< N, Tuple >::type >::PointeeType DiscreteFunction;

    template< class DataIO >
    static void apply ( const Tuple &tuple, DataIO &dataio, const std::string &name, const int &n )
    {
      std::stringstream dataname;
      dataname << name << "_" << N; 
      if( Parameter :: verbose() ) 
        std::cout << "    Dataset from " << dataname.str() << std::endl;

      DiscreteFunction *df = get< N >( tuple );
      assert( df );
      
      // check if lock file exists, and if exit
      FileIOCheckError check( dataname.str() );

      // make sure dof compression is applied 
      df->enableDofCompression();

      // read data 
      dataio.readData( *df, dataname.str(), n );
    }
  };


  template< class Tuple >
  template< int N >
  struct IOTuple< Tuple >::Output
  {
    typedef typename TypeTraits< typename tuple_element< N, Tuple >::type >::PointeeType DiscreteFunction;

    template< class DataIO >
    static void apply ( DataIO &dataio, const std::string &name, const int &n, const Tuple &tuple )
    {
      std::stringstream dataname;
      dataname << name << "_" << N;

      // create lock file which is removed after sucessful backup 
      FileIOLock lock( dataname.str() );

      const DiscreteFunction *df = get< N >( tuple );
      dataio.writeData( *df, xdr, dataname.str(), n );
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
      assert( dinf->comp );
      std::cout << "adding to display " << dinf->name << std::endl;
      disp.addData( *df, dinf, time );
    }

    template< class Disp >
    static void apply ( Disp &disp, Tuple &tuple )
    {
      DiscreteFunction *df = get< N >( tuple );
      disp.addData( *df );
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
    template <class DataIO,class GridType>
    static Tuple *input ( DataIO &dataio, GridType *&grid, double &t, int n,
                          const std::string &path, const std::string &name )
    {
      // true if grid has to be read 
      const bool newGrid = (grid == 0);

      if( newGrid ) 
      {
        // create and read grid 
        grid = IOTupleBase::restoreGrid(dataio,t,n,path,name);
      }
      
      std::string dname( dataName(path,name) );
      std::cout << "Reading data from " << dname << std::endl;

      // create all data
      Tuple *ret = new Tuple;

      if( newGrid ) 
      {
        // now read dofmanager and index sets 
        IOTupleBase::restoreDofManager(*grid,n,path,name);
      }

      if( newGrid ) 
      {
        typedef DofManager<GridType> DofManagerType;

        // get dof manager 
        DofManagerType& dm = DofManagerType :: instance(*grid);
     
        // compress all data 
        dm.compress();
      }
      
      std::cout << "    FINISHED!" << std::endl;
      return ret;
    }

    //! restore all data in tupel 
    template< class DataIO, class GridType >
    static void restoreData ( Tuple &data,
                              DataIO &dataio, const GridType &grid, int n,
                              const std::string &path,
                              const std::string &name ) 
    {
      std::string dname( dataName(path,name) );
      if( Parameter :: verbose() )
      {
        std::cout << "P["<< grid.comm().rank()<< "] Reading data from " << dname << std::endl;
      }

      // read dofmanager and index sets 
      IOTupleBase::restoreDofManager(grid,n,path,name);

      typedef DofManager<GridType> DofManagerType;

      // get dof manager 
      DofManagerType& dm = DofManagerType :: instance(grid);
      
      // compress all data 
      dm.compress();

      if( Parameter :: verbose() ) 
      {
        std::cout << "P["<<grid.comm().rank()<< "]  FINISHED!" << std::endl;
      }
    }

    //! write grid and data to given directory 
    template< class DataIO, class GridType >
    static void output ( DataIO &dataio, GridType &grid,
                         double t, int n,
                         const std::string &path,
                         const std::string &name,
                         const Tuple &tuple )
    {
      // write grid first 
      writeGrid( dataio, grid, t, n, path, name );

      std::string dname( dataName(path, name ) );

      if( Parameter::verbose() ) 
        std::cout << "Writing data to " << dname << std::endl;
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
