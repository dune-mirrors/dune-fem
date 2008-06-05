#ifndef DUNE_INPUTOUTPUTTUPLES_HH
#define DUNE_INPUTOUTPUTTUPLES_HH

//- system includes 
#include <string>

//- Dune includes 
#include <dune/common/typetraits.hh>

//- Dune fem includes 
#include <dune/fem/misc/femtuples.hh>
#include <dune/fem/misc/utility.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/io/file/iolock.hh>
#include <dune/fem/io/file/iointerface.hh>

namespace Dune {

template <int N,class DiscFuncType>
struct IOTupleCaller 
{
  template <class DataIO,class GridType>
  static DiscFuncType* createData(DataIO& dataio,std::string name,int n,
           GridType& grid) 
  {
    typedef typename DiscFuncType::DiscreteFunctionSpaceType SpaceType;
    typedef typename SpaceType::GridPartType GridPartType;
    
    GridPartType* gridPart = new GridPartType(grid);
    SpaceType* space = new SpaceType(*gridPart);
    DiscFuncType* df = new DiscFuncType ("",*space);
    return df;
  }
  
  static void removeData(DiscFuncType * df)
  {
    typedef typename DiscFuncType::DiscreteFunctionSpaceType SpaceType;
    typedef typename SpaceType::GridPartType GridPartType;
    
    const GridPartType* gridPart = &(df->space().gridPart());
    const SpaceType* space = &(df->space());
    delete df;
    delete space;
    delete gridPart;
  }
  
  template <class DataIO>
  static void restore(DiscFuncType* df, 
                      DataIO& dataio,std::string name,int n) 
  {
    std::stringstream dataname;
    dataname << name << "_" << N; 
    std::cout << "    Dataset from " << dataname.str() << std::endl;
    assert( df );
    
    // check if lock file exists, and if exit 
    FileIOCheckError check( dataname.str() );

    dataio.readData(*df, dataname.str().c_str(), n);
  }
  
  template <class DataIO>
  static void output(DataIO& dataio,std::string name,int n,
         const DiscFuncType& df) 
  {
    std::stringstream dataname;
    dataname << name << "_" << N;
    
    // create lock file which is removed after sucessful backup 
    FileIOLock lock( dataname.str() );
      
    dataio.writeData(df, xdr, dataname.str().c_str(), n);
  }
  
  template <class Disp,class DINFO>
  static void addToDisplay(Disp& disp,const DINFO* dinf,double time,
         DiscFuncType& df) 
  {
    assert( dinf->comp );
    std::cout << "adding to display " << dinf->name << std::endl;
    disp.addData(df,dinf,time);
  }
  
  template <class Disp>
  static void addToDisplay(Disp& disp, DiscFuncType& df) 
  {
    disp.addData(df);
  }
};


template <class T1,class T2,int N>
struct IOTupleHelper 
{
  typedef Pair<T1*,T2> ReturnType;
  typedef typename TypeTraits<typename T2::Type1>::PointeeType T21;
  typedef Pair<T1*,T2> ThisType;
  typedef IOTupleHelper<T21,typename T2::Type2,N+1> NextType;
  
  template <class DataIO,class GridType>
  static ReturnType createData(DataIO& dataio,std::string name,int n,
        GridType& grid) {
    T2 next = NextType::createData(dataio,name,n,grid);
    T1* df = IOTupleCaller<N,T1>::createData(dataio,name,n,grid);
    return ReturnType(df,next);
  }
  
  template <class DataIO>
  static void 
  restore(ReturnType ret,
          DataIO& dataio, std::string name, int n)
  {
    T2 next = ret.second();
    NextType::restore(next,dataio,name,n);
    
    T1* df = ret.first();
    IOTupleCaller<N,T1>::restore(df,dataio,name,n);
  }
  
  template <class DataIO>
  static void output(DataIO& dataio,std::string name,int n,
         const ThisType& tup) {
    IOTupleCaller<N,T1>::output(dataio,name,n,*(tup.first()));
    NextType::output(dataio,name,n,tup.second());
  }
  template <class Disp,class DINFO>
  static void addToDisplay(Disp& disp,const DINFO* dinf,double time,
         ThisType& tup) {
    NextType::addToDisplay(disp,dinf->next,time,tup.second());
    IOTupleCaller<N,T1>::addToDisplay(disp,dinf,time,*(tup.first()));
  }
  template <class Disp,class DINFO>
  static void addToDisplayOrRemove(Disp& disp,const DINFO* dinf,double time,
         ThisType& tup) 
  {
    NextType::addToDisplayOrRemove(disp,dinf->next,time,tup.second());

    // if comp is zero data set is not valid 
    if( dinf->comp )
    {
      IOTupleCaller<N,T1>::addToDisplay(disp,dinf,time,*(tup.first()));
    }
    else 
    {
      IOTupleCaller<N,T1>::removeData(tup.first());
      tup = ThisType(0,tup.second());
    }
  }
  template <class Disp>
  static void addToDisplay(Disp& disp, ThisType& tup) 
  {
    NextType::addToDisplay(disp,tup.second());
    IOTupleCaller<N,T1>::addToDisplay(disp,*(tup.first()));
  }
  
  static void removeData(ThisType& tup) 
  {
    NextType::removeData(tup.second());
    IOTupleCaller<N,T1>::removeData(tup.first());
  }
};

template <class T1,int N>
struct IOTupleHelper<T1,Nil,N> 
{
  typedef Pair<T1*,Nil> ReturnType;
  typedef Pair<T1*,Nil> ThisType;
  
  template <class DataIO,class GridType>
  static ReturnType createData(DataIO& dataio,std::string name,int n,
        GridType& grid) {
    return ReturnType(IOTupleCaller<N,T1>::
                createData(dataio,name,n,grid),nullType());
  }
  
  template <class DataIO>
  static void restore(ReturnType ret, 
        DataIO& dataio,std::string name,int n)
  {
    T1 * df = ret.first();
    IOTupleCaller<N,T1>::restore(df,dataio,name,n);
  }

  template <class DataIO>
  static void output(DataIO& dataio,std::string name,int n,
         const ThisType& tup) {
    IOTupleCaller<N,T1>::output(dataio,name,n,*(tup.first()));
  }
  template <class Disp,class DINFO>
  static void addToDisplay(Disp& disp,const DINFO* dinf,double time,
         ThisType& tup) {
    IOTupleCaller<N,T1>::addToDisplay(disp,dinf,time,*(tup.first()));
  }

  template <class Disp,class DINFO>
  static void addToDisplayOrRemove(
      Disp& disp,const DINFO* dinf,double time, ThisType& tup) 
  {
    // if comp is zero data set is not valid 
    if( dinf->comp )
    {
      IOTupleCaller<N,T1>::addToDisplay(disp,dinf,time,*(tup.first()));
    }
    else 
    {
      IOTupleCaller<N,T1>::removeData(tup.first());
      tup = ThisType(0,Nil());
    }
  }
  
  template <class Disp>
  static void addToDisplay(Disp& disp,ThisType& tup) 
  {
    IOTupleCaller<N,T1>::addToDisplay(disp,*(tup.first()));
  }
  
  static void removeData(ThisType& tup) 
  {
    IOTupleCaller<N,T1>::removeData(tup.first());
  }
};

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
};

template <class TupType>
struct IOTuple : public IOTupleBase 
{
  typedef typename TypeTraits<typename TupType::Type1>::PointeeType T1;
  typedef typename TupType::Type2 T2;
  typedef typename IOTupleHelper<T1,T2,0>::ReturnType ReturnType;

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
  static void restoreDofManager(const GridType& grid,
                                int n,
                                std::string path,
                                std::string name) 
  {
    std::cout << "Reading Dof Manager" << std::endl;
    typedef DofManager<GridType> DofManagerType;
    typedef DofManagerFactory<DofManagerType> DMFactoryType;
    std::string dmname;
    dmname = gridName(path,name) + "_dm";
    DofManagerType& dm = DMFactoryType::getDofManager(grid);
    std::cout << "    from file " << dmname << std::endl;
    // read dofmanager, i.e. read all index sets 
    DMFactoryType::readDofManager(grid,dmname,n);
    
    // resize all data because size of index set might have changed  
    // NOTE: avoid resize of index sets by using resizeForRestict 
    dm.resizeForRestrict();
  }
  
  template <class DataIO,class GridType>
  static ReturnType* input(DataIO& dataio,GridType*& grid,double& t,int n,
                           std::string path,
                           std::string name) 
  {
    // true if grid has to be read 
    const bool newGrid = (grid == 0);

    if( newGrid ) 
    {
      // create and read grid 
      grid = IOTuple<TupType>::restoreGrid(dataio,t,n,path,name);
    }
    
    std::string dname( dataName(path,name) );
    std::cout << "Reading data from " << dname << std::endl;

    // create all data
    ReturnType* ret =
      new ReturnType(IOTupleHelper<T1,T2,0>::createData(dataio,dname,n,*grid));
    
    if( newGrid ) 
    {
      // now read dofmanager and index sets 
      IOTuple<TupType>::restoreDofManager(*grid,n,path,name);
    }

    // now read all data 
    IOTupleHelper<T1,T2,0>::restore(*ret,dataio,dname,n);
    
    typedef DofManager<GridType> DofManagerType;
    typedef DofManagerFactory<DofManagerType> DMFactoryType;

    if( newGrid ) 
    {
      // get dof manager 
      DofManagerType& dm = DMFactoryType::getDofManager(*grid);
   
      // compress all data 
      dm.compress();
    }
    
    std::cout << "    FINISHED!" << std::endl;
    return ret;
  }

  //! restore all data in tupel 
  template <class DataIO,class GridType>
  static void restoreData(ReturnType& data, 
         DataIO& dataio,const GridType& grid,
         int n, std::string path, std::string name) 
  {
    std::string dname( dataName(path,name) );
    std::cout << "Reading data from " << dname << std::endl;
   
    // read dofmanager and index sets 
    IOTuple<TupType>::restoreDofManager(grid,n,path,name);

    // read all data now 
    IOTupleHelper<T1,T2,0>::restore(data,dataio,dname,n);

    typedef DofManager<GridType> DofManagerType;
    typedef DofManagerFactory<DofManagerType> DMFactoryType;

    // get dof manager 
    DofManagerType& dm = DMFactoryType::getDofManager(grid);
    
    // compress all data 
    dm.compress();
    
    std::cout << "    FINISHED!" << std::endl;
  }

  //! write grid and data to given directory 
  template <class DataIO,class GridType>
  static void output(DataIO& dataio,GridType& grid,
                     double t,int n,
                     std::string path,
                     std::string name, 
                     const Pair<T1*,T2>& tup, bool verbose = true ) 
  {
    std::string gname( gridName( path, name ) );
    
    if(verbose) 
    {
      std::cout << "Writing grid to " << gname << std::endl;
    }

    // extra braces to destroy lock before data is written 
    {
      // create lock file which is removed after sucessful backup 
      FileIOLock lock( gname );
      
      // write grid 
      dataio.writeGrid(grid, xdr, gname.c_str(), t, n);
    }

    std::string dname( dataName(path, name ) );

    if(verbose) 
    {
      std::cout << "Writing data to " << dname << std::endl;
    }
    
    // write data 
    IOTupleHelper<T1,T2,0>::output(dataio,dname,n,tup);
  }
  
  template <class Disp,class DINFO>
  static void addToDisplay(Disp& disp,const DINFO* dinf,double time,
                           Pair<T1*,T2>& tup) 
  {
    IOTupleHelper<T1,T2,0>::addToDisplay(disp,dinf,time,tup);
  }

  template <class Disp,class DINFO>
  static void addToDisplayOrRemove(Disp& disp,const DINFO* dinf,double time,
                                   Pair<T1*,T2>& tup) 
  {
    IOTupleHelper<T1,T2,0>::addToDisplayOrRemove(disp,dinf,time,tup);
  }

  template <class Disp>
  static void addToDisplay(Disp& disp, Pair<T1*,T2>& tup) 
  {
    IOTupleHelper<T1,T2,0>::addToDisplay(disp,tup);
  }
  
  static void removeData(Pair<T1*,T2>& tup) 
  {
    IOTupleHelper<T1,T2,0>::removeData(tup);
  }

};

} // end namespace Dune 
#endif
