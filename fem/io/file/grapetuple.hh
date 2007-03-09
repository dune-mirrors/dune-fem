#ifndef DUNE_GRAPETUPLES_HH
#define DUNE_GRAPETUPLES_HH

//- system includes 
#include <string>

//- Dune includes 
#include <dune/common/typetraits.hh>
#include <dune/common/tuples.hh>

//- Dune fem includes 
#include <dune/fem/pass/utility.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune {

template <int N,class DiscFuncType>
struct GrapeTupleCaller 
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
  
  template <class DataIO>
  static void restore(DiscFuncType* df, 
                      DataIO& dataio,std::string name,int n) 
  {
    std::stringstream dataname;
    dataname << name << "_" << N; 
    std::cout << "    Dataset from " << dataname.str() << std::endl;
    assert( df );
    dataio.readData(*df, dataname.str().c_str(), n);
  }
  
  template <class DataIO>
  static void output(DataIO& dataio,std::string name,int n,
		     const DiscFuncType& df) {
    std::stringstream dataname;
    dataname << name << "_" << N;
    dataio.writeData(df, xdr, dataname.str().c_str(), n);
  }
  template <class Disp,class DINFO>
  static void addToDisplay(Disp& disp,const DINFO* dinf,double time,
			   DiscFuncType& df) {
    std::cout << "adding to display " << df.name() << std::endl;
    disp.addData(df,dinf,time);
  }
  
  template <class Disp>
  static void addToDisplay(Disp& disp, DiscFuncType& df) 
  {
    disp.addData(df);
  }
};


template <class T1,class T2,int N>
struct GrapeTupleHelper 
{
  typedef Pair<T1*,T2> ReturnType;
  typedef typename TypeTraits<typename T2::Type1>::PointeeType T21;
  typedef Pair<T1*,T2> ThisType;
  typedef GrapeTupleHelper<T21,typename T2::Type2,N+1> NextType;
  
  template <class DataIO,class GridType>
  static ReturnType createData(DataIO& dataio,std::string name,int n,
			  GridType& grid) {
    T2 next = NextType::createData(dataio,name,n,grid);
    T1* df = GrapeTupleCaller<N,T1>::createData(dataio,name,n,grid);
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
    GrapeTupleCaller<N,T1>::restore(df,dataio,name,n);
  }
  
  template <class DataIO>
  static void output(DataIO& dataio,std::string name,int n,
		     const ThisType& tup) {
    GrapeTupleCaller<N,T1>::output(dataio,name,n,*(tup.first()));
    NextType::output(dataio,name,n,tup.second());
  }
  template <class Disp,class DINFO>
  static void addToDisplay(Disp& disp,const DINFO* dinf,double time,
			   ThisType& tup) {
    NextType::addToDisplay(disp,dinf->next,time,tup.second());
    GrapeTupleCaller<N,T1>::addToDisplay(disp,dinf,time,*(tup.first()));
  }
  template <class Disp>
  static void addToDisplay(Disp& disp, ThisType& tup) 
  {
    NextType::addToDisplay(disp,tup.second());
    GrapeTupleCaller<N,T1>::addToDisplay(disp,*(tup.first()));
  }
};


template <class T1,int N>
struct GrapeTupleHelper<T1,Nil,N> 
{
  typedef Pair<T1*,Nil> ReturnType;
  typedef Pair<T1*,Nil> ThisType;
  
  template <class DataIO,class GridType>
  static ReturnType createData(DataIO& dataio,std::string name,int n,
			  GridType& grid) {
    return ReturnType(GrapeTupleCaller<N,T1>::
                createData(dataio,name,n,grid),nullType());
  }
  
  template <class DataIO>
  static void restore(ReturnType ret, 
        DataIO& dataio,std::string name,int n)
  {
    T1 * df = ret.first();
    GrapeTupleCaller<N,T1>::restore(df,dataio,name,n);
  }

  template <class DataIO>
  static void output(DataIO& dataio,std::string name,int n,
		     const ThisType& tup) {
    GrapeTupleCaller<N,T1>::output(dataio,name,n,*(tup.first()));
  }
  template <class Disp,class DINFO>
  static void addToDisplay(Disp& disp,const DINFO* dinf,double time,
			   ThisType& tup) {
    GrapeTupleCaller<N,T1>::addToDisplay(disp,dinf,time,*(tup.first()));
  }
  template <class Disp>
  static void addToDisplay(Disp& disp,ThisType& tup) 
  {
    GrapeTupleCaller<N,T1>::addToDisplay(disp,*(tup.first()));
  }
};


template <class TupType>
struct GrapeTuple 
{
  typedef typename TypeTraits<typename TupType::Type1>::PointeeType T1;
  typedef typename TupType::Type2 T2;
  typedef typename GrapeTupleHelper<T1,T2,0>::ReturnType ReturnType;

  // return path/ prefix name as one string   
  static std::string pathAndName(const std::string& path, const std::string& name, const std::string prefix) 
  {
    std::string comboname;
    if (path != "") comboname += path;
    else comboname += ".";
    comboname += "/";
    comboname += prefix;
    comboname += name;
    return comboname;
  }
  
  // return grids name  
  static std::string gridName(const std::string& path, const std::string& name) 
  {
    return pathAndName(path,name,"g");
  }
  
  // return datas name  
  static std::string dataName(const std::string& path, const std::string& name) 
  {
    return pathAndName(path,name,"d");
  }

  template <class DataIO>
  static typename DataIO :: GridType*  
  restoreGrid(DataIO& dataio,
              double& t,int n,
			        std::string path,
			        std::string name) 
  {
    // grid is created inside of restore grid method
    return dataio.restoreGrid( gridName(path,name) , t, n); 
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
    DMFactoryType::getDofManager(grid);
    std::cout << "    from file " << dmname << std::endl;
    DMFactoryType::readDofManager(grid,dmname,n);
  }
  
  template <class DataIO,class GridType>
  static ReturnType* input(DataIO& dataio,GridType*& grid,double& t,int n,
			   std::string path,
			   std::string name) 
  {
    // read grid 
    grid = GrapeTuple<TupType>::restoreGrid(dataio,t,n,path,name);
    
    std::string dname( dataName(path,name) );
    std::cout << "Reading data from " << dname << std::endl;

    // create all data
    ReturnType* ret =
      new ReturnType(GrapeTupleHelper<T1,T2,0>::createData(dataio,dname,n,*grid));
    
    // read all data 
    GrapeTupleHelper<T1,T2,0>::restore(*ret,dataio,dname,n);
   
    // read dofmanager and index sets 
    GrapeTuple<TupType>::restoreDofManager(*grid,n,path,name);
    
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
   
    // read all data 
    GrapeTupleHelper<T1,T2,0>::restore(data,dataio,dname,n);
    
    // read dofmanager and index sets 
    GrapeTuple<TupType>::restoreDofManager(grid,n,path,name);

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

    // write grid 
    dataio.writeGrid(grid, xdr, gname.c_str(), t, n);

    std::string dname( dataName(path, name ) );

    if(verbose) 
    {
      std::cout << "Writing data to " << dname << std::endl;
    }
    
    // write data 
    GrapeTupleHelper<T1,T2,0>::output(dataio,dname,n,tup);
  }
  
  template <class Disp,class DINFO>
  static void addToDisplay(Disp& disp,const DINFO* dinf,double time,
			   Pair<T1*,T2>& tup) 
  {
    GrapeTupleHelper<T1,T2,0>::addToDisplay(disp,dinf,time,tup);
  }

  template <class Disp>
  static void addToDisplay(Disp& disp, Pair<T1*,T2>& tup) 
  {
    GrapeTupleHelper<T1,T2,0>::addToDisplay(disp,tup);
  }
};

} // end namespace Dune 
#endif
