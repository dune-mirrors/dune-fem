#include <string>
#include <dune/common/typetraits.hh>
#include <fem/pass/tuples.hh>
#include <dune/fem/space/common/dofmanager.hh>
namespace Dune {
template <int N,class DiscFuncType>
struct GrapeTupleCaller {
  template <class DataIO,class GridType>
  static DiscFuncType* input(DataIO& dataio,std::string name,int n,
			     GridType& grid) {
    std::stringstream dataname;
    dataname << name << "_" << N; 
    std::cout << "    Dataset from " << dataname.str() << std::endl;
    typedef typename DiscFuncType::DiscreteFunctionSpaceType SpaceType;
    typedef typename SpaceType::GridPartType GridPartType;
    GridPartType* gridPart = new GridPartType(grid);
    SpaceType* space = new SpaceType(*gridPart);
    DiscFuncType* df = new DiscFuncType (dataname.str().c_str(), *space);
    dataio.readData(*df, dataname.str().c_str(), n);
    return df;
  }
  template <class DataIO>
  static void output(DataIO& dataio,std::string name,int n,
		     DiscFuncType& df) {
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
};
template <class T1,class T2,int N>
struct GrapeTupleHelper {
  typedef Pair<T1*,T2> ReturnType;
  typedef typename TypeTraits<typename T2::Type1>::PointeeType T21;
  typedef Pair<T1*,T2> ThisType;
  typedef GrapeTupleHelper<T21,typename T2::Type2,N+1> NextType;
  template <class DataIO,class GridType>
  static ReturnType input(DataIO& dataio,std::string name,int n,
			  GridType& grid) {
    T2 next = NextType::input(dataio,name,n,grid);
    T1* df = GrapeTupleCaller<N,T1>::input(dataio,name,n,grid);
    return ReturnType(df,next);
  }
  template <class DataIO>
  static void output(DataIO& dataio,std::string name,int n,
		     ThisType& tup) {
    GrapeTupleCaller<N,T1>::output(dataio,name,n,*(tup.first()));
    NextType::output(dataio,name,n,tup.second());
  }
  template <class Disp,class DINFO>
  static void addToDisplay(Disp& disp,const DINFO* dinf,double time,
			   ThisType& tup) {
    NextType::addToDisplay(disp,dinf->next,time,tup.second());
    GrapeTupleCaller<N,T1>::addToDisplay(disp,dinf,time,*(tup.first()));
  }
};
template <class T1,int N>
struct GrapeTupleHelper<T1,Nil,N> {
  typedef Pair<T1*,Nil> ReturnType;
  typedef Pair<T1*,Nil> ThisType;
  template <class DataIO,class GridType>
  static ReturnType input(DataIO& dataio,std::string name,int n,
			  GridType& grid) {
    return ReturnType(GrapeTupleCaller<N,T1>::input(dataio,name,n,grid),nullType());
  }
  template <class DataIO>
  static void output(DataIO& dataio,std::string name,int n,
		     ThisType& tup) {
    GrapeTupleCaller<N,T1>::output(dataio,name,n,*(tup.first()));
  }
  template <class Disp,class DINFO>
  static void addToDisplay(Disp& disp,const DINFO* dinf,double time,
			   ThisType& tup) {
    GrapeTupleCaller<N,T1>::addToDisplay(disp,dinf,time,*(tup.first()));
  }
};
template <class TupType>
struct GrapeTuple {
  typedef typename TypeTraits<typename TupType::Type1>::PointeeType T1;
  typedef typename TupType::Type2 T2;
  typedef typename GrapeTupleHelper<T1,T2,0>::ReturnType ReturnType;
  template <class DataIO,class GridType>
  static ReturnType* input(DataIO& dataio,GridType*& grid,double& t,int n,
			   const char* path,
			   std::string name) {
    std::string gname;
    if (path) gname += path;
    else gname += ".";
    gname += "/g";
    gname += name;
    std::cout << "Reading grid from " << gname << std::endl;
    grid = new GridType();
    dataio.readGrid(*grid, gname.c_str(), t, n);
    std::string dname;
    if (path) dname += path;
    else dname += ".";
    dname += "/d";
    dname += name;
    std::cout << "Reading data from " << dname << std::endl;
    ReturnType* ret =
      new ReturnType(GrapeTupleHelper<T1,T2,0>::input(dataio,dname,n,*grid));
    std::cout << "Reading Dof Manager" << std::endl;
    typedef DofManager<GridType> DofManagerType;
    typedef DofManagerFactory<DofManagerType> DMFactoryType;
    std::string dmname;
    dmname = gname + "_dm";
    DMFactoryType::getDofManager(*grid);
    std::cout << "    from file " << dmname << std::endl;
    DMFactoryType::readDofManager(*grid,dmname,n);
    std::cout << "    FINISHED!" << std::endl;
    return ret;
  }
  template <class DataIO,class GridType>
  static void output(DataIO& dataio,GridType& grid,double t,int n,
		     const char* path,
		     std::string name, Pair<T1*,T2>& tup) {
    std::string gname;
    gname += name;
    gname += "/g";
    if (path) gname += path;
    else gname += ".";
    std::cout << "Writing grid to " << gname << std::endl;
    dataio.writeGrid(grid, xdr, gname.c_str(), t, n);
    std::string dname;
    dname += name;
    dname += "/d";
    if (path) dname += path;
    else dname += ".";
    std::cout << "Writing data to " << dname << std::endl;
    GrapeTupleHelper<T1,T2,0>::output(dataio,dname,n,tup);
  }
  template <class Disp,class DINFO>
  static void addToDisplay(Disp& disp,const DINFO* dinf,double time,
			   Pair<T1*,T2>& tup) {
    GrapeTupleHelper<T1,T2,0>::addToDisplay(disp,dinf,time,tup);
  }
};

};
