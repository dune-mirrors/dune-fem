#include <iostream>
#include <vector>
#include <assert.h>
#include <string.h>

#include <config.h>
#include <dune/common/stdstreams.cc>

using namespace Dune;

#ifndef DIM 
#define DIM 2
#endif
#ifndef DIM_OF_WORLD
#define DIM_OF_WORLD 2
#endif

#define AGRID 0
#define BGRID 0 
#define SGRID 1

//#include <dune/grid/sgrid.hh>

#define LARGE 1.0E308

#if AGRID 
#include <dune/grid/albertagrid.hh>

static const int dim = DIM; 
static const int dimworld = DIM_OF_WORLD; 

typedef AlbertaGrid<dim,dimworld> GR_GridType;
#endif

#if BGRID 
#include <dune/grid/alu3dgrid/includecc.cc>
#include <dune/grid/alu3dgrid.hh>

static const int dim = 3; 
static const int dimworld = 3; 

typedef ALU3dGrid<dim,dimworld,hexa> GR_GridType;
#endif


#if SGRID 
#include <dune/grid/sgrid.hh>

static const int dim = DIM; 
static const int dimworld = DIM_OF_WORLD; 

typedef SGrid <dim, dimworld> GR_GridType;
#endif

#include <dune/fem/dfadapt.hh>
#include <dune/fem/lagrangebase.hh>
#include <dune/fem/dfadapt.hh>
#include <dune/common/stack.hh>

#include <dune/grid/common/leafindexset.hh>
#include <dune/grid/common/gridpart.hh>

#include <dune/io/file/grapedataio.hh>


typedef double REAL;

#include <dune/io/visual/grapedatadisplay.hh>
#include <dune/io/visual/combinedgrapedisplay.hh>
#include "printhelp.cc"

//typedef FunctionSpace <double ,double , dim, dim+2 >  GR_FunctionSpaceType;

typedef FunctionSpace <double ,double , dim, 1 >        GR_FunctionSpaceType;


typedef DofManager<GR_GridType>                         GR_DofManagerType;
typedef DofManagerFactory <GR_DofManagerType>           GR_DofManagerFactoryType;

#if SGRID 
  typedef DefaultGridIndexSet<GR_GridType,LevelIndex>   GR_IndexSetType;
#else
  typedef GR_GridType :: Traits :: LeafIndexSet         GR_IndexSetType;
#endif


typedef DefaultGridPart<GR_GridType,GR_IndexSetType>    GR_GridPartType;

typedef LagrangeDiscreteFunctionSpace<GR_FunctionSpaceType,
				     GR_GridPartType,0> GR_DiscFuncSpaceType;

typedef DFAdapt < GR_DiscFuncSpaceType >                GR_DiscFuncType;
typedef GrapeDataDisplay<GR_GridType , GR_DiscFuncType> GrapeDispType;

#include "readdata.cc"
#include "readparams.cc" 
#include "readfile.cc" 


int main(int argc, char **argv)
{

  if (argc < 2)
  {
    print_help("datadisp");
    return(0);
  }   

  if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "-help"))
  {
      print_help("dunedisp");
      return(0);
  }



  if(argc > 2) 
    return readParameterList(argc,argv);
  else 
    return readParameterFile(argc,argv);
 
}
