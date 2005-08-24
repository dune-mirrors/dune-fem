#include <iostream>
#include <vector>
#include <assert.h>
#include <string.h>

#include <config.h>
#include <dune/common/stdstreams.cc>

using namespace Dune;

#ifndef DIM 
#define DIM 3
#endif
#ifndef DIM_OF_WORLD
#define DIM_OF_WORLD 3
#endif

#define AGRID 0 
#define BGRID 1 
#define SGRID 0

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

typedef ALU3dGrid<dim,dimworld,tetra> GR_GridType;
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

typedef FunctionSpace <double ,double , dim, dim+2 >  GR_FunctionSpaceType;
typedef DofManager<GR_GridType,DataCollectorInterface<GR_GridType,GR_GridType::ObjectStreamType> > GR_DofManagerType;
typedef DofManagerFactory <GR_DofManagerType> GR_DofManagerFactoryType;

//typedef GR_GridType :: LeafIndexSetType GR_IndexSetType;
//typedef DefaultGridIndexSet<GR_GridType, GlobalIndex > GR_IndexSetType;
//typedef AdaptiveLeafIndexSet<GR_GridType> GR_IndexSetType;

typedef LeafGridPart<GR_GridType> GR_GridPartType;
typedef GR_GridPartType::IndexSetType GR_IndexSetType;

typedef LagrangeDiscreteFunctionSpace<GR_FunctionSpaceType,GR_GridPartType,0, GR_DofManagerType> GR_DiscFuncSpaceType;

typedef DFAdapt < GR_DiscFuncSpaceType > GR_DiscFuncType;
typedef GrapeDataDisplay<GR_GridType , GR_DiscFuncType > GrapeDispType;

#include "readdata.cc"

int main(int argc, char **argv)
{
  int   i, i_start, i_end;
  INFO * info = 0;
  int    n = 0, n_info = 10;

  int    i_delta = 1;
  const  char *path = 0;
  bool   time_bar = false;
  int    parallel = 1;
  
  REAL   timestep = 1.0e-3;

  info = (INFO *) malloc(n_info*sizeof(INFO));
  assert(info != 0);
  
  info[0].datinf = 0;
  info[0].name = "grid";
  info[0].fix_mesh = 0;

  if (argc < 3)
  {
    print_help("dunedisp");
    return(0);
  }

  i_start = atoi(argv[1]);
  i_end = atoi(argv[2]);

  i = 3;
  while (i < argc)
  {
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help"))
    {
      print_help("dunedisp");
    }
    else if (!strcmp(argv[i], "-i"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -i `increment'\n");
      i_delta = atoi(argv[i+1]);
      i += 2;
    }
    else if (!strcmp(argv[i], "-s"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -s `dataprefix'\n");
      
      DATAINFO * dinf = (DATAINFO *) std::malloc(sizeof(DATAINFO));
      assert(dinf);
      dinf->name = argv[i+1];
      dinf->vector = 0;
      /* seems wrong order, but grape truns it arround, we can do nothing else here */
      dinf->next = info[n].datinf; 
      info[n].datinf = dinf;

      i += 2;
    }
    else if (!strcmp(argv[i], "-v"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -v `vectorprefix'\n");
      
      DATAINFO * dinf = (DATAINFO *) std::malloc(sizeof(DATAINFO));
      assert(dinf);
      dinf->name = argv[i+1];
      dinf->vector = 1;
      /* seems wrong order, but grape truns it arround, we can do nothing else here */
      dinf->next = info[n].datinf; 
      info[n].datinf = dinf;

      i += 2;
    }
    else if (!strcmp(argv[i], "-t"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -t `time step size'\n");
      timestep = atof(argv[i+1]);
      i += 2;
    }
    else if (!strcmp(argv[i], "-m"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -m `gridprefix'\n");
      assert(n < n_info);
      info[n].name = argv[i+1];
      info[n].datinf = 0;
      info[n].fix_mesh = 0;
      n++;
      i += 2;
    }
    else if (!strcmp(argv[i], "-b"))
    {
      time_bar = true;
      i += 1;
    }
    else if (!strcmp(argv[i], "-f"))
    {
      info[n].fix_mesh = 1;
      i += 1;
    }
    else if (!strcmp(argv[i], "-pg"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -pg `number of procs'\n");
      parallel += atoi(argv[i+1]);
      i += 2;
    }
    else if (!strcmp(argv[i], "-p"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -p `path'\n");
      path = argv[i+1];
      i += 2;
    }
    else
    {
      std::cerr << "unknown option " << argv[i] << std::endl;
      exit(1);
    }
    printf("i = %d, argc = %d\n", i, argc);
  }

  timeSceneInit(info, n , parallel , time_bar);
  readData(info, path,i_start,i_end,i_delta,n,timestep,parallel);
  
  // run grape 
  displayTimeScene(info,parallel);

  deleteAllObjects();
  return (EXIT_SUCCESS);
}
