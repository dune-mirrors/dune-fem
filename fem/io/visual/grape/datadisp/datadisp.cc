//************************************************************
//
//  (C) written and directed by Robert Kloefkorn 
//
//************************************************************
#include <iostream>
#include <vector>
#include <assert.h>
#include <string.h>

#include <config.h>
#include <dune/common/stdstreams.cc>

using namespace Dune;

#define AGRID 1
#define BGRID 0
#define SGRID 0

#include <dune/grid/sgrid.hh>

#define LARGE 1.0E308

#if AGRID 
#include <dune/grid/albertagrid.hh>

static const int dim = DUNE_PROBLEM_DIM; 
static const int dimworld = DUNE_WORLD_DIM; 

typedef AlbertaGrid<dim,dimworld> GR_GridType;
#endif

#if BGRID 

//#include <dune/grid/alu3dgrid/includecc.cc>
//#include <dune/grid/alu3dgrid.hh>
#include <dune/grid/alugrid.hh>

static const int dim = 3; 
static const int dimworld = 3; 

typedef ALU3dGrid<dim,dimworld,hexa> GR_GridType;
#endif


#if SGRID 
#include <dune/grid/sgrid.hh>

static const int dim = 2; 
static const int dimworld = 2; 

typedef SGrid <dim, dimworld> GR_GridType;
#endif

#include <dune/fem/discretefunction/dfadapt.hh>
#include <dune/fem/space/lagrangespace/lagrange.hh>
#include <dune/common/stack.hh>

#include <dune/grid/common/gridpart.hh>

//#include <dune/io/file/grapedataio.hh>


#include "../../../space/dgspace.hh"

typedef double REAL;


#include <dune/grid/io/visual/grapedatadisplay.hh>
#include <dune/grid/io/visual/combinedgrapedisplay.hh>
#include "printhelp.cc"

//typedef FunctionSpace <double ,double , dim, dim+2 >  GR_FunctionSpaceType;
typedef FunctionSpace <double ,double , dim, 1 >  GR_FunctionSpaceType;


typedef DofManager<GR_GridType>                         GR_DofManagerType;
typedef DofManagerFactory <GR_DofManagerType>           GR_DofManagerFactoryType;

//typedef DefaultGridPart<GR_GridType,GR_IndexSetType>    GR_GridPartType;
typedef LeafGridPart < GR_GridType > GR_GridPartType; 
typedef GR_GridPartType :: IndexSetType GR_IndexSetType;


//typedef DiscontinuousGalerkinSpace<GR_FunctionSpaceType, GR_GridPartType, 0> GR_DiscFuncSpaceType;
typedef LagrangeDiscreteFunctionSpace<GR_FunctionSpaceType,GR_GridPartType,0> GR_DiscFuncSpaceType;

typedef DFAdapt < GR_DiscFuncSpaceType >                GR_DiscFuncType;
typedef GrapeDataDisplay<GR_GridType> GrapeDispType;

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
