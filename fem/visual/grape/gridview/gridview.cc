#include <iostream>

#include "config.h"

#ifndef DIM 
#define DIM 2 
#endif

#ifndef DIM_OF_WORLD
#define DIM_OF_WORLD 2
#endif

//#define UGGRID_WITH_INDEX_SETS

#define PROBLEMDIM DIM 
#define WORLDDIM DIM_OF_WORLD

#define SGRID 0
#define UGRID 0
#define AGRID 1

#define BGRID 0

#define YGRID 0

#if SGRID
#include <dune/grid/sgrid.hh>
using namespace Dune;

#undef DIM 
#undef DIM_OF_WORLD

static const int DIM = 3;
static const int DIM_OF_WORLD = 3;


typedef SGrid <DIM,DIM_OF_WORLD> GridType;
#endif

#if UGRID

#include <dune/grid/uggrid.hh>
#include <dune/io/file/amirameshreader.hh>
#include <dune/io/file/amirameshwriter.hh>

#undef DIM 
#undef DIM_OF_WORLD

#ifdef _2 
static const int DIM = 2;
static const int DIM_OF_WORLD = 2;
#else 
static const int DIM = 3;
static const int DIM_OF_WORLD = 3;
#endif
using namespace Dune;
typedef UGGrid< DIM,DIM_OF_WORLD> GridType;
#endif

#if YGRID
#include <dune/grid/yaspgrid.hh>
using namespace Dune;
typedef YaspGrid <DIM,DIM_OF_WORLD> GridType;
#endif

#if AGRID
#include <dune/grid/albertagrid.hh>
using namespace Dune;
typedef AlbertaGrid<DIM,DIM_OF_WORLD> GridType;
#define _2  
#endif

#if BGRID

#undef DIM
#undef DIM_OF_WORLD
#undef PROBLEMDIM 
#undef WORLDDIM

#define DIM 3 
#define DIM_OF_WORLD 3

#define PROBLEMDIM 3
#define WORLDDIM 3

#include <dune/grid/alu3dgrid.hh>
#include <dune/grid/alu3dgrid/includecc.cc>
using namespace Dune;
typedef ALU3dGrid<DIM,DIM_OF_WORLD,tetra> GridType;
#endif

#include <dune/fem/discfuncarray.hh>
#include <dune/fem/dfadapt.hh>
#include <dune/fem/lagrangebase.hh>

typedef FunctionSpace < double , double, DIM , 1 > FuncSpace;
typedef DefaultGridIndexSet < GridType , GlobalIndex > IndexSetType;
typedef LagrangeDiscreteFunctionSpace < FuncSpace , GridType, IndexSetType , 0 > FuncSpaceType ;
typedef DiscFuncArray < FuncSpaceType > DiscFuncType;
typedef DFAdapt < FuncSpaceType > DFType;

#include <dune/io/visual/grapedatadisplay.hh>
#include <dune/common/stdstreams.cc>

int main (int argc, char **argv)
{
  // this leads to the same number of points for SGrid and AlbertGrid
#if SGRID
  int n[DIM];
  double h[DIM];
  for(int i=0; i<DIM; i++)
  {
    n[i] = 3; h[i] = 1.0;
  }
  GridType grid ((int *) &n, (double *) &h );
  int level = 0;
  if(argc == 2) level = atoi(argv[1]);
  grid.globalRefine (level);
#endif
 
#if UGRID
  GridType grid(200,30);
  if(argc < 2)
  {
    fprintf(stderr,"usage: %s <AmiraMesh> \n",argv[0]);
    abort();
  }
  int level = 0;
  if(argc == 3)
    level = atoi(argv[2]);

  AmiraMeshReader< GridType > am;
  am.read ( grid, argv[1] );
  //grid.globalRefine (level);
  AmiraMeshWriter< GridType > amw;
  std::string fn (  "testgrid.am" );
  amw.writeGrid ( grid ,fn );
#endif

#if YGRID
  // initialize MPI
  MPI_Init(&argc,&argv);
       
  std::cout << "We are using YaspGrid! \n";

  Vec<double,DIM> lang = 1.0;
  Vec<int,DIM>  anz = 4;
  Vec<bool,DIM> per; 
  for(int i=0; i<DIM; i++) per(i) = true;
  
  GridType grid (MPI_COMM_WORLD, lang, anz, per ,0 );
  for(int i=0; i<6; i++)  grid.globalRefine (1);
#endif

#if AGRID
  if(argc < 2)
  {
    fprintf(stderr,"usage: %s <Albert Macro Triang> \n",argv[0]);
    abort();
  }
  int level = 0;
  if(argc == 3)
    level = atoi(argv[2]);

  GridType grid( argv[1] );
  grid.globalRefine (level);
  
  //AmiraMeshWriter< GridType > am;
  //std::string fn (  "bucklev" );
  //am.writeGrid ( grid ,fn );
#endif

#if BGRID
  if(argc < 2)
  {
    fprintf(stderr,"usage: %s <BSGrid Macro Triang> \n",argv[0]);
    abort();
  }
  int level = 0;
  if(argc == 3)
    level = atoi(argv[2]);
  
 
#ifdef _ALU3DGRID_PARALLEL_
  GridType grid( argv[1], MPI_COMM_WORLD );
#else
  GridType grid( argv[1] );
#endif
  grid.globalRefine (level);
 
  GrapeDataIO< GridType > dataIO;
  dataIO.writeGrid(grid,ascii,"grid",0.0,0);
  //AmiraMeshWriter< GridType > am;
  //std::string fn (  "forward3d.am" );
  //am.writeGrid ( grid ,fn );
#endif


#if 0
#if !BRGID
  for(int i=0; i<=grid.maxlevel(); i++)
  {
    std::cout << "Print size of level " << i << "\n";
    std::cout << grid.size(i,0) << " number of Elements! \n";
    std::cout << grid.size(i,1) << " number of Faces! \n";
    std::cout << grid.size(i,2) << " number of Edges! \n";
    std::cout << grid.size(i,DIM) << " number of Points! \n";
  }
#endif
#endif

  // can be GrapeDisplay or GrapDataDisplay
  GrapeGridDisplay < GridType > grape(grid,-1);
  grape.display();

  /*
  {
    GrapeDataDisplay < GridType , DFType > grape(grid,-1);
    typedef DofManager< GridType > DofManagerType;
    DofManagerType dm ( grid );
 
    IndexSetType iset ( grid );
    FuncSpaceType  linFuncSpace ( grid , iset , dm , grid.maxlevel() ); 
    DiscFuncType df ( linFuncSpace );
    DFType df2 ( linFuncSpace );
    df.set(1.0);
    df2.set(0.5);
    grape.dataDisplay(df2);
  }
  */

#if YGRID 
  MPI_Finalize();
#endif

  return 0;
  
}

