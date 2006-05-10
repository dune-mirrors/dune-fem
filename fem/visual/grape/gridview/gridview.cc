//***************************************************************
//
//  (C) written and directed by Robert Kloefkorn 
//
//***************************************************************
#include <iostream>
#include <config.h>

// set exactly one of the following to 1, the others 0
#define SGRID 1
#define UGRID 0
#define AGRID 0
#define BGRID 0
#define YGRID 0

#if SGRID
#include <dune/grid/sgrid.hh>
using namespace Dune;

#undef dim 
#undef dim_OF_WORLD

static const int dim = 3;
static const int dimworld = 3;

typedef SGrid <dim,dimworld> GridType;
#endif

#if UGRID

#include <dune/grid/uggrid.hh>
#include <dune/io/file/amirameshreader.hh>
#include <dune/io/file/amirameshwriter.hh>

#undef dim 
#undef dimworld

#ifdef _2 
static const int dim = 2;
static const int dimworld = 2;
#else 
static const int dim = 3;
static const int dimworld = 3;
#endif
using namespace Dune;
typedef UGGrid< dim,dimworld> GridType;
#endif

#if YGRID
#include <dune/grid/yaspgrid.hh>
static const int dim = 3;
static const int dimworld = 3;
using namespace Dune;
typedef YaspGrid <dim,dimworld> GridType;
#endif

#if AGRID
#include <dune/grid/albertagrid.hh>
using namespace Dune;
static const int dim = DUNE_PROBLEM_DIM;
static const int dimworld = DUNE_WORLD_DIM;
typedef AlbertaGrid<dim,dimworld> GridType;
#endif

#if BGRID

static const int dim = 3;
static const int dimworld = 3;

#define PROBLEMdim 3
#define WORLDdim 3

// the following file does not exist... renamed or what?
//#include <dune/grid/alu3dgrid/includecc.cc>
//#include <dune/grid/alu3dgrid.hh>
#include <dune/grid/alugrid.hh>
using namespace Dune;
//typedef ALU3dGrid<dim,dimworld,tetra> GridType;
typedef ALU3dGrid<dim,dimworld,hexa> GridType;
// the new syntax should be the following, but incurs error
// typedef ALUCubeGrid<dim,dimworld> GridType; 

#endif

#include <dune/grid/io/visual/grapegriddisplay.hh>
#include <dune/common/stdstreams.cc>

int main (int argc, char **argv)
{
  // this leads to the same number of points for SGrid and AlbertGrid
#if SGRID
  int n[dim];
  double h[dim];
  for(int i=0; i<dim; i++)
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

//  Vec<double,dim> lang = 1.0;
//  Vec<int,dim>  anz = 4;
//  Vec<bool,dim> per; 
//  for(int i=0; i<dim; i++) per(i) = true;
  FieldVector<double,dim> lang(1.0);
  FieldVector<int,dim>  anz(4);
  FieldVector<bool,dim> per(true); 
  
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
 
  //GrapeDataIO< GridType > dataIO;
  //dataIO.writeGrid(grid,xdr,"grid",0.0,0);
  //AmiraMeshWriter< GridType > am;
  //std::string fn (  "forward3d.am" );
  //am.writeGrid ( grid ,fn );
#endif

  {
    GrapeGridDisplay < GridType > grape(grid,-1);
    grape.display();

//    GrapeDataDisplay < GridType > grape(grid,-1);
//    typedef DofManager< GridType > DofManagerType;
//    typedef DofManagerFactory < DofManagerType> DofManagerFactoryType; 
//    IndexSetType iset ( grid );    
//    LeafGridPart < GridType > gpart ( grid );
//    FuncSpaceType  space( gpart ); 
//    DFType df ( space );
//    setFunc(grid,df);
//    grape.dataDisplay(df,true);

  }

#if YGRID 
  MPI_Finalize();
#endif

  return 0;
  
}

