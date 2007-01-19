/**************************************************************************
**       Title: functionview: program visualizing alu3d grid with simple
**              scalar function on it. Slight extension of roberts gridview
**    $RCSfile$
**   $Revision$$Name$
**       $Date$
**   Copyright: GPL $Author$
** Description:
**
**    
**
**-------------------------------------------------------------------------
**
**  $Log$
**  Revision 1.2  2006/06/19 14:19:05  haasdonk
**  adopted to new dune-structure
**
**  Revision 1.1  2006/05/10 14:32:19  haasdonk
**  modified for latest dune-modularization
**
**
**************************************************************************/


#include <iostream>
#include <config.h>
#include <math.h>

// the following is required if we do not link to alberta
//#undef HAVE_ALBERTA 
//#define HAVE_ALBERTA 0

//#pragma interface
//#pragma implementation

#define SGRID 0
#define UGRID 0
#define AGRID 1
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
using namespace Dune;
static const int dim = 3;
static const int dimworld = 3;
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

#include <dune/grid/alugrid.hh>
using namespace Dune;
//typedef ALU3dGrid<dim,dimworld,tetra> GridType;
//typedef ALU3dGrid<dim,dimworld,hexa> GridType;
typedef ALUCubeGrid<dim,dimworld> GridType;
#endif

//#include <dune/quadrature/barycenter.hh>
//#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/quadrature/elementquadrature.hh>
//#include <dune/fem/discfuncarray.hh>
#include <dune/fem/discretefunction/dfadapt.hh>
#include <dune/fem/space/lagrangespace/lagrange.hh>
#include <dune/grid/common/gridpart.hh>

#if AGRID 
//typedef FunctionSpace < double , double, dim , dim > FuncSpace;
typedef FunctionSpace < double , double, dim , 1 > FuncSpace;
#else 
typedef FunctionSpace < double , double, dim , 1 > FuncSpace;
#endif


//typedef DefaultGridIndexSet < GridType , GlobalIndex > IndexSetType;
//#if SGRID
//  typedef LevelGridPart < GridType > GridPartType; 
//#else // SGRID has no hierarchic iterator
   typedef LeafGridPart < GridType > GridPartType; 
//#endif

// higher degrees than 1 currently not supported!!
typedef LagrangeDiscreteFunctionSpace < FuncSpace , GridPartType, 1 > FuncSpaceType ; 
typedef DFAdapt < FuncSpaceType > DiscFuncType;
typedef DiscFuncType DFType;

//typedef DFAdapt < FuncSpaceType > DFType;
typedef FieldVector<double, dim> FVType;

#include <dune/grid/io/visual/grapedatadisplay.hh>
#include <dune/fem/io/file/grapedataio.hh>
#include <dune/common/stdstreams.cc>

//#include "/usr/people/haasdonk/fuelcell/probs/FuelCell/3_FuelCellDiscretization/fc_common_discr.hh"

void global_function(const double* global_pt, double *res )
{
//   if (global_pt[0] < g_bip_width)
//       *res = 0.0;
//   else if (global_pt[0] < g_bip_width + g_gdl_width) 
//       *res = 1.0;
//   else if (global_pt[0] < g_bip_width + g_gdl_width +g_cat_width) 
//       *res = 2.0;
//   else if (global_pt[0] < g_bip_width + g_gdl_width +g_cat_width 
//            + g_mem_width) 
//       *res = 3.0;
//   else if (global_pt[0] < g_bip_width + g_gdl_width +2 * g_cat_width 
//            + g_mem_width) 
//       *res = 2.0;
//   else if (global_pt[0] < g_bip_width + 2 * g_gdl_width +2 * g_cat_width 
//            + g_mem_width) 
//       *res = 1.0;
//   else 
//       *res = 0.0;  

  if (dim>=3)
      * res = (1-cos(10*global_pt[0]-1.1)) * (1-cos(10*global_pt[1]+0.25))* 
          (1-cos(5*global_pt[2]+1.5));
  else
      * res = (1-cos(10*global_pt[0]-1.1)) * (1-cos(10*global_pt[1]+0.25));
  
}


//! velocity of the transport problem
// template <class FieldVectorType> 
// void rotatingPulse(const FieldVectorType &x, FieldVectorType & velo)
// {
//   assert(velo.dim() >= 2);
//   // rotating pulse 
//   velo[0] = -4.0*x[1];
//   velo[1] =  4.0*x[0];
//   if(velo.dim() > 2) velo[2] = x[2];
//   return;
// }


template <class GridType, class DiscFuncType> 
void setFunc (GridType & grid, DiscFuncType & df )
{
  typedef typename GridType :: template Codim<0> :: LeafIterator LeafIterator;
  typedef typename DiscFuncType :: LocalFunctionType LFType;
  
//  typedef BaryCenterQuad < 
//      typename DiscFuncType :: RangeFieldType ,
//      typename DiscFuncType::DomainType, 0 > QuadratureType;

//  typedef Quadrature <typename DiscFuncType::DomainType, dim > 
//      QuadratureType;

  typedef ElementQuadrature <GridType, 0> 
      QuadratureType;

  LeafIterator endit = grid.template leafend<0>();
  LeafIterator it = grid.template leafbegin<0> (); 

  std::cout << "dimrange:" << DiscFuncType :: DiscreteFunctionSpaceType :: DimRange << "\n";
//  std::cout << "dimrange:" << DiscFuncType :: DimRange << "\n";
  
  if( it == endit ) 
  {
    std::cerr << "ERROR, empty grid! \n";
    abort();
  }
  
  // 0th-order quadrature (is hopefully a barycenter-evaluation)
//  GeometryType geo = it->geometry().type();
  QuadratureType quad( *it, 0 );
  
  std::cout << "total nops of quadrature " << quad.nop() << "\n";
  std::cout << "point 0 of quadrature" << quad.point(0) << "\n";
    
  for( ; it != endit; ++it)
  {
    FieldVector<double,1> domnr;
//    FieldVector<double,dim> velo;
    
//    std::cout <<  "dim = " << dim << "\n" ;
    
//    FieldVector<double,3> bary = it->geometry().global( quad.point(0)); 
    FVType bary = it->geometry().global( quad.point(0)); 
    std::cout << "trying to access local function ..." << std::flush;
    
    LFType lf = df.localFunction(*it); 
    std::cout << " ... success\n";

//    rotatingPulse(bary,velo);
    global_function(&bary[0], &domnr[0]);
//    for(int i=0; i<lf.numDofs(); i++) 
//      lf[i] = velo[i];
    lf[0] = domnr[0];
  }
}


int main (int argc, char **argv)
{
  // this leads to the same number of points for SGrid and AlbertGrid
#if SGRID
  int n[dim];
  double h[dim];
  for(int i=0; i<dim; i++)
  {
    n[i] = 10; h[i] = 1.0;
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

  Vec<double,dim> lang = 1.0;
  Vec<int,dim>  anz = 4;
  Vec<bool,dim> per; 
  for(int i=0; i<dim; i++) per(i) = true;
  
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
 
//  GrapeDataIO< GridType > dataIO;
//  dataIO.writeGrid(grid,xdr,"grid",0.0,0);
  //AmiraMeshWriter< GridType > am;
  //std::string fn (  "forward3d.am" );
  //am.writeGrid ( grid ,fn );
#endif

  {
    GrapeDataDisplay < GridType > grape(grid,-1);
//    typedef DofManager< GridType > DofManagerType;
//    typedef DofManagerFactory < DofManagerType> DofManagerFactoryType; 

//    IndexSetType iset ( grid );
    
//#if SGRID
//    LevelGridPart < GridType > gpart ( grid, level );
//#else
    LeafGridPart < GridType > gpart ( grid );
//#endif

    FuncSpaceType  space( gpart );     
    DFType df ( space );
    setFunc(grid,df);
    grape.dataDisplay(df,false);  // no vector data!
  }

#if YGRID 
  MPI_Finalize();
#endif

  return 0;
}

