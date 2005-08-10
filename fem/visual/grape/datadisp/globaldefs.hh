#ifndef __DUNE_GLOBALDEFS_HH__
#define __DUNE_GLOBALDEFS_HH__

// here we have to define what grid and functionspace, functions we use 

#ifndef DIM 
#define DIM 3
#endif
#ifndef DIM_OF_WORLD
#define DIM_OF_WORLD 3
#endif

typedef double REAL;

//**************************************************************************

//#define DEBUG_LEAFINDEXSET

#if HAVE_ALBERTA
#define GRIDNAME AlbertaGrid
typedef GRIDNAME < DIM , DIM_OF_WORLD > GR_GridType; 
#endif

#if HAVE_ALUGRID 
typedef ALU3dGrid < DIM , DIM_OF_WORLD , tetra > GR_GridType; 
#endif

typedef FunctionSpace <double ,double , DIM , DIM+2 >  GR_FunctionSpaceType;
//typedef DofManager<GR_GridType , AdaptiveLeafIndexSet<GR_GridType> > GR_DofManagerType;
typedef DofManager<GR_GridType,DataCollectorInterface<GR_GridType,GR_GridType::ObjectStreamType> > GR_DofManagerType;
typedef DofManagerFactory <GR_GridType,DataCollectorInterface<GR_GridType,GR_GridType::ObjectStreamType> > GR_DofManagerFactoryType;

//typedef GR_GridType :: LeafIndexSetType GR_IndexSetType;
typedef DefaultGridIndexSet<GR_GridType, GlobalIndex > GR_IndexSetType;
//typedef AdaptiveLeafIndexSet<GR_GridType> GR_IndexSetType;
//typedef DofManager<GR_GridType , DefaultGridIndexSet<GR_GridType, LevelIndex > > GR_DofManagerType;

typedef LagrangeDiscreteFunctionSpace < GR_FunctionSpaceType , GR_GridType , GR_IndexSetType, 0, GR_DofManagerType> 
        GR_DiscFuncSpaceType;
        
//typedef DGDiscreteFunctionSpace < GR_FunctionSpaceType , GR_GridType , 1> 
//        GR_DiscFuncSpaceType;

//typedef RaviartThomasSpace < GR_FunctionSpaceType , GR_GridType , 0> 
//        GR_DiscFuncSpaceType;

typedef DFAdapt < GR_DiscFuncSpaceType > GR_DiscFuncType;

typedef GrapeDataDisplay<GR_GridType , GR_DiscFuncType > GrapeDispType;

#endif
