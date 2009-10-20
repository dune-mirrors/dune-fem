//**************************************************************
//  (C) written and directecd by Robert Kloefkorn 
//**************************************************************
#ifndef DATADISP_READTUPLE_CC
#define DATADISP_READTUPLE_CC

//- system includes 
#include <stack> 

//- include grape io stuff 
#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/file/grapedataio.hh>
#include <dune/fem/io/file/iotuple.hh>
#include <dune/fem/io/file/iointerface.hh>

#include <dune/fem/io/visual/grape/datadisp/grcommon.hh>

//! type of used grid 
typedef GridType GR_GridType;

//! type of discrete function tuple 
typedef GR_InputType GR_DiscFuncType;

static std::stack<GR_GridType *> gridStack;

static GrapeDataIO <GR_GridType> dataIO;

template <class T> 
inline void deleteObjects(std::stack<T *> & stack);

inline void dataDispErrorExit(std::string msg) 
{
  std::cerr << msg << std::endl; 
  std::cerr.flush();
  exit(EXIT_FAILURE);
}

#if !HAVE_GRAPE
inline void removeWarnings() 
{
  // this is to remove compiler warnings 
  maxlevelButton = 0;
  minMaxColorbar = 0;
  iteratorButton = 0;
  defaultIteratorValue = 0;
  partitionTypeButton = 0;
  globalTsc = 0;
}
#endif

inline void readTupleData(const char * path, const char * filename, 
            double & time , 
            const int n, 
            const int timestep, 
            const int myRank,
            const int numProcs, 
            INFO* info)
{
  // check whether we use a fixed mesh 
  const bool fixedMesh = (info->fix_mesh == 1) ? true : false;

  // the grid is new if not a fixed mesh or the stack is empty 
  const bool newGrid = ( ! fixedMesh || gridStack.empty() );

  // init grid 
  GR_GridType * grid = ( ! newGrid ) ? gridStack.top() : 0;

  assert(filename);
  std::string fn (filename);

  IOTuple<GR_DiscFuncType>::ReturnType* tup = 
    IOTuple<GR_DiscFuncType>::input(dataIO,grid,time,timestep,path,fn);
  std::cout << "Finished reading grid" << std::endl;

  // push all new grids to grid stack 
  if( newGrid ) gridStack.push(grid);
  
  // do some processing 
  process(*grid,*tup,time,timestep,myRank,numProcs);

  IOTuple<GR_DiscFuncType>::removeData(*tup);
  delete tup;
  
  if( ! fixedMesh ) delete grid;
}

// setup the hole data tree for grape 
inline INFO * readData(INFO * info , const char * path, int i_start, int i_end, 
    int i_delta, int n, double timestep, int numProcs) 
{
  double t_start = 1e308;
  double t_end = -t_start, t_act = 0.0;
  int  ntime, n_step = 0;

  for (ntime = i_start; ntime <= i_end; ntime += i_delta)
  {
    printf("timestep = %d | last timestep = %d | stepsize = %d\n", ntime, i_end, i_delta);
    {
      int anzProcs = numProcs;
      
      for(int proc=0; proc<anzProcs; ++proc)
      {
        assert(path || numProcs <= 1); 

        int anz = n; // (n > 0) ? n : 1;
        for(int i=0; i<anz; ++i)
        {
          // use standard procedure to create path name 
          std::string newpath = 
            IOInterface::createRecoverPath(path,proc,info[i].name,ntime);

          readTupleData(newpath.c_str(), info[i].name, 
                        t_act , i , ntime, proc, numProcs, info);
        }
      }
    }
    printf("actual time: %f (timestep size: %e)\n\n",t_act,timestep);


    if (ntime == i_start) t_start = t_end = t_act;
    t_start = std::min(t_start, t_act);
    t_end = std::max(t_end, t_act);

    if (timestep > 0) t_act += timestep*i_delta;
    n_step++;
  }
#if !HAVE_GRAPE
  removeWarnings();
#endif
  return info;
}

template <class T> 
inline void deleteObjects(std::stack<T *> & stack) 
{
  while(! stack.empty() )
  {
    T * obj = stack.top();
    stack.pop();
    delete obj;
  }
  return;
}
 
inline void deleteAllObjects() 
{
  deleteObjects(gridStack);
  return ;
}
#endif
