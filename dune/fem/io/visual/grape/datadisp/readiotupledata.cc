//**************************************************************
//  (C) written and directecd by Robert Kloefkorn
//**************************************************************
#ifndef DATADISP_READIOTUPLE_CC
#define DATADISP_READIOTUPLE_CC

// needed in readiotparams.cc
#define USE_GRAPE_DISPLAY

//- system includes
#include <stack>

//- include grape io stuff
#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/file/iotuple.hh>
#include <dune/fem/io/file/iointerface.hh>

//! type of used grid
typedef Dune::GridSelector::GridType GR_GridType;

//! type of GrapeDisplay
typedef GrapeDataDisplay<GR_GridType> GrapeDispType;

//! type of discrete function tuple
typedef GR_InputType GR_DiscFuncType;

static std::stack<GR_GridType *> gridStack;
static std::stack<GrapeDispType *> dispStack;

template <class T>
inline void deleteObjects(std::stack<T *> & stack);

inline void dataDispErrorExit(std::string msg)
{
  std::cerr << msg << std::endl;
  std::cerr.flush();
  exit(EXIT_FAILURE);
}

inline GrapeDispType * readTupleData(const char * path, const char * filename,
            double & time , int n,
            int timestep, int myRank, int mySize,
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
  DATAINFO * dinf = info[n].datinf;

  Fem::IOTuple<GR_DiscFuncType>::ReturnType* tup =
    Fem::IOTuple<GR_DiscFuncType>::input(grid,time,myRank,mySize,path,fn);
  std::cout << "Finished reading grid" << std::endl;

  // push all new grids to grid stack
  if( newGrid ) gridStack.push(grid);

  GrapeDispType * disp = new GrapeDispType ( *grid, myRank );
  dispStack.push(disp);

  // discrete functions of non-valid data are removed
  Fem::IOTuple<GR_DiscFuncType>::addToDisplayOrRemove(*disp,dinf,time,*tup);

  //Element<0> :: get( *tup );
  //    *(Element<0> :: get( *tup ))
  // do some post processing
  postProcessing(*disp,*grid,time,timestep, *tup);
  return disp;
}

// setup the hole data tree for grape
inline INFO * readData(INFO * info , const char * path, int i_start, int i_end,
    int i_delta, int n, double timestep, int numProcs)
{
  double t_start = 1e308;
  double t_end = -t_start, t_act = 0.0;
  typedef CombinedGrapeDisplay < GrapeDispType > CombinedDisplayType;

  CombinedGrapeDisplay < GrapeDispType > * comdisp = new CombinedDisplayType ();

  int  ntime, n_step = 0;

  for (ntime = i_start; ntime <= i_end; ntime += i_delta)
  {
    printf("timestep = %d | last timestep = %d | stepsize = %d\n", ntime, i_end, i_delta);
    {
      int anzProcs = numProcs;

      for(int proc=0; proc<anzProcs; ++proc)
      {
        assert(path || numProcs <= 1);
        GrapeDispType *newdisp = 0;

        int anz = n; // (n > 0) ? n : 1;
        for(int i=0; i<anz; ++i)
        {
          // use standard procedure to create path name
          std::string newpath =
            Fem::IOInterface::createRecoverPath(path,proc,info[i].name,ntime);

          newdisp = readTupleData(newpath.c_str(), info[i].name,
                                  t_act , i , ntime, proc, anzProcs, info);

          if( comdisp )
          {
            assert(newdisp != 0);
            assert( comdisp );
            comdisp->addDisplay( *newdisp );
          }
        }
        assert(newdisp);
        newdisp->addMyMeshToTimeScene(info[0].tsc,t_act,proc);
      }

      if( comdisp )
      {
        // this is the combine object, which is put to the last time scene
        comdisp->addMyMeshToGlobalTimeScene(t_act,0);
      }
    }
    printf("actual time: %f (timestep size: %e)\n\n",t_act,timestep);


    if (ntime == i_start) t_start = t_end = t_act;
    t_start = std::min(t_start, t_act);
    t_end = std::max(t_end, t_act);

    if (timestep > 0) t_act += timestep*i_delta;
    n_step++;
  }
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
  deleteObjects(dispStack);
  deleteObjects(gridStack);
  return ;
}
#endif
