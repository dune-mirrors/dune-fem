#include <iostream>
#include <vector>
#include <assert.h>
#include <string.h>

#include <config.h>

#include <dune/common/stdstreams.cc>

//#include <dune/grid/sgrid.hh>

#define LARGE 1.0E308
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif

#if HAVE_ALUGRID
#include "dune/grid/alu3dgrid.hh"
//#include "dune/grid/alu3dgrid/includecc.cc"
#endif

#include <dune/fem/dfadapt.hh>
#include <dune/fem/lagrangebase.hh>
#include <dune/fem/discfuncarray.hh>
#include <dune/common/stack.hh>

#include <dune/grid/common/leafindexset.hh>

#include <dune/io/file/grapedataio.hh>

using namespace Dune;

#ifndef DIM 
#define DIM 3 
#endif

#ifndef DIM_OF_WORLD
#define DIM_OF_WORLD 3 
#endif

#include <dune/io/visual/grapedatadisplay.hh>
#include "globaldefs.hh"

static const int char_space = 2;

void dunedispErrorExit(const char * msg);

static GR_DiscFuncSpaceType * globalSpace = 0;
static GR_IndexSetType * indexSet = 0;

static Stack<GR_GridType *> gridStack;
static Stack<GrapeDispType *> dispStack;
static Stack<GR_DiscFuncSpaceType *> fsStack;
static Stack<GR_DiscFuncType *> funcStack;
static Stack<GR_IndexSetType *> indexStack;

static GrapeDataIO < GR_GridType> dataIO;

template <class T> 
void deleteObjects(Stack<T *> & stack);

// read data from file and generate discrete function 
void readFuncData ( GrapeDispType& disp, GR_DiscFuncSpaceType &fspace, 
     const char * path, const char * filename , double time , int ntime , int proc )
{
  char * fn = 0;

  if(path) 
  {
    fn = new char [strlen(path) + strlen(filename) + char_space]; assert(fn);
    sprintf(fn,"%s/%s",path,filename);
  }
  else 
  {
    fn = new char [strlen(filename) + char_space]; assert(fn);
    sprintf(fn,"%s",filename);
  }

  GR_DiscFuncType *df = new GR_DiscFuncType ( filename, fspace );
  funcStack.push(df);
  
  dataIO.readData (*df, fn , ntime ) ;
  disp.addData(*df,filename,time);

  return ;
}

// read Grid and create GrapeDataDisplay with Hmesh 
GrapeDispType * readGrid(const char * path, const char * filename, 
                         double & time , int ntime, int myRank )
{
  GR_GridType * grid = new GR_GridType ();
  gridStack.push(grid);
  
  assert(filename);
  char * fn = 0;
  if(path)
  {
    fn = new char[ strlen(path) + strlen(filename) + char_space];
    assert(fn);
    sprintf(fn,"%s/%s",path,filename);
  }
  else 
  {
    fn = new char[ strlen(filename) + char_space];
    assert(fn);
    sprintf(fn,"%s",filename);
  }
  
  std::cout << "Make new Grapedisplay for grid = " << fn << "\n";
  dataIO.readGrid( *grid, fn , time , ntime );
  
  GrapeDispType * disp = new GrapeDispType ( *grid, myRank );  
  dispStack.push(disp);
  if(fn) delete [] fn;
  return disp;
}

// read dof manager from file 
void readDofManager(GR_DofManagerType & dm, const char * path, int ntime) 
{
  // generate dof manager name 
  int length = 0;
  if(path) length = strlen(path);
  char * fn = new char [length + 10];
  assert(fn);
  if(path) sprintf(fn,"%s/dm",path);
  else sprintf(fn,"dm");
  dm.read(fn,ntime);
  if(fn) delete [] fn;
}

// read all data that belong to grid with name info[n].name 
INFO *makeData( GrapeDispType * disp, INFO * info , const char * path, 
    const char * filename, double & time ,int n, int ntime, bool fix_mesh,
    int proc=0)
{
  if(info[n].datinf)
  {
    GR_DiscFuncSpaceType * space = 0;
    if(!info[n].fix_mesh)
    {
      GR_DofManagerType * dm = & GR_DofManagerFactoryType::getDofManager (disp->getGrid());
       
      GR_IndexSetType * iSet = new GR_IndexSetType ( disp->getGrid() );
      indexStack.push(iSet);
      space  = new GR_DiscFuncSpaceType ( disp->getGrid() , *iSet, *dm, disp->getGrid().maxlevel() );
      readDofManager(*dm,path,ntime); 
      
      fsStack.push(space);
    }
    else 
    {
      std::cout << "We use the same space because uniform Grids! \n";
      if(!globalSpace) 
      {
        GR_DofManagerType * dm = & GR_DofManagerFactoryType::getDofManager (disp->getGrid());

        indexSet = new GR_IndexSetType ( disp->getGrid() );
        globalSpace = new GR_DiscFuncSpaceType ( disp->getGrid() , *indexSet, *dm, disp->getGrid().maxlevel());
        readDofManager(*dm,path,ntime); 
      }
      space = globalSpace;
    }
    assert(space != 0);

    DATAINFO * dinf = info[n].datinf;
    while (dinf) 
    {
      readFuncData ( *disp , *space , path, dinf->name ,time, ntime,proc);
      dinf = dinf->next;
    }
  }    
  else 
  {
    std::cout << "***  No Data, displaying only grid! ***\n";
  }
  // store mesh and data in timescene tree 
  return info;
}

// setup the hole data tree for grape 
INFO * readData(INFO * info , const char * path, int i_start, int i_end, 
    int i_delta, int n, double timestep, int numProcs) 
{
  double f_t_start = LARGE;
  double t_start = LARGE;
  double t_end = -LARGE, t_act=0.0;
  GrapeDispType *disp = 0;
  
  int  ntime, n_step = 0;
  
  if(info[n].fix_mesh)
  {
    disp = readGrid( path, info[0].name, t_act , 0, 0);
    f_t_start = t_act;
  }

  for (ntime = i_start; ntime <= i_end; ntime += i_delta)
  {
    printf("timestep = %d | last timestep = %d | stepsize = %d\n", ntime, i_end, i_delta);
    if (!info[n].fix_mesh)
    {
      int anzProcs = numProcs;
      if(numProcs > 1) anzProcs--;
      
      for(int proc=0; proc<anzProcs; proc++)
      {
        assert(path || numProcs <= 1); 
        char * newpath = new char[strlen(path)+5];
        sprintf(newpath,"%s_%d",path,proc); 
        if(numProcs <= 1)
        {
          sprintf(newpath,"%s",path); 
        }
      
        GrapeDispType *newdisp = 0;
        int anz = (n > 0) ? n : 1;
        for(int i=0; i<anz; i++)
        {
          newdisp = readGrid( newpath, info[i].name, t_act , ntime, proc );
          assert(newdisp != 0);
          info = makeData(newdisp,info,newpath,info[i].name,t_act, n ,ntime,proc);
        }
        newdisp->addMyMeshToTimeScene(info[0].tsc,t_act,proc);
        if(newpath) delete [] newpath;
      }
    }
    else
    {
      info = makeData(disp,info,path,info[0].name,t_act,n,ntime,numProcs);
      disp->addMyMeshToTimeScene(info[0].tsc,t_act, -1 );//proc);
      t_act = f_t_start+ntime*timestep;
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

#include "printhelp.cc"

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
        dunedispErrorExit("usage: -i `increment'\n");
      i_delta = atoi(argv[i+1]);
      i += 2;
    }
    else if (!strcmp(argv[i], "-s"))
    {
      if (i+1 == argc)
        dunedispErrorExit("usage: -s `dataprefix'\n");
      
      DATAINFO * dinf = (DATAINFO *) std::malloc(sizeof(DATAINFO));
      assert(dinf);
      dinf->name = argv[i+1];
      dinf->vector = 0;
      /* seems wrong order, but grape truns it arround, we can do nothing else here */
      dinf->next = info[n].datinf; 
      info[n].datinf = dinf;

      i += 2;
    }
    else if (!strcmp(argv[i], "-t"))
    {
      if (i+1 == argc)
        dunedispErrorExit("usage: -t `time step size'\n");
      timestep = atof(argv[i+1]);
      i += 2;
    }
    else if (!strcmp(argv[i], "-m"))
    {
      if (i+1 == argc)
        dunedispErrorExit("usage: -m `gridprefix'\n");
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
        dunedispErrorExit("usage: -pg `number of procs'\n");
      parallel += atoi(argv[i+1]);
      i += 2;
    }
    else if (!strcmp(argv[i], "-p"))
    {
      if (i+1 == argc)
        dunedispErrorExit("usage: -p `path'\n");
      path = argv[i+1];
      i += 2;
    }
    else
    {
      fprintf(stderr,"unknow option %s\n", argv[i]);
      exit(1);
    }
    printf("i = %d, argc = %d\n", i, argc);
  }

  timeSceneInit(info, n , parallel , time_bar);
  readData(info, path,i_start,i_end,i_delta,n,timestep,parallel);
  
  // run grape 
  displayTimeScene(info,parallel);
 
  deleteObjects(funcStack);
  if(globalSpace) delete globalSpace;
  deleteObjects(fsStack);
  deleteObjects(indexStack);

  deleteObjects(dispStack);
  deleteObjects(gridStack);

  return (EXIT_SUCCESS);
}

template <class T> 
void deleteObjects(Stack<T *> & stack) 
{
  while(! stack.empty() )
  {
    T * obj = stack.pop();
    delete obj;
  }
  return;
}
 
void dunedispErrorExit(const char * msg) 
{
  std::cerr << msg << std::endl; 
  std::cerr.flush();
  exit(EXIT_FAILURE);
}

