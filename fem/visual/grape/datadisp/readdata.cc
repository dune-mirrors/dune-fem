//**************************************************************
//  (C) written and directecd by Robert Kloefkorn 
//**************************************************************
#ifndef __DATADISP_READDATA_CC__
#define __DATADISP_READDATA_CC__

#define LARGE 1.0E308

#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/file/grapedataio.hh>


void dataDispErrorExit(std::string msg);

static GR_DiscFuncSpaceType * globalSpace = 0;

static Stack<GR_GridType *> gridStack;
static Stack<GrapeDispType *> dispStack;
static Stack<GR_DiscFuncSpaceType *> fsStack;
static std::list <GR_DiscFuncType *> funcStack;
//static Stack<GR_IndexSetType *> indexStack;
static Stack<GR_GridPartType *> gridPartStack;

static GrapeDataIO < GR_GridType> dataIO;

template <class T> 
void deleteObjects(Stack<T *> & stack);

// read data from file and generate discrete function 
void readFuncData ( GrapeDispType& disp, GR_DiscFuncSpaceType &fspace, 
     const char * path, const DATAINFO * dinf, double time , int ntime )
{
  const char * name = dinf->base_name;
  assert(name);

  std::string fakename ( name );

  typedef std::list <GR_DiscFuncType *> :: iterator ListIterator; 
  for(ListIterator it = funcStack.begin(); it != funcStack.end(); ++it)
  {
    std::string tmpname ( (*it)->name() );
    if(fakename == tmpname) 
    {
      //GR_DiscFuncType * df = (*it);
      //disp.addData(*df,dinf,time); 
      break;//return ;
    }
  }

  std::string fn (path);
  if(path) fn += "/"; 
  fn += name;

  GR_DiscFuncType *df = new GR_DiscFuncType (name , fspace );
  funcStack.insert(funcStack.begin(),df);
  dataIO.readData (*df, fn , ntime ) ;
  disp.addData(*df,dinf,time);

  return ;
}

// read Grid and create GrapeDataDisplay with Hmesh 
GrapeDispType * readGrid(const char * path, const char * filename, 
                         double & time , int ntime, int myRank )
{
#if SGRID
  GR_GridType * grid = new GR_GridType ();
  grid->globalRefine(9);
#else

  GR_GridType * grid = 0 ;
  /*
  GR_GridType * grid = new GR_GridType (
#ifdef _ALU3DGRID_PARALLEL_
  MPI_COMM_WORLD 
#endif
      );
      */

#endif

  gridStack.push(grid);

  assert(filename);
  std::string fn (path);
  if(path) fn += "/"; 
  fn += filename;
  
  std::cout << "Make new Grapedisplay for grid = " << fn << "\n";
#if SGRID

#else
  dataIO.readGrid( *grid, fn , time , ntime );
#endif
  
  GrapeDispType * disp = new GrapeDispType ( *grid, myRank );  
  dispStack.push(disp);
  
  return disp;
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
//      GR_DofManagerType * dm = & GR_DofManagerFactoryType::getDofManager (disp->getGrid());
       
      //GR_IndexSetType * iSet = new GR_IndexSetType ( disp->getGrid() );
      //indexStack.push(iSet);
      
      GR_GridPartType* gridPart = new GR_GridPartType(disp->getGrid());//,*iSet);
      gridPartStack.push(gridPart);
      space  = new GR_DiscFuncSpaceType (*gridPart);
      
      fsStack.push(space);
    }
    else 
    {
      std::cout << "We use the same space because uniform Grids! \n";
      if(!globalSpace) 
      {
//        GR_DofManagerType * dm = & GR_DofManagerFactoryType::getDofManager (disp->getGrid());

        //GR_IndexSetType * indexSet = new GR_IndexSetType ( disp->getGrid() );
        //indexStack.push(indexSet);
         
        GR_GridPartType* gridPart = new GR_GridPartType(disp->getGrid());//,*indexSet);
        gridPartStack.push(gridPart);
        globalSpace = new GR_DiscFuncSpaceType (*gridPart);
      }
      space = globalSpace;
    }
    assert(space != 0);

    DATAINFO * dinf = info[n].datinf;
    while (dinf) 
    {
      readFuncData( *disp , *space , path, dinf , time, ntime);
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
  typedef CombinedGrapeDisplay < GrapeDispType > CombinedDisplayType; 
  CombinedGrapeDisplay < GrapeDispType > * comdisp = new CombinedDisplayType ();
  
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
        std::string newpath (path);

        if(numProcs > 1) 
        {
          char procstr[128]; 
          sprintf(procstr,"%d",proc);
          newpath += "_"; 
          newpath += procstr; 
        }
        else
        {
          //newpath += "_-1";
        }

      	std::cout << "NewPath = "<<newpath << std::endl;

        GrapeDispType *newdisp = 0;
        int anz = (n > 0) ? n : 1;
        for(int i=0; i<anz; i++)
        {
          newdisp = readGrid( newpath.c_str(), info[i].name, t_act , ntime, proc );
          assert(newdisp != 0);
          assert( comdisp );
          info = makeData(newdisp,info,newpath.c_str(),info[i].name,t_act, n ,ntime,proc);
          comdisp->addDisplay( *newdisp );
        }
        newdisp->addMyMeshToTimeScene(info[0].tsc,t_act,proc);
        assert( comdisp );
      }

      // this is the combine object, which is put to the last time scene 
      comdisp->addMyMeshToGlobalTimeScene(t_act,0);
    }
    else
    {
      info = makeData(disp,info,path,info[0].name,t_act,n,ntime,numProcs);
      disp->addMyMeshToTimeScene(info[0].tsc,t_act, -1 );
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
 
void dataDispErrorExit(std::string msg) 
{
  std::cerr << msg << std::endl; 
  std::cerr.flush();
  exit(EXIT_FAILURE);
}

void deleteAllObjects() 
{
  //deleteObjects(funcStack);
  if(globalSpace) delete globalSpace;
  deleteObjects(fsStack);
//  deleteObjects(indexStack);
  deleteObjects(gridPartStack);
  deleteObjects(dispStack);
  deleteObjects(gridStack);
  return ;
}

void readDataInfo(std::string path, DATAINFO * dinf, bool parallel = false) 
{
  char dummy[2048];
  std::cout << "Reading data base for " << dinf->name << "! \n";
  std::string dataname(path); 
  if(parallel) 
  { 
    sprintf(dummy,"_0"); 
    dataname += dummy;
  }
  dataname += "/"; 
  dataname += dinf->name;

  int fakedata = 1;
  bool fake = readParameter(dataname,"Fake_data",fakedata);
  if((!fake) || (!fakedata) )
  {
    dinf->base_name = dinf->name; 
    dinf->dimVal = 1;
    dinf->comp = new int [1];
    dinf->comp[0] = 0;
  }
  else
  {
    readParameter(dataname,"DataBase",dummy);
    std::string * basename = new std::string (dummy);
    dinf->base_name = basename->c_str();

    int dimrange;
    readParameter(dataname,"Dim_Range",dimrange);
    if(dimrange <= 0) dataDispErrorExit("wrong dimrange");

    int dimVal = 1;
    readParameter(dataname,"DimVal",dimVal);
    if((dimVal <= 0) || (dimVal > dimrange)) dataDispErrorExit("wrong DimVal");
    dinf->dimVal = dimVal;

    int * comp = new int [dimVal];
    for(int k=0; k<dimVal; k++)
    {
      sprintf(dummy,"%d",k);
      std::string compkey ("comp_");
      compkey += dummy;
      bool couldread = readParameter(dataname,compkey.c_str(),comp[k]);
      if(!couldread) dataDispErrorExit("wrong " + compkey);
    }
    dinf->comp = comp;
  }
  return ;
}

#endif
