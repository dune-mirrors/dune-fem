//**************************************************************
//  (C) written and directecd by Robert Kloefkorn 
//**************************************************************
#ifndef DATADISP_READIOTUPLE_CC
#define DATADISP_READIOTUPLE_CC

#include <stack> 

#define LARGE 1.0E308

#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/file/grapedataio.hh>
#include <dune/fem/io/file/iotuple.hh>
#include <dune/fem/io/file/iointerface.hh>

void dataDispErrorExit(std::string msg);

static std::stack<GR_GridType *> gridStack;
static std::stack<GrapeDispType *> dispStack;
static std::stack<GR_GridPartType *> gridPartStack;

static GrapeDataIO < GR_GridType> dataIO;

template <class T> 
void deleteObjects(std::stack<T *> & stack);


GrapeDispType * readTupleData(const char * path, const char * filename, 
            double & time , int n, 
            int ntime, int myRank,
            INFO* info)
{
  GR_GridType * grid;

  assert(filename);
  std::string fn (filename);
  DATAINFO * dinf = info[n].datinf;

  IOTuple<GR_DiscFuncType>::ReturnType* tup = 
    IOTuple<GR_DiscFuncType>::input(dataIO,grid,time,ntime,path,fn);
  std::cout << "Finished reading grid" << std::endl;

  gridStack.push(grid);
  
  GrapeDispType * disp = new GrapeDispType ( *grid, myRank );  
  dispStack.push(disp);
  IOTuple<GR_DiscFuncType>::addToDisplay(*disp,dinf,time,*tup);

  addError(*disp,*grid,time,*(tup->first()));
  return disp;
}

// setup the hole data tree for grape 
INFO * readData(INFO * info , const char * path, int i_start, int i_end, 
    int i_delta, int n, double timestep, int numProcs) 
{
  double t_start = LARGE;
  double t_end = -LARGE, t_act = 0.0;
  typedef CombinedGrapeDisplay < GrapeDispType > CombinedDisplayType; 

  
  CombinedGrapeDisplay < GrapeDispType > * comdisp = 0 ;
  if( numProcs > 1 ) comdisp = new CombinedDisplayType ();
  
  int  ntime, n_step = 0;
  
  for (ntime = i_start; ntime <= i_end; ntime += i_delta)
  {
    printf("timestep = %d | last timestep = %d | stepsize = %d\n", ntime, i_end, i_delta);
    {
      int anzProcs = numProcs;
      //if(numProcs > 1) anzProcs--;
      
      for(int proc=0; proc<anzProcs; ++proc)
      {
        assert(path || numProcs <= 1); 
        GrapeDispType *newdisp = 0;

        int anz = n; // (n > 0) ? n : 1;
        for(int i=0; i<anz; ++i)
        {
          // use standard procedure to create path name 
          std::string newpath = 
            IOInterface::createRecoverPath(path,proc,info[i].name,ntime);

          newdisp = readTupleData(newpath.c_str(), info[i].name, 
                                t_act , i , ntime, proc,   info);

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
void deleteObjects(std::stack<T *> & stack) 
{
  while(! stack.empty() )
  {
    T * obj = stack.top();
    stack.pop();
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
  //if(globalSpace) delete globalSpace;
  //deleteObjects(fsStack);
  //deleteObjects(indexStack);
  deleteObjects(gridPartStack);
  deleteObjects(dispStack);
  deleteObjects(gridStack);
  return ;
}

void readDataInfo(std::string path, DATAINFO * dinf, 
    const int dataSet, const int timestamp, const int k) 
{
  char dummy[2048];
  std::cout << "Reading data base for " << dinf->name << "! \n";
  std::string dataname = 
    IOTupleBase::dataName( 
      IOInterface::createRecoverPath(path,0, dinf->name, timestamp),
      dinf->name);

  {
    std::stringstream dummy; 
    dummy << dataSet; 
    dataname += "_";
    dataname += dummy.str();
  }
  
  std::cerr << "reading dofs from: " << dataname << std::endl;

  int fakedata = 1;
  bool fake = readParameter(dataname,"Fake_data",fakedata);
  std::cerr << "FAKE: " << fake << " " << fakedata << std::endl;
  if( (!fake) || (!fakedata) )
  {
    readParameter(dataname,"DataBase",dummy);
    std::string * basename = new std::string (dummy);
    std::cout << "Read Function: " << *basename << std::endl;
    dinf->base_name = basename->c_str();
    dinf->name = basename->c_str();
    dinf->dimVal = 1;
    if (!dinf->comp)
        dinf->comp = new int [1];
    dinf->comp[0] = 0;
  }
  else
  {
    readParameter(dataname,"DataBase",dummy);
    std::string * basename = new std::string (dummy);
    std::cout << "Read Function: " << *basename << std::endl;
    dinf->base_name = basename->c_str();

    int dimrange;
    readParameter(dataname,"Dim_Range",dimrange);
    if(dimrange <= 0) dataDispErrorExit("wrong dimrange");

    int dimVal = 1;
    readParameter(dataname,"Dim_Domain",dimVal);
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
