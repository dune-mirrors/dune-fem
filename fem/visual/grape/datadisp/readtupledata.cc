//**************************************************************
//  (C) written and directecd by Robert Kloefkorn 
//**************************************************************
#ifndef __DATADISP_READDATA_CC__
#define __DATADISP_READDATA_CC__

#include <stack> 

#define LARGE 1.0E308

#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/file/grapedataio.hh>
#include "../grapetuple.hh"

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
  gridStack.push(grid);

  assert(filename);
  std::string fn (filename);
  DATAINFO * dinf = info[n].datinf;

  GrapeTuple<GR_DiscFuncType>::ReturnType* tup = 
    GrapeTuple<GR_DiscFuncType>::input(dataIO,grid,time,ntime,path,fn);
  std::cout << "Finished reading grid" << std::endl;
  
  GrapeDispType * disp = new GrapeDispType ( *grid, myRank );  
  //GrapeDispType * disp = 
  //  new GrapeDispType ( tup->first()->getFunctionSpace().gridPart(), myRank );  
  dispStack.push(disp);
  GrapeTuple<GR_DiscFuncType>::addToDisplay(*disp,dinf,time,*tup);

  addError(*disp,*grid,time,*(tup->first()));
  return disp;
}

// setup the hole data tree for grape 
INFO * readData(INFO * info , const char * path, int i_start, int i_end, 
    int i_delta, int n, double timestep, int numProcs) 
{
  double t_start = LARGE;
  double t_end = -LARGE, t_act=0.0;
  typedef CombinedGrapeDisplay < GrapeDispType > CombinedDisplayType; 
  CombinedGrapeDisplay < GrapeDispType > * comdisp = new CombinedDisplayType ();
  
  int  ntime, n_step = 0;
  
  for (ntime = i_start; ntime <= i_end; ntime += i_delta)
  {
    printf("timestep = %d | last timestep = %d | stepsize = %d\n", ntime, i_end, i_delta);
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
        int anz = n; // (n > 0) ? n : 1;
        for(int i=0; i<anz; i++)
        {
      	  newdisp = readTupleData(newpath.c_str(), info[i].name, 
                                t_act , i , ntime, proc,   info);
          assert(newdisp != 0);
          assert( comdisp );
          comdisp->addDisplay( *newdisp );
        }
	      assert(newdisp);
        newdisp->addMyMeshToTimeScene(info[0].tsc,t_act,proc);
        assert( comdisp );
      }

      // this is the combine object, which is put to the last time scene 
      comdisp->addMyMeshToGlobalTimeScene(t_act,0);
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
  /*
  //deleteObjects(funcStack);
  //if(globalSpace) delete globalSpace;
  //deleteObjects(fsStack);
  //deleteObjects(indexStack);
  deleteObjects(gridPartStack);
  deleteObjects(dispStack);
  deleteObjects(gridStack);
  return ;
  */
}

void readDataInfo(std::string path, DATAINFO * dinf, int k, bool parallel = false) 
{
  char dummy[2048];
  std::cout << "Reading data base for " << dinf->name << "! \n";
  std::stringstream dataname; 
  dataname << path;
  if(parallel) 
  { 
    dataname << "_0";
  }
  dataname << "/d"; 
  dataname << dinf->name;
  dataname << "_" << k;

  std::cerr << "reading dofs from: " << dataname.str() << std::endl;

  int fakedata = 1;
  bool fake = readParameter(dataname.str(),"Fake_data",fakedata);
  std::cerr << "FAKE: " << fake << " " << fakedata << std::endl;
  if( (!fake) || (!fakedata) )
  {
    readParameter(dataname.str(),"DataBase",dummy);
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
    readParameter(dataname.str(),"DataBase",dummy);
    std::string * basename = new std::string (dummy);
    std::cout << "Read Function: " << *basename << std::endl;
    dinf->base_name = basename->c_str();

    int dimrange;
    readParameter(dataname.str(),"Dim_Range",dimrange);
    if(dimrange <= 0) dataDispErrorExit("wrong dimrange");

    int dimVal = 1;
    readParameter(dataname.str(),"Dim_Domain",dimVal);
    if((dimVal <= 0) || (dimVal > dimrange)) dataDispErrorExit("wrong DimVal");
    dinf->dimVal = dimVal;

    int * comp = new int [dimVal];
    for(int k=0; k<dimVal; k++)
    {
      sprintf(dummy,"%d",k);
      std::string compkey ("comp_");
      compkey += dummy;
      bool couldread = readParameter(dataname.str(),compkey.c_str(),comp[k]);
      if(!couldread) dataDispErrorExit("wrong " + compkey);
    }
    dinf->comp = comp;
  }
  return ;
}

#endif
