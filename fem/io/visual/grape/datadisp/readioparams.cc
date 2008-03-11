//**************************************************************
//  (C) written and directecd by Robert Kloefkorn 
//**************************************************************
#ifndef GRAPE_READ_PARAMS_CC
#define GRAPE_READ_PARAMS_CC

//- system includes 
#include <string>
#include <sys/types.h>
#include <dirent.h>

//- Dune includes 
#include <dune/fem/misc/femtuples.hh>
#include <dune/fem/io/file/iointerface.hh>


// return number of procs of data set 
int scanProcsPaths(const std::string globalPath, 
                   const std::string dataPrefix,
                   int step)
{
  int procs = 0;
  while ( true )
  {
    std::string path( IOInterface::
        createRecoverPath(globalPath,procs,dataPrefix,step));
    // if directory does not exist
    // return number of procs 
    if( ! opendir(path.c_str()) )
    {
      return procs;
    }
    ++procs;
  }

  return procs;
}

  std::string path,solprefix;

int readParameterList (int argc, char **argv)
{
  int   i, i_start, i_end;
  INFO * info = 0;
  int    n = 0, n_info = 10;

  int    i_delta = 1;
  const  char *replay = 0;
  REAL   timestep = 1.0e-3;
  
  info = (INFO *) malloc(n_info*sizeof(INFO));
  assert(info != 0);
  
  info[0].datinf = 0;
  info[0].name = "grid";
  info[0].fix_mesh = 0;

  if (argc < 3)
  {
    print_help("datadisp");
    return(0);
  }
  if (argc == 3) {
    path = Parameter::prefix();
    std::string dummyfile;
    Parameter::get("fem.io.datafileprefix",solprefix);
    info[n].name = solprefix.c_str();
    info[n].datinf = 0;
    info[n].fix_mesh = 0;
    for (int df=0;df<TupleLength<GR_DiscFuncType>::value;++df)  {
      DATAINFO * dinf = (DATAINFO *) std::malloc(sizeof(DATAINFO));
      assert(dinf);
      dinf->name = solprefix.c_str();
      dinf->base_name = 0; 
      dinf->comp = 0;
      dinf->dimVal = 0;
      dinf->next = info[n].datinf; 
      info[n].datinf = dinf;
    }
    n++;
  }
  
  i_start = atoi(argv[1]);
  i_end = atoi(argv[2]);

  std::cout << "Read paramaeterlist ... "<< std::endl;

  i = 3;
  while (i < argc)
  {
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help"))
    {
      print_help("datadisp");
    }
    else if (!strcmp(argv[i], "-i"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -i `increment'\n");
      i_delta = atoi(argv[i+1]);
      i += 2;
    }
    else if (!strcmp(argv[i], "-v"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -v `vectorprefix'\n");
      
      DATAINFO * dinf = (DATAINFO *) std::malloc(sizeof(DATAINFO));
      assert(dinf);
      dinf->name = argv[i+1];
      dinf->base_name = argv[i+1];

      dinf->comp = 0;
      dinf->dimVal = 0;

      /* seems wrong order, but grape truns it arround, we can do nothing else here */
      dinf->next = info[n].datinf; 
      info[n].datinf = dinf;

      i += 2;
    }
    else if (!strcmp(argv[i], "-t"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -t `time step size'\n");
      timestep = atof(argv[i+1]);
      i += 2;
    }
    else if (!strcmp(argv[i], "-m"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -m `gridprefix'\n");
      assert(n < n_info);
      info[n].name = argv[i+1];
      info[n].datinf = 0;
      info[n].fix_mesh = 0;
      for (int df=0;df<TupleLength<GR_DiscFuncType>::value;++df) 
      {
        DATAINFO * dinf = (DATAINFO *) std::malloc(sizeof(DATAINFO));
        assert(dinf);
        dinf->name = argv[i+1];
        dinf->base_name = 0; 
        dinf->comp = 0;
        dinf->dimVal = 0;
        dinf->next = info[n].datinf; 
        info[n].datinf = dinf;
      }
      n++;
      i += 2;
    }
    else if (!strcmp(argv[i], "-f"))
    {
      info[n].fix_mesh = 1;
      i += 1;
    }
    else if (!strcmp(argv[i], "-p"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -p `path'\n");
      path = argv[i+1];
      i += 2;
    }
    else if (!strcmp(argv[i], "-replay"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -replay `manager.replay file'\n");
      replay = argv[i+1];
      i += 2;
    }
    else
    {
      std::cerr << "unknown option " << argv[i] << std::endl;
      exit(1);
    }
    printf("i = %d, argc = %d\n", i, argc);
  }
 
  if(replay)
  {
    std::string replayfile(replay);
    // if strcmp > 0 then strins not equal 
    if(replayfile != "manager.replay")
    {
      std::string cmd("ln -s ");
      cmd += replayfile;
      cmd += " manager.replay";

      //std::cout << "call : " << cmd << "\n";
      int result = system(cmd.c_str());

      if(result != 0) replay = 0;
    }
  }
 
  // scan for max number of processor paths  
  int numberProcessors = 0;
  for(int k=0; k<n; k++) 
  {
    // scan for max number of processor paths  
    int para = scanProcsPaths(path,info[k].name,i_start);
    // should be at least 1 
    if( para <= 0 )
    {
      std::cerr << "ERROR: not a valid data path! \n";
      abort();
    }   
    numberProcessors = std::max(numberProcessors,para);
    
    int df = 0;
    DATAINFO * dinf = info[k].datinf; 
    while ( dinf ) 
    {
      if( path == "") path = "./";
      readDataInfo(path, dinf, i_start, df);
      dinf = dinf->next;
      ++df;
    }
  }
  
  // initialize time scenes
  timeSceneInit(info, n , numberProcessors);

  // read all data 
  readData(info, path.c_str(),i_start,i_end,i_delta,n,timestep,numberProcessors);
  
  std::cout << "Displaying data of " << numberProcessors << " processors! \n";
  // run grape 
  displayTimeScene(info,numberProcessors);
  
  if(replay) 
  {
    std::string cmd("rm manager.replay");
    system(cmd.c_str());
  }

  return (EXIT_SUCCESS);
}
#endif
