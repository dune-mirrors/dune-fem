//**************************************************************
//  (C) written and directecd by Robert Kloefkorn 
//**************************************************************
#ifndef __READ_PARAMS_CC__
#define __READ_PARAMS_CC__

#include "leng.hh"

int readParameterList (int argc, char **argv)
{
  int   i, i_start, i_end;
  INFO * info = 0;
  int    n = 0, n_info = 10;

  int    i_delta = 1;
  const  char *path = 0;
  const  char *replay = 0;
  int    parallel = 1;
  bool   paravis = false;
  
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
      for (int df=0;df<Dune::length<GR_DiscFuncType>::value;++df) {
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
    else if (!strcmp(argv[i], "-pg"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -pg `number of procs'\n");
      parallel = atoi(argv[i+1]);
      i += 2;
      paravis = true;
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
        dataDispErrorExit("usage: -p `path'\n");
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
  
  for(int k=0; k<n; k++) 
  {
    int df = 0;
    DATAINFO * dinf = info[k].datinf; 
    while ( dinf ) 
    {
      if(!path) path = "./";
      readDataInfo(path,dinf,df,paravis);
      assert(dinf->comp);
      dinf = dinf->next;
      ++df;
    }
  }
  
  timeSceneInit(info, n , parallel );
  readData(info, path,i_start,i_end,i_delta,n,timestep,parallel);
  

  // run grape 
  displayTimeScene(info,parallel);
  
  if(replay) 
  {
    std::string cmd("rm manager.replay");
    system(cmd.c_str());
  }

  //deleteAllObjects();

  return (EXIT_SUCCESS);
}
#endif
