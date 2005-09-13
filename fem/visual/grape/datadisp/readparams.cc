#ifndef __READ_PARAMS_CC__
#define __READ_PARAMS_CC__

int readParameterList (int argc, char **argv)
{
  int   i, i_start, i_end;
  INFO * info = 0;
  int    n = 0, n_info = 10;

  int    i_delta = 1;
  const  char *path = 0;
  bool   time_bar = false;
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
    else if (!strcmp(argv[i], "-s"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -s `dataprefix'\n");
      
      DATAINFO * dinf = (DATAINFO *) std::malloc(sizeof(DATAINFO));
      assert(dinf);
      dinf->name = argv[i+1];
      dinf->base_name = argv[i+1]; 

      dinf->comp = 0;
      dinf->dimVal = 0;

      dinf->next = info[n].datinf; 
      info[n].datinf = dinf;

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
        dataDispErrorExit("usage: -pg `number of procs'\n");
      parallel += atoi(argv[i+1]);
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
    else
    {
      std::cerr << "unknown option " << argv[i] << std::endl;
      exit(1);
    }
    printf("i = %d, argc = %d\n", i, argc);
  }
  
  
  for(int k=0; k<n+1; k++) 
  {
    DATAINFO * dinf = info[k].datinf; 
    assert(dinf); 
    while ( dinf ) 
    {
      if(!path) path = "./";
      readDataInfo(path,dinf,paravis);
      assert(dinf->comp);
      dinf = dinf->next;
    }
  }
  

  std::cout << "Path = "<< path << std::endl;

  timeSceneInit(info, n , parallel , time_bar);
  readData(info, path,i_start,i_end,i_delta,n,timestep,parallel);
  
  // run grape 
  displayTimeScene(info,parallel);

  deleteAllObjects();
  return (EXIT_SUCCESS);
}
#endif
