//**************************************************************
//  (C) written and directecd by Robert Kloefkorn
//**************************************************************
#ifndef GRAPE_READ_PARAMS_CC
#define GRAPE_READ_PARAMS_CC

//- system includes
#include <string>
#include <sys/types.h>
#include <dirent.h>
#include <tuple>

//- Dune includes
#include <dune/fem/io/file/iointerface.hh>
#include <dune/fem/io/file/datawriter.hh>

inline bool readDataInfo(std::string path, DATAINFO * dinf,
                         const int timestamp, const int dataSet)
{
  static const bool useRankPath = Dune::Fem::DataWriterParameters().separateRankPath();
  std::cout << "Reading data base for " << dinf->name << "! \n";
  std::string dataname;
  if( useRankPath )
  {
    dataname = Fem::IOTupleBase::dataName(
                  Fem::IOInterface::createRecoverPath(path,0, dinf->name, timestamp, useRankPath ),
                  dinf->name);
  }
  else
  {
    return true ;
  }

  {
    std::stringstream dummy;
    dummy << dataSet;
    dataname += "_";
    dataname += dummy.str();
  }

  std::cerr << "reading dofs from: " << dataname << std::endl;

  std::ifstream check ( dataname.c_str() );
  if( ! check )
  {
    std::cerr << "Removing non-valid data set `" << dataname << "'\n";
    // comp = 0 marks non-valid data set
    dinf->comp = 0;
    return false;
  }

  int fakedata = 1;
  bool fake = Fem::readParameter(dataname,"Fake_data",fakedata);

  std::cerr << "FAKE: " << fake << " " << fakedata << std::endl;
  if( (!fake) || (!fakedata) )
  {
    std::string dummy;
    Fem::readParameter(dataname,"DataBase",dummy);
    std::string * basename = new std::string (dummy);
    std::cout << "Read Function: " << *basename << std::endl;
    dinf->base_name = basename->c_str();
    dinf->name = basename->c_str();
    dinf->dimVal = 1;
    if (!dinf->comp)
    {
       dinf->comp = new int [1];
    }
    dinf->comp[0] = 0;
  }
  else
  {
    std::string dummy;
    Fem::readParameter(dataname,"DataBase",dummy);
    std::string * basename = new std::string (dummy);
    std::cout << "Read Function: " << *basename << std::endl;
    dinf->base_name = basename->c_str();

    int dimrange;
    Fem::readParameter(dataname,"Dim_Range",dimrange);
    if(dimrange <= 0) dataDispErrorExit("wrong dimrange");

    int dimVal = 1;
    Fem::readParameter(dataname,"Dim_Domain",dimVal);
    if((dimVal <= 0) || (dimVal > dimrange)) dataDispErrorExit("wrong DimVal");
    dinf->dimVal = dimVal;

    int * comp = new int [dimVal];
    for(int k=0; k<dimVal; k++)
    {
      std::stringstream tmpDummy;
      tmpDummy << k;

      std::string compkey ("comp_");
      compkey += tmpDummy.str();

      bool couldread = Fem::readParameter(dataname,compkey.c_str(),comp[k]);
      if(!couldread) dataDispErrorExit("wrong " + compkey);
    }
    dinf->comp = comp;
  }
  return true;
}

// return number of procs of data set
inline int scanProcsPaths(const std::string globalPath,
                   const std::string dataPrefix,
                   int step)
{
  int procs = 0;
  while ( true )
  {
    std::string path( Fem::IOInterface::
        createRecoverPath(globalPath,procs,dataPrefix,step) );

    // check for directory
    if( ! Dune::Fem::directoryExists( path ) )
    {
      return procs;
    }
    ++procs;
  }

  return procs;
}

// return number of procs of data set
inline int scanProcsFiles(const std::string globalPath,
                          const std::string dataPrefix,
                          int step)
{
  const bool singleBackupRestoreFile = Dune::Fem::PersistenceManager :: singleBackupRestoreFile ;

  int procs = 0;
  std::string path( Fem::IOInterface::
      createRecoverPath(globalPath,procs,dataPrefix,step, false ) );

  while ( true )
  {
    std::stringstream filename;
    filename << path << "/" << dataPrefix << "." << procs;

    // check for file
    const bool fileOk = Dune::Fem::fileExists( filename.str() ) ;

    if( (singleBackupRestoreFile && fileOk) || (! singleBackupRestoreFile && ! fileOk ) )
    {
      return procs;
    }
    ++procs;
  }

  return procs;
}

std::string path,solprefix;

inline int readParameterList (int argc, char **argv, bool displayData = true )
{
  int   i, i_start, i_end;
  INFO * info = 0;
  int    n = 0, n_info = 10;

  int    i_delta = 1;
  double   timestep = 1.0e-3;
  //bool fixedMesh = false;
#ifdef USE_GRAPE_DISPLAY
  const  char *replay = 0;
#endif

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
    path = Dune::Fem::Parameter::commonOutputPath();
    std::string dummyfile;
    Dune::Fem::Parameter::get("fem.io.datafileprefix",solprefix);
    info[n].name = solprefix.c_str();
    info[n].datinf = 0;
    info[n].fix_mesh = 0;
    const int tupleSize = std::tuple_size<GR_DiscFuncType>::value ;
    for (int df = 0; df < tupleSize; ++df )
    {
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
      const int tupleSize = std::tuple_size<GR_DiscFuncType>::value ;
      for (int df = 0; df < tupleSize; ++df )
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
    /*
    else if (!strcmp(argv[i], "-f"))
    {
      fixedMesh = true;
      i += 1;
    }
    */
    else if (!strcmp(argv[i], "-p"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -p `path'\n");
      path = argv[i+1];
      i += 2;
    }
#ifdef USE_GRAPE_DISPLAY
    else if (!strcmp(argv[i], "-replay"))
    {
      if (i+1 == argc)
        dataDispErrorExit("usage: -replay `manager.replay file'\n");
      replay = argv[i+1];
      i += 2;
    }
#endif
    else
    {
      std::cerr << "unknown option " << argv[i] << std::endl;
      exit(1);
    }
    printf("i = %d, argc = %d\n", i, argc);
  }

  // defined in readiotupledata.cc
#ifdef USE_GRAPE_DISPLAY
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
/*
  if( fixedMesh )
  {
    for(int j=0; j<n; ++j)
    {
      info[j].fix_mesh = 1;
    }
  }
*/
#endif

  // scan for max number of processor paths
  int numberProcessors = 0;
  static const bool useRankPath = Dune::Fem::DataWriterParameters().separateRankPath();
  for(int k=0; k<n; k++)
  {
    // scan for max number of processor paths
    int para = ( useRankPath ) ?
          scanProcsPaths(path,info[k].name,i_start) :
          scanProcsFiles(path,info[k].name,i_start) ;

    std::cout << "*******************************************" << std::endl;
    std::cout << "***   Start reading data for " << para << " procs." << std::endl;
    std::cout << "*******************************************" << std::endl;

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

#ifdef USE_GRAPE_DISPLAY
  // initialize time scenes
  timeSceneInit(info, n , numberProcessors);
#endif

  // read all data
  readData(info, path.c_str(),i_start,i_end,i_delta,n,timestep,numberProcessors);

#ifdef USE_GRAPE_DISPLAY
  if( displayData )
  {
    std::cout << "Displaying data of " << numberProcessors << " processors! \n";
    // run grape
    displayTimeScene(info,numberProcessors);
  }

  if(replay)
  {
    std::string cmd("rm manager.replay");
    int result = system(cmd.c_str());
    if( result != 0 )
      std::cerr << "System call possibly failed!" << std::endl;
  }
#endif

  return (EXIT_SUCCESS);
}
#endif
