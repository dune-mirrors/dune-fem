//**************************************************************
//  (C) written and directecd by Robert Kloefkorn 
//**************************************************************
#ifndef __READ_FILE_CC__
#define __READ_FILE_CC__

#include <dune/fem/io/file/asciiparser.hh>
int readParameterFile (int argc, char **argv)
{
  if(argc != 2) 
  {
    std::cerr << "ERROR: wrong parameters !\n";
    abort();
  }

  FILE * file = fopen (argv[1],"r");
  assert( file );
  fclose(file);

  // name of parameter file 
  std::string filename(argv[1]);
 
  char dummy[2048]; 
  readParameter(filename,"Path",dummy);
  std::string path (dummy);

  int i_start, i_end;
  int i_delta = 1;
  int parallel = 0;
  
  readParameter(filename,"StartStep",i_start);
  readParameter(filename,"EndStep",i_end);
  readParameter(filename,"SkipStep",i_delta);
  bool procs = readParameter(filename,"Processors",parallel);
  // proc + 1 
  parallel++ ;
  
  if(!procs) 
    parallel = 0;
  
  readParameter(filename,"GridPrefix",dummy);
  std::string gridpref (dummy);
  std::vector < std::string > datapref(100); 
  int count = 0;
  {
    bool success = readParameter(filename,"DataPrefix_0",dummy);
    // we might have more than one data prefix     
    while ( success )
    {
      datapref[count] += dummy; 
      count++;
      
      std::string keyword ("DataPrefix_");
      sprintf(dummy,"%d",count);
      keyword += dummy;

      success = readParameter(filename,keyword.c_str(),dummy);
    }
  }
  
  INFO * info = 0;
  int    n = 0, n_info = 10;

  REAL   timestep = 1.0e-3;

  info = (INFO *) malloc(n_info*sizeof(INFO));
  assert(info != 0);
  
  info[0].datinf = 0;
  info[0].name = "grid";
  info[0].fix_mesh = 0;

  info[n].name = gridpref.c_str();
  info[n].datinf = 0;
  info[n].fix_mesh = 0;
  n++;
  
  for(int i=0; i<count; i++)
  {
    DATAINFO * dinf = (DATAINFO *) std::malloc(sizeof(DATAINFO));
    assert(dinf);

    std::string dataname(path); 
    dataname += "/"; dataname += datapref[i];
    dinf->name = datapref[i].c_str();
    dinf->comp = 0;
    dinf->dimVal = 0;
    
    readDataInfo(path,dinf, (parallel > 0));
    assert(dinf->comp);

    /* seems wrong order, but grape truns it arround, we can do nothing else here */
    dinf->next = info[n].datinf; 
    info[n].datinf = dinf;
  }

  timeSceneInit(info, n , parallel);
  readData(info, path.c_str(),i_start,i_end,i_delta,n,timestep,parallel);
  
  // run grape 
  displayTimeScene(info,parallel);

  deleteAllObjects();
  return (EXIT_SUCCESS);
}
#endif
