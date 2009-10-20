/**************************************************************************
**       Title: 
**    $RCSfile$
**   $Revision$$Name$
**       $Date$
**   Copyright: GPL Author: robertk
** Description: View of an arbitrary grid.
**              Either a DGF-File is required as parameter or a 
**              Alberta/UG/ALUGrid macrofile. An addditional integer can be passed indicating the 
**              additional refinement levels.  
**
**              Select gridtype at compile time by
**
**              make GRIDTYPE=ALBERTAGRID GRIDDIM=2
**
**              or one of UGGRID,ALUGRID_SIMPLEX,ALUGRID_CUBE, YASPGRID, SGRID
**
**              GRIDTYPE=YASPGRID:
**                 compiles and displays dgf-file
**              GRIDTYPE=SGRID:
**                 compiles and displays dgf-file
**              GRIDTYPE=ALBERTAGRID:
**                 compiles and displays dgf-file and ALberta-macrofile
**              GRIDTYPE=ALUGRID_SIMPLEX:
**                 compiles and displays dgf-file and ALugrid-macrofile
**              GRIDTYPE=ALUGRID_CUBE GRIDDIM=3:
**                 compiles, but throws error in display of dgf-file, successful ALugrid-macrofile
**              GRIDTYPE=UGGRID:
**                 error in compilation
** 
**              parallel grids are not tested/supported
**              
**-------------------------------------------------------------------------
**
**  $Log$
**  Revision 1.10  2006/09/25 13:42:57  haasdonk
**  adopted to dfg-parser and grid-selection as make-arguments
**
**
**************************************************************************/

#include <iostream>
#include <string>
#include <config.h>
#include <dune/common/stdstreams.cc>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapegriddisplay.hh>
#endif

using namespace Dune;
using namespace std;

int main (int argc, char **argv)
{
// check whether dgf-file
#if HAVE_GRAPE
  if(argc < 2)
  {
    fprintf(stderr,"usage: %s <MacroFile> \n",argv[0]);
    abort();
  }

  bool dgf_argument = false;
  std::string macroGridName(argv[1]);
  
  string::size_type pos;
  pos = macroGridName.find(".dgf",0);
  if (pos!=string::npos)
      dgf_argument = true;
  
  
  if (dgf_argument) // all gridtypes allowed for dgf-file
  { 
    int level = 0;
    if(argc == 3) level = atoi(argv[2]);
    GridPtr<GridType> gridptr(macroGridName.c_str()); 
    gridptr->globalRefine (level);    
    GrapeGridDisplay < GridType > grape(*gridptr,-1);
    grape.display();
  }
  
  else // specialized macrofiles expected
  {
    
#if defined SGRID || YASPGRID
    fprintf(stderr,"usage: %s <DGF-Macrofile> in ase of SGRID or YASPGRID\n",
            argv[0]);
    abort();   
#endif

#if defined UGGRID
    GridType grid(200,30);
    int level = 0;
    if(argc == 3)
        level = atoi(argv[2]);
    
    AmiraMeshReader< GridType > am;
    am.read ( grid, argv[1] );
    grid.globalRefine (level);
    AmiraMeshWriter< GridType > amw;
    std::string fn (  "testgrid.am" );
    amw.writeGrid ( grid ,fn );
#endif
    
#if defined ALBERTAGRID
    int level = 0;
    if(argc == 3)
        level = atoi(argv[2]);
    GridType grid( argv[1] );
    grid.globalRefine (level);
#endif
    
#if defined ALUGRID_SIMPLEX || ALUGRID_CUBE
    int level = 0;
    if(argc == 3)
        level = atoi(argv[2]); 
#ifdef _ALU3DGRID_PARALLEL_
    GridType grid( argv[1], MPI_COMM_WORLD );
#else
    GridType grid( argv[1] );
#endif
    grid.globalRefine (level);
#endif
    
#if defined ALBERTAGRID || UGGRID || ALUGRID_SIMPLEX || ALUGRID_CUBE
    GrapeGridDisplay < GridType > grape(grid,-1);
    grape.display();
#endif
  }    
#else 
    std::cerr << "Grape was not found! Reconfigure! \n";
#endif // HAVE_GRAPE
  return 0;
}

