#ifndef DUNE_GRCOMMON_HH
#define DUNE_GRCOMMON_HH

#if HAVE_GRAPE
#include <dune/grid/io/visual/grape/grapewrapper.hh>
// now include the necessary file
#include <dune/grid/io/visual/grape/grapecommon.hh>
#else

// if grape is not available do some dummy typedefing
typedef int BUTTON;
typedef int COMBOBUTTON;
typedef int TIMESCENE;
typedef int MANAGER;

/* info about data on one mesh */
typedef struct datainfo DATAINFO;
struct datainfo
{
  const char * name;
  const char * base_name;
  DATAINFO *next;

  int dimVal; /* length of vector (dimVal = 1 --> scalar, otherwise vector  */
  int * comp; /* number of each component */
};

/* info about one mesh */
typedef struct info INFO;
struct info
{
  int fix_mesh; /* if no dynamic grid 1 : else 0 */
  const char  *name;
  DATAINFO * datinf;
  void  *tsc;
};
#endif

#endif
