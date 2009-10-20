#ifndef DUNE_GRCOMMON_HH
#define DUNE_GRCOMMON_HH

#if HAVE_GRAPE 
#include <dune/grid/io/visual/grape/grapewrapper.hh>
#else 

// if grape is not available do some dummy typedefing
typedef int BUTTON;
typedef int COMBOBUTTON;
typedef int TIMESCENE;
typedef int MANAGER;

#endif

// now include the necessary file 
#include <dune/grid/io/visual/grape/grapecommon.hh>
#endif
