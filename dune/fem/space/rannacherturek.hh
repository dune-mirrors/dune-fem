#if ! HAVE_DUNE_LOCALFUNCTIONS

#error RannacherTurekDiscreteFunctionSpace depends \
       on dune-localfunctions. Rebuild your Dune \
       with dune-localfunctions.

#else
#include <dune/fem/space/rannacherturek/space.hh>
#endif
