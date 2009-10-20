#ifndef DUNE_ADAPTCAPS_HH
#define DUNE_ADAPTCAPS_HH

namespace Dune
{

  // Forward Declaration of Grids
  // ----------------------------

  template< int dim, int dimworld >
  class ALUSimplexGrid;

  template< int dim, int dimworld >
  class ALUCubeGrid;

  template< int dim, int dimworld >
  class ALUConformGrid;

#if HAVE_ALBERTA
  template< int dim, int dimworld >
  class AlbertaGrid;
#endif


  // Additional Adaptation Capabilities
  // ----------------------------------

  namespace Capabilities
  {

    template< class Grid >
    struct supportsCallbackAdaptation
    {
      static const bool v = false;
    };

    template< class Grid >
    struct supportsCallbackAdaptation< const Grid >
    {
      static const bool v = Dune :: Capabilities :: supportsCallbackAdaptation< Grid > :: v;
    };

    template< int dim, int dimworld >
    struct supportsCallbackAdaptation< ALUSimplexGrid< dim, dimworld > >
    {
      static const bool v = true;
    };

    template< int dim, int dimworld >
    struct supportsCallbackAdaptation< ALUCubeGrid< dim, dimworld > >
    {
      static const bool v = true;
    };

    template< int dim, int dimworld >
    struct supportsCallbackAdaptation< ALUConformGrid< dim, dimworld > >
    {
      static const bool v = true;
    };


#if HAVE_ALBERTA
    template< int dim, int dimworld >
    struct supportsCallbackAdaptation< AlbertaGrid< dim, dimworld > >
    {
      static const bool v = true;
    };
#endif

  }

}

#endif
