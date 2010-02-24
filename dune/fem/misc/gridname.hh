#ifndef DUNE_GRIDNAME_HH
#define DUNE_GRIDNAME_HH

#include <string>



namespace Dune 
{

#ifdef ENABLE_ALBERTA
  template< int dim, int dimW >
  class AlbertaGrid;
#endif
  
#ifdef ENABLE_ALUGRID
  template< int dim, int dimW >
  class ALUCubeGrid;
  template< int dim, int dimW >
  class ALUSimplexGrid;
  template< int dim, int dimW >
  class ALUConformGrid;
#endif

  namespace Fem 
  {
    
    template <class GridImp> 
    struct GridName
    {
      static std::string name(const GridImp& grid) 
      {
        //return std::string( "unknown" );
        return grid.name();
      }
    };

#ifdef ENABLE_ALBERTA
    template <int dim, int dimworld> 
    struct GridName< AlbertaGrid<dim, dimworld> >
    {
      static std::string name(const AlbertaGrid<dim, dimworld>& ) 
      {
        return std::string( "AlbertaGrid" );
      }
    };
#endif

#ifdef ENABLE_ALUGRID
    template <int dim, int dimworld> 
    struct GridName< ALUCubeGrid<dim, dimworld> >
    {
      static std::string name(const ALUCubeGrid<dim, dimworld>& ) 
      {
        return std::string( "ALUCubeGrid" );
      }
    };
    template <int dim, int dimworld> 
    struct GridName< ALUSimplexGrid<dim, dimworld> >
    {
      static std::string name(const ALUSimplexGrid<dim, dimworld>& ) 
      {
        return std::string( "ALUSimplexGrid" );
      }
    };
    template <int dim, int dimworld> 
    struct GridName< ALUConformGrid<dim, dimworld> >
    {
      static std::string name(const ALUConformGrid<dim, dimworld>& ) 
      {
        return std::string( "ALUConformGrid" );
      }
    };
#endif

    template <class GridImp>
    inline std::string gridName(const GridImp& grid) 
    {
      return GridName< GridImp > :: name( grid ); 
    }
  }
}
#endif
