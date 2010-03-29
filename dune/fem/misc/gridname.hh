#ifndef DUNE_GRIDNAME_HH
#define DUNE_GRIDNAME_HH

#include <string>



namespace Dune 
{

  template< int dimW >
  class YaspGrid;

  template< int dim, int dimW , class ctype >
  class SGrid;

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
      static std::string computeName(const GridImp& grid) 
      {
        std::string name ( typeid( GridImp ).name() ); 

        size_t dunePos = name.find( "Dune" );
        name.erase( 0, dunePos+4 );

        char* endptr = 0;
        // get position of strings that are not numbers 
        strtol( name.c_str(), &endptr, 0 );

        if( endptr ) 
          name = std::string( endptr );

        // Grid seems to be followed by IL 
        size_t pos = name.find( "GridIL" );
        pos += 4; // add length of Grid to get pos of IL 

        if( pos < name.size() ) 
          name.erase( pos, name.size() - pos );

        return name;
      }

      static std::string name(const GridImp& grid) 
      {
        static std::string name ( computeName( grid ) ); 
        return name;
      }
    };

#if 0
    template < int dimworld > 
    struct GridName< YaspGrid<dimworld> >
    {
      static std::string name(const YaspGrid<dimworld>& ) 
      {
        std::string name ( typeid( YaspGrid<dimworld> ).name() ); 

        //return std::string( "YaspGrid" );
      }
    };

    template <int dim, int dimworld, class ctype > 
    struct GridName< SGrid<dim, dimworld, ctype> >
    {
      static std::string name(const SGrid<dim, dimworld, ctype>& ) 
      {
        return std::string( "SGrid" );
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
#endif

    template <class GridImp>
    inline std::string gridName(const GridImp& grid) 
    {
      return GridName< GridImp > :: name( grid ); 
    }
  }
}
#endif
