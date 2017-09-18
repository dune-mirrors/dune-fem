#ifndef DUNE_FEM_MISC_GRIDNAME_HH
#define DUNE_FEM_MISC_GRIDNAME_HH

// C++ includes
#include <cstdlib>
#include <climits>
#include <iostream>
#include <string>
#include <typeinfo>
#include <vector>

// dune-common includes
#include <dune/common/exceptions.hh>


namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------

    template< class Grid > struct GridName;



    // gridName
    // --------

    template< class GridImp >
    static const std::string &gridName ()
    {
      return GridName< GridImp >::str();
    }

    template< class GridImp >
    static const std::string &gridName ( const GridImp &grid )
    {
      return gridName< GridImp >();
    }



    // UnknownGridException
    // --------------------

    class UnknownGridException : public Exception {};



    // GridName
    // --------

    template< class GridImp >
    struct GridName
    {
      static const std::string &str ()
      {
        static std::string str = computeString();
        return str;
      }

    private:
      static std::string computeString ()
      {
        std::string name( typeid( GridImp ).name() );

        size_t dunePos = name.find( "Dune" );
        name.erase( 0, dunePos+4 );

        char *endptr = 0;
        // get position of strings that are not numbers
        long int result = std::strtol( name.c_str(), &endptr, 0 );
        if( result == LONG_MAX || result == LONG_MIN )
          DUNE_THROW( UnknownGridException, "GridName: faild to determine name of grid!" );

        if( endptr )
          name = std::string( endptr );

        // Grid seems to be followed by IL
        size_t pos = name.find( "GridI" );
        pos += 4; // add length of Grid to get pos of IL

        if( pos < name.size() )
          name.erase( pos, name.size() - pos );
#ifndef NDEBUG
        std::vector< std::string > knownGrids;
        knownGrids.push_back( "AlbertaGrid" );
        knownGrids.push_back( "ALUConformGrid" );
        knownGrids.push_back( "ALUCubeGrid" );
        knownGrids.push_back( "ALUGrid" );
        knownGrids.push_back( "ALUSimplexGrid" );
        knownGrids.push_back( "CacheItGrid" );
        knownGrids.push_back( "CartesianGrid" );
        knownGrids.push_back( "GeometryGrid" );
        knownGrids.push_back( "OneDGrid" );
        knownGrids.push_back( "P4estGrid" );
        knownGrids.push_back( "ParallelGrid" );
        knownGrids.push_back( "ParallelSimplexGrid" );
        knownGrids.push_back( "PrismGrid" );
        knownGrids.push_back( "SPGrid" );
        knownGrids.push_back( "UGGrid" );
        knownGrids.push_back( "YaspGrid" );

        bool found = false ;
        for( size_t i=0; i < knownGrids.size(); ++i )
        {
          if( name == knownGrids[ i ] )
          {
            found = true;
            break;
          }
        }

        if( !found )
        {
          std::cerr << "WARNING: Grid name `" << name << "' not found in list of known grids! Please add in file " << __FILE__ << std::endl;
        }
#endif
        return name;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MISC_GRIDNAME_HH
