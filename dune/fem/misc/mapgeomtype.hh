#ifndef DUNE_FEM_MISC_MAPGEOMTYPE_HH
#define DUNE_FEM_MISC_MAPGEOMTYPE_HH

//- C++ includes
#include <vector>

//- dune-geometry includes
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>


namespace Dune
{

  namespace Fem
  {

    // MapGeometryType
    // ---------------

    template< int dim, class T >
    struct MapGeometryType
    {
      /** \brief dimension */
      static const int dimension = dim;

      /** \brief constructor */
      MapGeometryType ( const T &value = T() )
      {
        const size_t size = Dune::LocalGeometryTypeIndex::size( dimension - 1 );
        data_.resize( size, 0 );
      }

      /** \brief return const reference to data */
      const int &operator[] ( const GeometryType &type ) const
      {
        const size_t index = Dune::LocalGeometryTypeIndex::index( type );
        return data_[ index ];
      }

      /** \brief return const reference to data */
      int &operator[] ( const GeometryType &type )
      {
        return const_cast< int & >( static_cast< const MapGeometryType & >( *this ).operator[]( type ) );
      }

    private:
      std::vector< T > data_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MISC_MAPGEOMTYPE_HH
