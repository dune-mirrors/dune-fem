#ifndef DUNE_FEM_ALLGEOMTYPES_HH
#define DUNE_FEM_ALLGEOMTYPES_HH

//- system includes
#include <vector>
#include <set>

//- Dune includes
#include <dune/fem/common/forloop.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/capabilities.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  namespace Fem
  {

    // GeometryInformation
    // -------------------

    /**  \brief ReferenceVolume and local bary center keeper class.
     */
    template< class GridImp, int codim >
    class GeometryInformation
    {
      typedef GeometryInformation< GridImp, codim > ThisType;

    public:
      //! grid type
      typedef GridImp GridType;

      //! dimension
      static const int dim = GridType::dimension - codim;

      //! coordinate type
      typedef typename GridType::ctype ctype;

    protected:
      typedef Dune::ReferenceElements< ctype, dim > ReferenceElementsType;

    public:
      //! type of reference element
      typedef std::decay_t< decltype( ReferenceElementsType::general( std::declval< const GeometryType & >() ) ) > ReferenceElementType;

      //! type of domain vector
      typedef FieldVector<ctype, dim> DomainType;

    protected:
      //! constructor creating empty geometry information
      GeometryInformation () : isNoneLocalCenter_( 0.5 )
      {}

    public:
      //! creating geometry information due to given geometry types list
      explicit GeometryInformation( const std::vector< GeometryType > &geomTypes )
      {
        buildMaps( geomTypes );
      }

      //! return local bary center for geometry of type type
      const DomainType localCenter ( const GeometryType &type ) const
      {
        return type.isNone() ? isNoneLocalCenter_ : referenceElement( type ).position( 0, 0 );
      }

      //! return volume of reference element for geometry of type type
      ctype referenceVolume ( const GeometryType &type ) const
      {
        return type.isNone() ? ctype(1.0) : referenceElement( type ).volume();
      }

      //! return reference element for type
      static auto referenceElement ( const GeometryType &type )
        -> decltype( ReferenceElementsType::general( type ) )
      {
        assert( ! type.isNone() );
        return ReferenceElementsType::general( type );
      }

    protected:
      DomainType isNoneLocalCenter_;

      //! build maps
      void buildMaps ( const std::vector< GeometryType > &geomTypes )
      {}
    };


    /**  \brief default implementation uses method geomTypes of given index
         set. Used in DiscreteFunctionSpaces.
     */
    template< class IndexSetImp, class GridImp >
    class AllGeomTypes : public GeometryInformation< GridImp , 0>
    {
    public:
      typedef IndexSetImp IndexSetType;
      typedef GridImp     GridType;

    private:
      typedef AllGeomTypes< IndexSetType, GridType > ThisType;
      static const unsigned int ncodim = GridType :: dimension + 1;

    protected:
      std::vector< std::vector< GeometryType > > geomTypes_;

    public:
      // insert all types of the reference element into the geomTypes list
      template <int dim>
      struct InsertGeometryTypes
      {
        static void apply( std::vector< std::set< GeometryType > >& geomTypes )
        {
          static const int codim = GridType :: dimension - dim ;
          typedef Dune::ReferenceElements< typename GridType :: ctype, dim > ReferenceElementContainer ;
          typedef typename ReferenceElementContainer :: Iterator Iterator ;
          for( Iterator it = ReferenceElementContainer::begin(),
                       end = ReferenceElementContainer::end(); it != end; ++it )
          {
            geomTypes[ codim ].insert( it->type() );
          }
        }
      };

      //! constructor storing index set reference
      inline explicit AllGeomTypes( const IndexSetType &indexSet )
        : geomTypes_( ncodim )
      {
        if( multipleGeomTypes() )
        {
          std::vector< std::set< GeometryType > > geomTypes( ncodim );

          // store all possible geom types
          Fem::ForLoop< InsertGeometryTypes, 0, GridType::dimension > :: apply( geomTypes );

          auto types = indexSet.types( 0 );
          for( const auto type : types )
          {
            geomTypes[ 0 ].insert( type );
          }

          for( unsigned int cd=0; cd < ncodim; ++cd )
          {
            geomTypes_[ cd ].reserve( geomTypes[ cd ].size() );
            for( const auto& type : geomTypes[ cd ] )
            {
              geomTypes_[ cd ].push_back( type );
            }
          }
        }
        else
        {
          // single geometry type
          for( int codim=0; codim<GridType::dimension+1; ++codim )
          {
            typename IndexSetType::Types types = indexSet.types( codim );
            const int size = types.size();
            geomTypes_[ codim ].resize( size );
            std::copy_n( types.begin(), size, geomTypes_[ codim ].begin() );
          }
          // build geometry information for codim 0
          this->buildMaps( geomTypes_[ 0 ] );
        }
      }

      //! returns vector with geometry tpyes this index set has indices for
      const std :: vector< GeometryType > &geomTypes ( unsigned int codim ) const
      {
        assert( codim < ncodim );
        return geomTypes_[ codim ];
      }

      //! UGGrid might have different geom types
      static bool multipleGeomTypes ()
      {
        return !Dune::Capabilities :: hasSingleGeometryType < GridType > :: v;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALLGEOMTYPES_HH
