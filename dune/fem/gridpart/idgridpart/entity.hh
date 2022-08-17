#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ENTITY_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_ENTITY_HH

#include <type_traits>
#include <utility>

//- dune-grid includes
#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>

//- dune-fem includes
#include <dune/fem/gridpart/common/defaultgridpartentity.hh>
#include <dune/fem/gridpart/idgridpart/geometry.hh>

namespace Dune
{
  namespace Fem {

    // IdEntityBasic
    // -------------

    template< int codim, int dim, class GridFamily >
    class IdEntityBasic
      : public DefaultGridPartEntity < codim, dim, GridFamily >
    {
    protected:
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

    public:
      /** \name Attributes
       *  \{ */

      //! codimensioon of the entity
      static const int codimension = codim;
      //! dimension of the grid
      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      //! dimension of the entity
      static const int mydimension = dimension - codimension;
      //! dimension of the world
      static const int dimensionworld = std::remove_const< GridFamily >::type::dimensionworld;

      /** \} */

      /** \name Types Required by DUNE
       *  \{ */

      //! coordinate type of the grid
      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      //! type of corresponding entity seed
      typedef typename GridFamily::template Codim< codimension >::EntitySeed EntitySeedType;
      //! type of corresponding geometry
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;

      /** \} */

      // type of the host grid
      typedef typename Traits::HostGridPartType  HostGridPartType;

    protected:
      // type of extra data, e.g. a pointer to grid (here empty)
      typedef typename Traits::ExtraData ExtraData;

    public:
      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;

      /** \} */

      /** \name Construction, Initialization and Destruction
       *  \{ */

      /** \brief construct a null entity */
      IdEntityBasic () = default;

      /** \brief construct an initialized entity
       *
       *  \param[in]  hostEntity  corresponding entity in the host grid
       */
      IdEntityBasic ( ExtraData data, HostEntityType hostEntity )
      : data_( std::move( data ) ),
        hostEntity_( std::move( hostEntity ) )
      {}

      /** \} */

      /** \name Methods Shared by Entities of All Codimensions
       *  \{ */

      /** \brief obtain the name of the corresponding reference element
       *
       *  This type can be used to access the DUNE reference element.
       */
      GeometryType type () const
      {
        return hostEntity().type();
      }

      /** \brief obtain the level of this entity */
      int level () const
      {
        return hostEntity().level();
      }

      /** \brief obtain the partition type of this entity */
      PartitionType partitionType () const
      {
        return hostEntity().partitionType();
      }

      /** obtain the geometry of this entity */
      Geometry geometry () const
      {
        return Geometry( hostEntity().geometry() );
      }

      /** \brief return EntitySeed of host grid entity */
      EntitySeedType seed () const { return hostEntity().seed(); }

      /** \brief check for equality */
      bool equals ( const IdEntityBasic &rhs ) const
      {
        return hostEntity() == rhs.hostEntity();
      }

      /** \} */


      /** \name Methods Supporting the GridPart Implementation
       *  \{ */

      const HostEntityType &hostEntity () const
      {
        return hostEntity_;
      }

      const ExtraData &data () const { return data_; }

      /** \} */

    protected:
      ExtraData data_;
      HostEntityType hostEntity_;
    };



    // IdGridEntity
    // ------------

    template< int codim, int dim, class GridFamily >
    class IdEntity : public IdEntityBasic< codim, dim, GridFamily >
    {
      typedef IdEntityBasic< codim, dim, GridFamily > BaseType ;
    protected:
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

    public:
      // type of the host grid
      typedef typename Traits::HostGridPartType  HostGridPartType;
    protected:
      // type of extra data, e.g. a pointer to grid (here empty)
      typedef typename Traits::ExtraData ExtraData;

    public:
      using BaseType :: codimension ;

      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;

      /** \} */

      /** \brief construct a null entity */
      IdEntity () = default;

      IdEntity ( ExtraData data, HostEntityType hostEntity )
      : BaseType( data, hostEntity )
      {}
    };


    // IdGridEntity for codimension 0
    // ----------------------------------

    /** \copydoc IdGridEntity
     *
     *  \nosubgrouping
     */
    template< int dim, class GridFamily >
    class IdEntity< 0, dim, GridFamily > : public IdEntityBasic< 0, dim, GridFamily >
    {
      typedef IdEntityBasic< 0, dim, GridFamily > BaseType ;
    protected:
      typedef typename BaseType::Traits Traits;

      // type of extra data, e.g. a pointer to grid (here empty)
      typedef typename BaseType::ExtraData ExtraData;

    public:
      typedef typename BaseType::HostGridPartType HostGridPartType;

      using BaseType::codimension ;
      using BaseType::data ;
      using BaseType::hostEntity ;
      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;
      /** \} */

    public:
      /** \name Types Required by DUNE
       *  \{ */

      //! type of corresponding local geometry
      typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;

      /** \} */

      /** \name Construction, Initialization and Destruction
       *  \{ */

      /** \brief construct a null entity
       *  \param[in]  data  data pointer (here empty)
       */
      IdEntity () = default;

      /** \brief construct an initialized entity
       *
       *  \param[in]  data        data pointer (here empty)
       *  \param[in]  hostEntity  corresponding entity in the host grid
       */
      IdEntity ( ExtraData data, HostEntityType hostEntity )
      : BaseType( data, hostEntity )
      {}

      /** \} */

      unsigned int subEntities( const unsigned int codim ) const
      {
        return hostEntity().subEntities( codim );
      }

      template< int codim >
      int count () const
      {
        return hostEntity().template count< codim >();
      }

      template< int codim >
      typename Traits::template Codim< codim >::Entity
      subEntity ( int i ) const
      {
        typedef typename Traits::template Codim< codim >::Entity::Implementation EntityImpl;
        return EntityImpl( data(), hostEntity().template subEntity< codim >( i ) );
      }

      bool hasBoundaryIntersections () const
      {
        return hostEntity().hasBoundaryIntersections();
      }

      /** \} */

    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_IDGRID_ENTITY_HH
