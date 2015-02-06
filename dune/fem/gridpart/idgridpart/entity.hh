#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_ENTITY_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_ENTITY_HH

//- dune-common includes
#include <dune/common/nullptr.hh>

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
      typedef typename remove_const< GridFamily >::type::Traits Traits;

    public:
      /** \name Attributes
       *  \{ */

      //! codimensioon of the entity
      static const int codimension = codim;
      //! dimension of the grid
      static const int dimension = Traits::dimension;
      //! dimension of the entity
      static const int mydimension = dimension - codimension;
      //! dimension of the world
      static const int dimensionworld = Traits::dimensionworld;

      /** \} */

      /** \name Types Required by DUNE
       *  \{ */

      //! coordinate type of the grid
      typedef typename remove_const< GridFamily >::type::ctype ctype;

      //! type of corresponding entity seed
      typedef typename GridFamily::template Codim< codimension >::EntitySeed EntitySeedType;
      //! type of corresponding geometry
      typedef typename Traits::template Codim< codimension >::Geometry Geometry;

      /** \} */

    protected:
      // type of the host grid
      typedef typename Traits::HostGridPartType  HostGridPartType;

      // type of extra data, e.g. a pointer to grid (here empty)
      typedef typename Traits::ExtraData ExtraData;

    public:
      /** \name Host Types
       *  \{ */

      //! type of corresponding host entity
      typedef typename HostGridPartType::template Codim< codimension >::EntityType HostEntityType;
      //! type of corresponding host entity pointer
      typedef typename HostGridPartType::template Codim< codimension >::EntityPointerType HostEntityPointerType;
      /** \} */

      /** \name Construction, Initialization and Destruction
       *  \{ */

      /** \brief construct a null entity */
      explicit IdEntityBasic ( ExtraData data )
      : hostEntity_( nullptr ),
        data_ ( data )
      {}

      /** \brief construct an initialized entity
       *
       *  \param[in]  hostEntity  corresponding entity in the host grid
       *
       *  \note The reference to the host entity must remain valid  as long as
       *        this entity is in use.
       */
      IdEntityBasic ( ExtraData data, const HostEntityType &hostEntity )
      : hostEntity_( &hostEntity ),
        data_( data )
      {}

      /** \} */

      /** \brief return true if entity hold a vaild host entity */
      operator bool () const { return bool( hostEntity_ ); }

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

      /** \} */


      /** \name Methods Supporting the GridPart Implementation
       *  \{ */

      const HostEntityType &hostEntity () const
      {
        assert( *this );
        return *hostEntity_;
      }

      ExtraData data() const { return data_; }

      /** \} */

    protected:
      const HostEntityType *hostEntity_;
      ExtraData             data_;
    };



    // IdGridEntity
    // ------------

    template< int codim, int dim, class GridFamily >
    class IdEntity : public IdEntityBasic< codim, dim, GridFamily >
    {
      typedef IdEntityBasic< codim, dim, GridFamily > BaseType ;
    protected:
      typedef typename remove_const< GridFamily >::type::Traits Traits;

    protected:
      // type of the host grid
      typedef typename Traits::HostGridPartType  HostGridPartType;

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
      explicit IdEntity ( ExtraData data )
      : BaseType( data )
      {}

      IdEntity ( ExtraData data, const HostEntityType &hostEntity )
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
      typedef typename BaseType::HostGridPartType HostGridPartType;

      // type of extra data, e.g. a pointer to grid (here empty)
      typedef typename BaseType::ExtraData ExtraData;

    public:
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
      //! type of corresponding entity pointer
      typedef typename Traits::template Codim< codimension >::EntityPointer EntityPointer;

      /** \} */

      /** \name Construction, Initialization and Destruction
       *  \{ */

      /** \brief construct a null entity
       *  \param[in]  data  data pointer (here empty)
       */
      explicit IdEntity ( ExtraData data )
      : BaseType( data )
      {}

      /** \brief construct an initialized entity
       *
       *  \param[in]  data        data pointer (here empty)
       *  \param[in]  hostEntity  corresponding entity in the host grid
       *
       *  \note The reference to the host entity must remain valid as long as
       *        this entity is in use.
       */
      IdEntity ( ExtraData data, const HostEntityType &hostEntity )
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
      typename Traits::template Codim< codim >::EntityPointer
      subEntity ( int i ) const
      {
        typedef typename Traits::template Codim< codim >::EntityPointerImpl EntityPointerImpl;
        return EntityPointerImpl( data(), hostEntity().template subEntity< codim >( i ) );
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
