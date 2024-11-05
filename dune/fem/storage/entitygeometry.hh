#ifndef DUNE_FEM_STORAGE_ENTITYGEOMETRY_HH
#define DUNE_FEM_STORAGE_ENTITYGEOMETRY_HH

// C++ includes
#include <cassert>
#include <cstddef>
#include <memory>
#include <optional>

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{

  namespace Fem
  {
    /**
     * \brief implementation of entity and geometry storage for
     *        basis function set and local functions
     *
     * \tparam  Entity            entity type
     * \tparam  storeGeometry     if true, geometry is also stored
     *
     * \note EntityGeometryStorageImpl must be a copyable object.
     */
    template< class Entity, bool storeGeometry = true >
    class EntityGeometryStorageImpl
    {

    public:
      //! \brief entity type
      typedef Entity EntityType;

      //! \brief type of geometry
      typedef typename EntityType::Geometry Geometry ;

      //! \brief type of coordinate field
      typedef typename Geometry::ctype ctype;

      //! \brief type of reference element
      typedef std::decay_t< decltype( Dune::ReferenceElements< ctype, Geometry::coorddimension >::
          general( std::declval< const Dune::GeometryType & >() ) ) > ReferenceElementType;

    protected:
      typedef std::optional< EntityType > EntityStorageType;

      struct Empty
      {
        operator bool () const { return false; }
      };
      typedef typename std::conditional< storeGeometry,
                                         std::optional< Geometry >,
                                         Empty > :: type GeometryStorageType;

    public:
      //! \brief constructor
      EntityGeometryStorageImpl ()
#ifdef TESTTHREADING
        : thread_(-1)
#endif
      {}

      //! \brief constructor
      explicit EntityGeometryStorageImpl( const EntityType &entity )
#ifdef TESTTHREADING
        : thread_( MPIManager::thread() )
#endif
      {
        bind( entity );
      }

      //! \brief copy constructor
      EntityGeometryStorageImpl ( const EntityGeometryStorageImpl &other )
        : entity_( other.entity_ )
#ifdef TESTTHREADING
        , thread_( MPIManager::thread() )
#endif
      {
        copyGeometry( other );
      }

      //! \brief assignment operator
      EntityGeometryStorageImpl &operator= ( const EntityGeometryStorageImpl &other )
      {
        entity_ = other.entity_;
        copyGeometry( other );
#ifdef TESTTHREADING
        thread_ = MPIManager::thread();
#endif
        return *this;
      }

      //! \brief return entity
      const Entity &entity () const
      {
        assert( valid() );
        return entity_.value();
      }

      //! \brief return true if entity pointer is set
      bool valid () const { return bool(entity_); }

      //! \brief return geometry
      const Geometry& geometry () const
      {
        if constexpr ( storeGeometry )
        {
          assert( geometry_ );
          return geometry_.value();
        }
        else
        {
          DUNE_THROW(InvalidStateException,"EntityGeometryStorageImpl::geometry not available when storeGeometry is false!");
          return *((Geometry *) nullptr);
        }
      }

      //! \brief return geometry type
      Dune::GeometryType type () const { return entity().type(); }

      //! \brief return reference element
      const ReferenceElementType& referenceElement () const
      {
        return Dune::ReferenceElements< ctype, Geometry::coorddimension >::general( type() );
      }

      //! set new entity object and geometry if enabled
      void bind( const EntityType& entity )
      {
#ifdef TESTTHREADING
        if (thread_==-1) thread_ = MPIManager::thread();
        if (thread_ != MPIManager::thread())
        {
          std::cout << "EntityGeometryStorageImpl::bind: wrong thread number!" std::endl;
          assert(0);
          std::abort();
        }
        if (entity_ || (storeGeometry && geometry_) )
        {
          std::cout << "EntityGeometryStorageImpl::bind: bind called on object before unbind was called" << std::endl;
          std::abort();
        }
        assert(!entity_ && !geometry_); // this will fail with dune-fem-dg
#endif

        entity_.emplace( entity );
        if constexpr ( storeGeometry )
        {
          // Note that this should be geometry_ = entity.geometry()
          // But Dune::Geometries are not assignable ...
          // geometry_.reset();
          geometry_.emplace( entity.geometry() );
        }
      }

      //! release entity and geometry object
      void unbind()
      {
#ifdef TESTTHREADING
        if (thread_ != MPIManager::thread())
        {
          std::cout << "EntityGeometryStorageImpl::unbind: wrong thread number" << std::endl;
          assert(0);
          std::abort();
        }
#endif
        entity_.reset();
        if constexpr ( storeGeometry )
        {
          geometry_.reset();
        }
      }

    protected:
      void copyGeometry( const EntityGeometryStorageImpl& other )
      {
        if constexpr ( storeGeometry )
        {
          // Note that this should be geometry_ = entity.geometry()
          // But Dune::Geometries are not assignable ...
          geometry_.reset();
          if( other.geometry_ )
            geometry_.emplace( other.geometry_.value() );
        }
      }

    protected:
      EntityStorageType    entity_;
      GeometryStorageType  geometry_;
#ifdef TESTTHREADING
      int thread_;
#endif
    };

    template <class Entity>
    using EntityGeometryStorage = EntityGeometryStorageImpl< Entity, true >;

    template <class Entity>
    using EntityStorage = EntityGeometryStorageImpl< Entity, false >;
  } // end namespace Fem

} // end namespace Dune
#endif
