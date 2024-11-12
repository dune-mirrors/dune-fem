#ifndef DUNE_FEM_QUADPROVIDER_HH
#define DUNE_FEM_QUADPROVIDER_HH

#include <iostream>
#include <memory>
#include <map>
#include <vector>

#include <dune/common/shared_ptr.hh>
#include <dune/fem/quadrature/quadratureimp.hh>
#include <dune/fem/quadrature/idprovider.hh>
#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{

  namespace Fem
  {

    /*! \class FemQuadratureKey
     *  \ingroup Quadrature
     *  \brief A simple quadrature key class for use FemPy
     *
     */
    class FemQuadratureKey
    {
      int id_;

      // something like maxOrder
      static const int maxFirst  = 25 ;
      // something like different quadratures of the same order
      static const int maxSecond = 10 ;
    public:
      //! empty constructor
      FemQuadratureKey() : id_( -1 ) {}

      //! copy constructor
      FemQuadratureKey( const FemQuadratureKey& key ) = default;

      // this may need to be adjusted in case more than 100 quadratures
      // need to be stored
      static const int highest_order = maxFirst * maxSecond ;

      //! constructor taking to ids, like std::pair
      FemQuadratureKey( const int first, const int second )
        : id_( second * maxFirst + first )
      {
        assert( first  < maxFirst );
        assert( second < maxSecond );
      }

      //! constructor taking only order (fallback for standard Fem quadratures)
      FemQuadratureKey( const int first )
        : id_( first )
      {
        assert( first  < maxFirst );
      }

      //! cast into fast int identifier based storage of quadratures
      operator int () const { return id_; }

      //! return first component
      int first () const  { return id_ % highest_order ; }
      //! return second component
      int second () const { return id_ / highest_order ; }
    };


    /*! \class QuadCreator
     *  \ingroup Quadrature
     *  \brief the actual quadrature storage
     *
     *  QuadCreator is a utility class providing the actual quadrature storage.
     *
     *  The template argument is used to distinguish classes for different geometry
     *  types (maybe GeometryType :: BasicType would be a better choice).
     */
    template< unsigned int dummy >
    class QuadCreator
    {
      //! class holding map with pointer to quadrature objects
      template< class QuadImp, class QuadratureKey, unsigned int >
      class QuadratureStorage;
    private:
      //! class holding map with pointer to quadrature objects
      template< class QuadImp, class QuadratureKey >
      class QuadratureStorage< QuadImp, QuadratureKey, dummy >
      {
      public:
        typedef QuadImp QuadType;

      protected:
        struct Storage
        {
          typedef std :: map< QuadratureKey, std::unique_ptr< QuadType > >  StorageType;
          std::mutex mutex_;
          StorageType storage_;

          static QuadImp* create( const GeometryType &geometry, const QuadratureKey& key )
          {
            return instance().createImpl( geometry, key );
          }

        private:
          static Storage& instance()
          {
            return Singleton< Storage >::instance();
          }

          QuadImp* createImpl( const GeometryType &geometry, const QuadratureKey& key )
          {
            std::lock_guard< std::mutex > guard( mutex_ );

            auto& quadPtr = storage_[ key ];
            if( ! quadPtr )
            {
              quadPtr.reset( new QuadImp( geometry, key, IdProvider :: instance().newId() ) );
            }

            assert( quadPtr );
            return quadPtr.operator->();
          }
        };

        typedef std :: map< QuadratureKey, std::unique_ptr< QuadType, Dune::null_deleter<QuadType> > > PointerStorageType;
        PointerStorageType quadPtrs_;

      public:
        QuadratureStorage () {}

        QuadType &getQuadrature( const GeometryType &geometry, const QuadratureKey& key )
        {
          auto& quadPtr = quadPtrs_[ key ];
          if( ! quadPtr )
          {
            // create is mutex protected
            quadPtr.reset( Storage::create( geometry, key ) );
          }

          assert( quadPtr );
          return *quadPtr;
        }
      };

      //! class holding vector with pointer to quadrature objects
      template< class QuadImp>
      class QuadratureStorage< QuadImp, int, dummy > // QuadratureKey == int
      {
      public:
        typedef QuadImp QuadType;

      protected:
        struct Storage
        {
          std::mutex mutex_;

          // vector storing the actual objects
          std::vector< std::unique_ptr< QuadType > > storage_;
          Storage() : storage_( QuadType :: maxOrder() + 1 )
          {}

          static QuadType* create( const GeometryType &geometry, int order )
          {
            return instance().createImpl( geometry, order );
          }

      private:
          static Storage& instance()
          {
            return Singleton< Storage >::instance();
          }

          QuadType* createImpl( const GeometryType &geometry, int order )
          {
            std::lock_guard< std::mutex > guard( mutex_ );

            auto& quadPtr = storage_[ order ];
            if( ! quadPtr )
            {
              quadPtr.reset( new QuadImp( geometry, order, IdProvider :: instance().newId() ) );
            }
            assert( quadPtr );
            return quadPtr.operator ->();
          }
        };

        typedef std::vector< QuadType* > QuadPointerVecType;
        // store all the pointers for all threads separately to avoid race conditions
        QuadPointerVecType quadPtrs_;

      public:
        QuadratureStorage ()
          : quadPtrs_( QuadPointerVecType( QuadType :: maxOrder() + 1, nullptr ))
        {
        }

      public:
        QuadImp &getQuadrature( const GeometryType &geometry, int order )
        {
          if(order >= int(quadPtrs_.size()) )
          {
#ifndef NDEBUG
            static thread_local bool showMessage = true ;
            if( showMessage )
            {
              std::cerr << "WARNING: QuadratureStorage::getQuadrature: A quadrature of order " << order
                        << " is not implemented!" << std::endl
                        << "Choosing maximum order: " << quadPtrs_.size()-1 << std::endl << std::endl;
              showMessage = false;
            }
#endif
            order = quadPtrs_.size() - 1;
          }

          auto& quadPtr = quadPtrs_[ order ];
          if( ! quadPtr )
          {
            // create is mutex protected
            quadPtr = Storage::create( geometry, order );
          }

          assert( quadPtr );
          return *( quadPtr );
        }
      }; // end class QuadratureStorage

      //! class holding vector with pointer to quadrature objects
      template< class QuadImp >
      class QuadratureStorage< QuadImp, FemQuadratureKey, dummy >
        : public QuadratureStorage< QuadImp, int, dummy >
      {
      };

    public:
      /*! \brief provide quadrature
       *
       *  \param[in]  geometry  type of geometry, the quadrature is requested for
       *  \param[in]  key       key to identify quadrature, i.e. minimal order of the requested quadrature
       */
      template< class QuadImp, class QuadratureKey >
      static const QuadImp &provideQuad( const GeometryType&  geometry,
                                         const QuadratureKey& key )
      {
        // QuadratureStorage stores only pointers to the quadratures that are
        // created and stored by a singleton storage
        static thread_local QuadratureStorage< QuadImp, QuadratureKey, dummy > storage;
        return storage.getQuadrature( geometry, key );
      }

      /*! \brief provide quadrature
       *
       *  \param[in]  geometry  type of geometry, the quadrature is requested for
       *  \param[in]  key       key to identify quadrature, i.e. minimal order of the requested quadrature
       *  \param[in]  defaultOrder  to identify polyhedral geometries
       */
      template< class QuadImp,  class QuadratureKey >
      static const QuadImp &provideQuad( const GeometryType&  geometry,
                                         const QuadratureKey& key,
                                         const int defaultOrder )
      {
        // this function should only be called for geometry types equal to none
        assert( geometry.isNone() );
        DUNE_THROW(NotImplemented,"provideQuad for polyhedral cells (defaultOrder = 0) not implemented for arbitrary QuadratureKey!");
        QuadImp* ptr = nullptr;
        return *ptr;
      }

      /*! \brief provide quadrature
       *
       *  \param[in]  geometry  type of geometry, the quadrature is requested for
       *  \param[in]  order     minimal order of the requested quadrature
       *  \param[in]  defaultOrder  to identify polyhedral geometries
       */
      template< class QuadImp >
      static const QuadImp &provideQuad( const GeometryType&  geometry,
                                         const int ,
                                         const int defaultOrder )
      {
        assert( geometry.isNone() );
        // QuadratureStorage stores only pointers to the quadratures that are
        // created and stored by a singleton storage
        static thread_local QuadratureStorage< QuadImp, int, dummy > storage;
        return storage.getQuadrature( geometry, defaultOrder );
      }
    };

    /*! \class QuadratureProvider
     *  \ingroup Quadrature
     *  \brief provide a single instance pool of quadratures
     *
     *  QuadratureProvider follows the monostate pattern. It provides a single
     *  point of access (and storage) for the actual implementation of
     *  quadratures. Hence, the expensive creations of quadratures should be
     *  reduced to a minimum.
     *
     *  There are the following specializations:
     *  - QuadratureProvider<FieldImp,0,QuadratureTraits>
     *  - QuadratureProvider<FieldImp,1,QuadratureTraits>
     *  - QuadratureProvider<FieldImp,2,QuadratureTraits>
     *  - QuadratureProvider<FieldImp,3,QuadratureTraits>
     */
    template< typename FieldImp, int dim, template< class, int > class IntegrationTraits >
    class QuadratureProvider
    {
    public:
      typedef FieldImp FieldType;

      static const int dimension = dim;

    private:
      typedef QuadratureProvider< FieldType, dimension, IntegrationTraits > ThisType;

      typedef IntegrationTraits< FieldType, dimension > QuadratureTraits;

      // to be removed
      typedef QuadratureTraits FactoryTraits;

      template <class PointQuadrature, class QuadratureKey>
      class PointQuadratureStorage : public PointQuadrature
      {
      public:
        // only call IdProvider ::instance().newId() when object is created
        PointQuadratureStorage( const GeometryType &geometry, const QuadratureKey& quadKey )
          : PointQuadrature( geometry, quadKey, IdProvider::instance().newId() )
        {}
      };

    public:
      //! type of integration point list implementation
      typedef typename QuadratureTraits::IntegrationPointListType  QuadratureImplementationType;

      //! Access to the quadrature implementations.
#if 0
      template <class FactoryTraits>
      static const QuadratureImplementationType &getQuadrature( const FactoryTraits traits,
                                                            const GeometryType &geometry,
                                                            const typename FactoryTraits::QuadratureKeyType& quadKey )
#else
      static const QuadratureImplementationType &getQuadrature( const GeometryType &geometry,
                                                                const typename FactoryTraits::QuadratureKeyType& quadKey )
#endif
      {
        static_assert( std::is_same< QuadratureImplementationType,
                                     typename FactoryTraits::IntegrationPointListType> :: value );

        // for 0d point quadratures
        if constexpr ( dimension == 0 )
        {
          typedef typename FactoryTraits :: PointQuadratureType PointQuadratureType;
          typedef typename FactoryTraits :: QuadratureKeyType  QuadratureKeyType;

          typedef PointQuadratureStorage< PointQuadratureType, QuadratureKeyType > PointQuadratureStorageType;

          assert( geometry.isCube() || geometry.isSimplex() );
          return Singleton< PointQuadratureStorageType > :: instance( geometry, quadKey );
        }
        else // all other cases
        {
          // for 1d return LineQuadrature
          if constexpr ( dimension == 1 )
          {
            typedef typename FactoryTraits::LineQuadratureType  LineQuadratureType;
            assert( geometry.isCube() || geometry.isSimplex() );
            return QuadCreator< 0 > :: template provideQuad< LineQuadratureType > ( geometry, quadKey );
          }

          // for 2d Simplex and Cube
          if constexpr ( dimension >= 2 )
          {
            typedef typename FactoryTraits::SimplexQuadratureType SimplexQuadratureType;
            typedef typename FactoryTraits::CubeQuadratureType    CubeQuadratureType;
            if( geometry.isSimplex() )
            {
              return QuadCreator< 0 > ::
                template provideQuad< SimplexQuadratureType > ( geometry, quadKey );
            }

            if( geometry.isCube() )
            {
              return QuadCreator< 1 > ::
                template provideQuad< CubeQuadratureType >    ( geometry, quadKey ) ;
            }

            if( geometry.isNone() )
            {
              // dummy return for polygonal grid cells, i.e. geometry type none
              return QuadCreator< 1 > :: template provideQuad< CubeQuadratureType > ( geometry, 0 );
            }

            if constexpr ( dimension == 3 )
            {
              typedef typename FactoryTraits::PrismQuadratureType    PrismQuadratureType;
              typedef typename FactoryTraits::PyramidQuadratureType  PyramidQuadratureType;

              if( geometry.isPrism() )
                return QuadCreator< 2 > :: template provideQuad< PrismQuadratureType >
                  ( geometry, quadKey );
              if( geometry.isPyramid() )
                return QuadCreator< 3 > :: template provideQuad< PyramidQuadratureType >
                  ( geometry, quadKey );
            }
          }

          DUNE_THROW( RangeError, "Element type not available for dimension " << dimension );
          // dummy return
          if constexpr ( dimension == 1 )
          {
            typedef typename FactoryTraits::LineQuadratureType  LineQuadratureType;
            return QuadCreator< 0 > :: template provideQuad< LineQuadratureType > ( geometry, quadKey, 0 );
          }
          else
          {
            typedef typename FactoryTraits::SimplexQuadratureType  SimplexQuadratureType;
            return QuadCreator< 0 > :: template provideQuad< SimplexQuadratureType > ( geometry, quadKey, 0 );
          }
        } // else not dimension 0
      }

#if 0
      //! Access to the quadrature implementations.
      template <class FactoryTraits>
      static const QuadratureImplementationType &getQuadrature( const FactoryTraits traits,
                                                            const GeometryType &geometry,
                                                            const GeometryType &elementGeometry,
                                                            const typename FactoryTraits::QuadratureKeyType& quadKey )
#else
      //! Access to the quadrature implementations.
      static const QuadratureImplementationType &getQuadrature( const GeometryType &geometry,
                                                                const GeometryType &elementGeometry,
                                                                const typename FactoryTraits::QuadratureKeyType& quadKey )
#endif
      {
        if constexpr ( dimension == 0 )
        {
          return getQuadrature(geometry, quadKey);
        }
        else
        {
          if constexpr ( dimension == 1 )
          {
            typedef typename FactoryTraits::LineQuadratureType  LineQuadratureType;
            assert( geometry.isCube() || geometry.isSimplex() );
            // we need here to distinguish between the basic types
            // otherwise the this won't work for UGGrid
            return ( elementGeometry.isSimplex() ) ?
              QuadCreator< 0 > :: template provideQuad< LineQuadratureType > ( geometry, quadKey ) :
              QuadCreator< 1 > :: template provideQuad< LineQuadratureType > ( geometry, quadKey ) ;
          }

          if constexpr ( dimension == 2 )
          {
            typedef typename FactoryTraits::SimplexQuadratureType  SimplexQuadratureType;
            typedef typename FactoryTraits::CubeQuadratureType     CubeQuadratureType;
            assert( geometry.isCube() || geometry.isSimplex() );

            // if geometry is simplex return simplex quadrature
            if ( geometry.isSimplex() )
            {
              // check element geometry to provide quadratures with different ids
              if( elementGeometry.isSimplex() )
                return QuadCreator< 0 > :: template provideQuad< SimplexQuadratureType > ( geometry, quadKey ) ;
              else if( elementGeometry.isCube() )
                return QuadCreator< 1 > :: template provideQuad< SimplexQuadratureType > ( geometry, quadKey ) ;
              else if( elementGeometry.isPrism() )
                return QuadCreator< 2 > :: template provideQuad< SimplexQuadratureType > ( geometry, quadKey ) ;
              else if( elementGeometry.isPyramid() )
                return QuadCreator< 3 > :: template provideQuad< SimplexQuadratureType > ( geometry, quadKey ) ;
              else
                DUNE_THROW( RangeError, "Element type not available for dimension 3" );
            }
            else
            {
              // return cube quadrature
              // check element geometry to provide quadratures with different ids
              if( elementGeometry.isSimplex() )
                return QuadCreator< 4 > :: template provideQuad< CubeQuadratureType > ( geometry, quadKey ) ;
              else if( elementGeometry.isCube() )
                return QuadCreator< 5 > :: template provideQuad< CubeQuadratureType > ( geometry, quadKey ) ;
              else if( elementGeometry.isPrism() )
                return QuadCreator< 6 > :: template provideQuad< CubeQuadratureType > ( geometry, quadKey ) ;
              else if( elementGeometry.isPyramid() )
                return QuadCreator< 7 > :: template provideQuad< CubeQuadratureType > ( geometry, quadKey ) ;
              else
                DUNE_THROW( RangeError, "Element type not available for dimension 3" );
            }
          }

          if constexpr ( dimension >= 3 )
          {
            typedef typename FactoryTraits::SimplexQuadratureType  SimplexQuadratureType;
            DUNE_THROW( RangeError, "QuadProvider::getQuadrature not implemented for 3d face quadratures!" );
            // dummy return
            return QuadCreator< 0 > :: template provideQuad< SimplexQuadratureType >
              ( geometry, quadKey, 0 );
          }
        }
      }

      QuadratureProvider() = delete;
      QuadratureProvider( const ThisType& ) = delete;
      QuadratureProvider &operator=( const ThisType& ) = delete;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_QUADPROVIDER_HH
