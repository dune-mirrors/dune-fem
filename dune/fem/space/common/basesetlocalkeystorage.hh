#ifndef DUNE_FEM_BASESETLOCALKEYSTORAGE_HH
#define DUNE_FEM_BASESETLOCALKEYSTORAGE_HH

#include <utility>

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/fem/common/forloop.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/storage/singletonlist.hh>


namespace Dune
{

  namespace Fem
  {

    /*! \brief storage class for base function set pointer
          and compiled local key pointers */
    template< class Entry >
    class BaseSetLocalKeyStorage
    {
      // interface class for factory
      class FactoryIF
      {
      protected:
        FactoryIF () {}
      public:
        virtual ~FactoryIF () {}
        virtual Entry* getObject( const GeometryType& geomType ) const = 0;
        virtual void removeObjects( std::vector<Entry*>& entryStorage ) const = 0;
        virtual FactoryIF* clone() const = 0;
      };

      // factory implementation depends on type of singleton provider
      template <class SingletonProvider>
      class FactoryImpl : public FactoryIF
      {
      public:
        FactoryImpl() {}

        Entry* getObject( const GeometryType& geomType ) const
        {
          return & SingletonProvider :: getObject( geomType );
        }

        void removeObjects( std::vector<Entry*>& entryStorage ) const
        {
          const size_t size = entryStorage.size();
          for( size_t i=0; i<size; ++i )
          {
            if( entryStorage[ i ] )
            {
              SingletonProvider :: removeObject( *(entryStorage[ i ]) );
              entryStorage[ i ] = 0 ;
            }
          }
        }

        virtual FactoryIF* clone() const { return new FactoryImpl<SingletonProvider> (); }
      };

      // pointer to apropriate factory
      const FactoryIF* factory_;
      // vector caching singleton pointers
      std::vector< Entry* > entryStorage_;
    public:
      // export value type confomring to std::vector
      typedef Entry value_type ;

      // default constructor
      BaseSetLocalKeyStorage()
        : factory_( 0 )
        , entryStorage_()
      {}

      //! copy constructor
      BaseSetLocalKeyStorage( const BaseSetLocalKeyStorage& other )
        : factory_( other.factory_ ? other.factory_->clone() : 0 )
        , entryStorage_( other.entryStorage_.size(), ( Entry * ) 0 )
      {
        // make a copy of the vector
        const size_t size = entryStorage_.size();
        for( size_t i=0; i<size; ++i )
        {
          Entry* otherEntry = other.entryStorage_[ i ];
          if( otherEntry )
          {
            // we need the interface method geometry
            // (on base function sets and compiled local keys )
            entryStorage_[ i ] = factory_->getObject( otherEntry->type() );
          }
        }
      }

      //! move constructor
      BaseSetLocalKeyStorage ( BaseSetLocalKeyStorage &&other )
        : factory_( other.factory_ ),
          entryStorage_( std::move( other.entryStorage_ ) )
      {
        other.factory_ = nullptr;
        other.entryStorage_.clear();
      }

      //! destructor
      ~BaseSetLocalKeyStorage()
      {
        if( entryStorage_.size() > 0 )
        {
          factory_->removeObjects( entryStorage_ );
        }

        delete factory_;
        factory_ = 0;
      }

      // get maxSize of compiled keys
      unsigned int maxSize() const
      {
        unsigned int maxSize = 0;
        const size_t size = entryStorage_.size() ;
        for( size_t i=0; i<size; ++i)
        {
          if( entryStorage_[ i ] )
          {
            unsigned int enSize = entryStorage_[ i ]->size();
            maxSize = std::max( enSize , maxSize );
          }
        }
        return maxSize;
      }

      //! insert entry to storage for given geometry type
      template <class SingletonProvider>
      bool insert( const GeometryType geomType )
      {
        // create factory if not existing yet
        if( factory_ == 0 )
        {
          factory_ = new FactoryImpl< SingletonProvider > ();
        }

        // check that type of factory is correct
        assert( dynamic_cast< const FactoryImpl< SingletonProvider >* > ( factory_ ) != 0 );

        // get geometry type index
        const size_t geomIndex = index( geomType ) ;

        if( entryStorage_.size() <= geomIndex )
          entryStorage_.resize( geomIndex + 1, (Entry* ) 0 );

        assert( geomIndex < entryStorage_.size() );

        // if entry is still not used, insert it
        if( entryStorage_[ geomIndex ] == 0 )
        {
          entryStorage_[ geomIndex ] = factory_->getObject( geomType );
          return true ;
        }
        return false ;
      }

      //! return true if an entry for this geometry type exists
      bool exists( const GeometryType& geomType ) const
      {
        if( index( geomType ) < static_cast< int >( entryStorage_.size() ) )
          return (entryStorage_[ index( geomType ) ] != 0) ;
        else
          return false;
      }

      //! access to stored entry with given geometry type
      const Entry& operator [] ( const GeometryType& geomType ) const
      {
        // assert( factory_ );
        assert( index( geomType ) < static_cast< int >( entryStorage_.size() ) );
        assert( entryStorage_[ index( geomType ) ] != 0 );
        return *( entryStorage_[ index( geomType ) ]);
      }
    protected:
      int index( const GeometryType& geomType ) const
      {
        return LocalGeometryTypeIndex::index( geomType );
      }
    };



    /** \brief class for storage local keys for a given range of polynomial order and
               available geometry type */
    template< class CompiledLocalKey, unsigned int minPolOrder, unsigned int maxPolOrder >
    class CompiledLocalKeyContainer
    {
    public:
      // type of compiled local key
      typedef CompiledLocalKey CompiledLocalKeyType;

      //! type of storage class for compiled local keys
      typedef BaseSetLocalKeyStorage< CompiledLocalKeyType > LocalKeyStorageType;

      // vector containing storages for each polynomial order
      typedef std::vector< LocalKeyStorageType > LocalKeyVectorType;

    protected:
      enum { numOrders = maxPolOrder - minPolOrder + 1 };

      template <int pOrd>
      struct ConstructCompiledLocalKeys
      {
        /** HelperClasses
           \brief
           CompiledLocalKeyFactory method createObject and
           deleteObject for the SingletonList
        */
        class CompiledLocalKeyFactory
        {
        public:
          //! create new BaseFunctionSet
          static CompiledLocalKeyType* createObject( const GeometryType& type )
          {
            return new CompiledLocalKeyType( type, pOrd );
          }

          //! delete BaseFunctionSet
          static void deleteObject( CompiledLocalKeyType* obj )
          {
            delete obj;
          }
        };

        static void apply( LocalKeyVectorType& compiledLocalKeys,
                           const GeometryType& geometryType )
        {
          const size_t k = pOrd ;

          //! type of singleton list (singleton provider) for compiled local keys
          typedef SingletonList
            < GeometryType, CompiledLocalKeyType, CompiledLocalKeyFactory >
            CompiledLocalKeySingletonProviderType;

          // insert compiled local key
          compiledLocalKeys[ k - minPolOrder ].template insert< CompiledLocalKeySingletonProviderType > ( geometryType );
        }
      };

    protected:
      // all lagrange point sets for available geometry types
      LocalKeyVectorType compiledLocalKeys_ ;

    private:
      CompiledLocalKeyContainer( const CompiledLocalKeyContainer& );

    public:
      template <class GridPart>
      CompiledLocalKeyContainer( const GridPart& gridPart )
        : compiledLocalKeys_( numOrders )
      {
        typedef typename GridPart :: IndexSetType IndexSetType ;
        typedef typename GridPart :: GridType  GridType ;
        const IndexSetType &indexSet = gridPart.indexSet();

        // get all available geometry types
        AllGeomTypes< IndexSetType, GridType > allGeometryTypes( indexSet );
        const std :: vector< GeometryType >& geometryTypes
          = allGeometryTypes.geomTypes( 0 );

        // create compiled local keys
        for( unsigned int i = 0; i < geometryTypes.size(); ++i )
        {
          Fem::ForLoop< ConstructCompiledLocalKeys, minPolOrder, maxPolOrder > ::
              apply( compiledLocalKeys_, geometryTypes[ i ] );
        }
      }

      /** \brief provide access to all compiled local keys for a given polynomial order
       *
       *  \param[in]  order polynomial order for given geometry type
       *
       *  \returns CompiledLocalKeys storage
       */
      inline const LocalKeyStorageType& compiledLocalKeys( const int order ) const
      {
        assert( order - minPolOrder >= 0 );
        assert( int( order - minPolOrder ) < int( compiledLocalKeys_.size() ) );
        return compiledLocalKeys_[ order - minPolOrder ];
      }

      /** \brief provide access to the compiled local keys for a geometry type and polynomial order
       *
       *  \param[in]  type  type of geometry the compiled local key is requested for
       *  \param[in]  order polynomial order for given geometry type
       *
       *  \returns CompiledLocalKey
       */
      inline const CompiledLocalKeyType &compiledLocalKey( const GeometryType& type, const int order ) const
      {
        return compiledLocalKeys( order )[ type ];
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // DUNE_FEM_BASESETLOCALKEYSTORAGE_HH
