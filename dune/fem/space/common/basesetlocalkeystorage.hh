#ifndef DUNE_FEM_BASESETLOCALKEYSTORAGE_HH
#define DUNE_FEM_BASESETLOCALKEYSTORAGE_HH

//- Dune includes 
#if HAVE_DUNE_GEOMETRY
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>
#else
#include <dune/common/geometrytype.hh>
#include <dune/geometry/geometrytypeindex.hh>
#endif
#include <dune/common/exceptions.hh>

//- Fem includes 
#include <dune/fem/storage/singletonlist.hh>

namespace Dune 
{

namespace Fem {

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
      
      FactoryIF* clone() const { return new FactoryImpl<SingletonProvider> (); }
    };

    // pointer to apropriate factory 
    const FactoryIF* factory_; 
    // vector caching singleton pointers 
    std::vector< Entry* > entryStorage_;
  public:
    // export value type confomring to std::vector 
    typedef Entry value_type ;

    BaseSetLocalKeyStorage()
      : factory_( 0 ) 
      , entryStorage_()
    {} 

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
          entryStorage_[ i ] = factory_->getObject( otherEntry->geometryType() );
        }
      }
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
      //assert( typeid( factory_ ) == typeid( FactoryImpl< SingletonProvider >* ) );

      // get geometry type index 
      const size_t geomIndex = GlobalGeometryTypeIndex :: index( geomType ) ;

      if( entryStorage_.size() <= geomIndex ) 
        entryStorage_.resize( geomIndex + 1, (Entry* ) 0 );

      // if entry is still not used, insert it  
      if( entryStorage_[ geomIndex ] == 0 )
      {
        entryStorage_[ geomIndex ] = factory_->getObject( geomType );
        return true ;
      }
      return false ;
    }

    //! access to stored entry with given geometry type 
    const Entry& operator [] ( const GeometryType geomType ) const 
    {
      assert( GlobalGeometryTypeIndex :: index( geomType ) < entryStorage_.size() );
      assert( entryStorage_[ GlobalGeometryTypeIndex :: index( geomType ) ] != 0 );
      return *( entryStorage_[ GlobalGeometryTypeIndex :: index( geomType ) ]);
    }
  };

  } // end namespace Fem 

} // end namespace Dune 

#endif // DUNE_FEM_BASESETLOCALKEYSTORAGE_HH
