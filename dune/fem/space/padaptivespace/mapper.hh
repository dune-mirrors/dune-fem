#ifndef DUNE_PADAPTIVELAGRANGESPACE_MAPPER_HH
#define DUNE_PADAPTIVELAGRANGESPACE_MAPPER_HH

//- Dune includes 
#include <dune/common/geometrytype.hh>
#include <dune/common/exceptions.hh>

//- Dune-Fem includes 
#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/misc/codimmap.hh>
#include <dune/fem/misc/metaprogramming.hh>
#include <dune/fem/space/common/dofmanager.hh>

#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>

//- local includes 
#include <dune/fem/space/lagrangespace/lagrangepoints.hh>

#include <dune/grid/utility/persistentcontainer.hh>

// include generic adaptive dof mapper 
#include <dune/fem/space/mapper/genericadaptivedofmapper.hh>

namespace Dune
{

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
      virtual void removeObject( Entry& entry ) const = 0;
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

      void removeObject( Entry& entry ) const 
      {
        SingletonProvider :: removeObject( entry );
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
      remove();
    }

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

  protected:  
    void remove() 
    {
      const size_t size = entryStorage_.size();
      if( size == 0 ) return ;

      assert( factory_ );
      for(size_t i=0; i<size; ++i)
      {
        if( entryStorage_[ i ] )
          factory_->removeObject( *(entryStorage_[ i ]) );
      } 
      delete factory_; 
      factory_ = 0 ;
    }
  };


  
  template< class GridPart, int polOrder >
  class PAdaptiveLagrangeMapper;

  template< class GridPart, int polOrder >
  struct PAdaptiveLagrangeMapperTraits
  {
    typedef GridPart GridPartType;
    
    static const int polynomialOrder = polOrder;

    typedef typename GridPartType::template Codim< 0 >::IteratorType::Entity EntityType;
    typedef PAdaptiveLagrangeMapper< GridPartType, polynomialOrder > DofMapperType;
    typedef DefaultDofMapIterator< EntityType, DofMapperType > DofMapIteratorType;

    //! type of the compiled local key 
    typedef LagrangePointSet< GridPartType, polynomialOrder >  CompiledLocalKeyType;
    typedef BaseSetLocalKeyStorage< CompiledLocalKeyType > BaseSetLocalKeyStorageType;

    typedef std::vector< BaseSetLocalKeyStorageType > CompiledLocalKeyVectorType ;
  };



  // First Order Lagrange Mapper
  // ---------------------------

  template< class GridPart >
  class PAdaptiveLagrangeMapper< GridPart, 1 >
  : public CodimensionMapper< GridPart, GridPart::GridType::dimension >
  {
    typedef PAdaptiveLagrangeMapper< GridPart, 1 > ThisType;
    typedef CodimensionMapper< GridPart, GridPart::GridType::dimension > BaseType;

  public:
    //! type of the grid part
    typedef typename BaseType::GridPartType GridPartType;

    //! type of entities (codim 0)
    typedef typename BaseType::EntityType EntityType;

    //! type of the underlying grid
    typedef typename GridPartType::GridType GridType;

    //! type of coordinates within the grid
    typedef typename GridType::ctype FieldType;

    //! dimension of the grid
    static const int dimension = GridType::dimension;

    //! order of the Lagrange polynoms
    static const int polynomialOrder = 1; 

    // my traits class 
    typedef PAdaptiveLagrangeMapperTraits< GridPart, polynomialOrder > Traits;

    typedef typename Traits :: CompiledLocalKeyVectorType  CompiledLocalKeyVectorType;

  public:
    //! constructor
    PAdaptiveLagrangeMapper ( const GridPartType &gridPart, CompiledLocalKeyVectorType& )
    : BaseType( gridPart )
    {}

    bool fixedDataSize ( const int codim ) const
    {
      return true;
    }

    int polynomOrder( const EntityType& entity ) const 
    {
      return 1;
    }

    void setPolynomOrder( const EntityType& entity, const int polOrd ) 
    {
    }
  };



  // Higher Order Lagrange Mapper
  // ----------------------------

  template< class GridPart, int polOrder >
  class PAdaptiveLagrangeMapper 
    : public GenericAdaptiveDofMapper< PAdaptiveLagrangeMapperTraits< GridPart, polOrder > >
  {
  public:   
    // my traits class 
    typedef PAdaptiveLagrangeMapperTraits< GridPart, polOrder > Traits;

  private:   
    typedef PAdaptiveLagrangeMapper< GridPart, polOrder > ThisType;
    typedef GenericAdaptiveDofMapper< Traits > BaseType;

  public:
    //! type of the grid part
    typedef typename Traits::GridPartType GridPartType;

    //! type of compiled local keys vector 
    typedef typename Traits :: CompiledLocalKeyVectorType  CompiledLocalKeyVectorType;

  public:
    //! constructor
    PAdaptiveLagrangeMapper ( const GridPartType &gridPart,
                              CompiledLocalKeyVectorType &compiledLocalKeys )
      : BaseType( gridPart, compiledLocalKeys )
    {
    }

    //! sort of copy constructor
    PAdaptiveLagrangeMapper ( const PAdaptiveLagrangeMapper& other,
                              CompiledLocalKeyVectorType &compiledLocalKeys )
      : BaseType( other, compiledLocalKeys ) 
    {} 
  };

} // end namespace Dune 

#endif // #ifndef DUNE_LAGRANGESPACE_MAPPER_HH
