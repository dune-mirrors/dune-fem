#ifndef DUNE_LAGRANGESPACE_LAGRANGEPOINTS_HH
#define DUNE_LAGRANGESPACE_LAGRANGEPOINTS_HH

#include <dune/grid/common/referenceelements.hh>

#include "genericgeometry.hh"
#include "genericlagrangepoints.hh"

namespace Dune
{

  template< GeometryType :: BasicType type,
            unsigned int dim,
            unsigned int order >
  class LagrangePoint
  : public GenericLagrangePoint
    < typename GeometryWrapper< type, dim > :: GenericGeometryType, order >
  {
  private:
    typedef GeometryWrapper< type, dim > GeometryWrapperType;
    typedef GenericLagrangePoint
      < typename GeometryWrapperType :: GenericGeometryType, order >
      BaseType;
    typedef LagrangePoint< type, dim, order > ThisType;

  public:
    enum { dimension = BaseType :: dimension };
    typedef typename BaseType :: DofCoordinateType DofCoordinateType;

    enum { polynomialOrder = BaseType :: polynomialOrder };

    enum { numLagrangePoints = BaseType :: numLagrangePoints };

  public:
    template< unsigned int codim >
    class Codim
    {
    public:
      static inline unsigned int maxDofs ()
      {
        return BaseType :: template Codim< codim > :: maxDofs();
      }
    };
    
  public:
    inline LagrangePoint ( unsigned int index )
    : BaseType( index )
    {
    }

    inline LagrangePoint ( const BaseType &point )
    : BaseType( point )
    {
    }

    inline void dofSubEntity ( unsigned int &codim,
                               unsigned int &subEntity )
    {
      BaseType :: dofSubEntity( codim, subEntity );
      GeometryWrapperType :: duneSubEntity( codim, subEntity );
    }

    inline void dofSubEntity ( unsigned int &codim,
                               unsigned int &subEntity,
                               unsigned int &dofNumber )
    {
      BaseType :: dofSubEntity( codim, subEntity, dofNumber );
      GeometryWrapperType :: duneSubEntity( codim, subEntity );
    }

    static inline unsigned int entityDofNumber ( unsigned int codim,
                                                 unsigned int subEntity,
                                                 unsigned int dof )
    {
      GeometryWrapperType :: duneSubEntity( codim, subEntity );
      return BaseType :: entityDofNumber( codim, subEntity, dof );
    }
  };


  
  template< class FieldImp, unsigned int dim, unsigned int polOrder >
  class LagrangePointSetInterface;

  

  template< class FieldImp,
            unsigned int dim,
            unsigned int codim,
            unsigned int polOrder >
  class SubEntityLagrangePointIterator
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = dim };
    enum { codimension = codim };

    enum { polynomialOrder = polOrder };

    typedef FieldVector< FieldType, dimension > pointType;

  private:
    typedef SubEntityLagrangePointIterator< FieldType,
                                            dimension,
                                            codimension,
                                            polynomialOrder >
      ThisType;

    typedef LagrangePointSetInterface< FieldType, dimension, polynomialOrder >
      LagrangePointSetType;

    typedef ReferenceElementContainer< FieldType, dimension >
      ReferenceElementContainerType;
    typedef typename ReferenceElementContainerType :: value_type
      ReferenceElementType;

  private:
    ReferenceElementContainerType refElementContainer_;
    
    const LagrangePointSetType *lagrangePointSet_;
    unsigned int subEntity_;

    unsigned int codim_;
    unsigned int subIndex_, numSubIndices_;
    unsigned int subSubEntity_;
    unsigned int dofNumber_, numDofs_;

  private:
    inline SubEntityLagrangePointIterator
      ( const LagrangePointSetType &lagrangePointSet,
        const unsigned int subEntity,
        const bool beginIterator )
    : refElementContainer_()
    {
      lagrangePointSet_ = &lagrangePointSet;
      subEntity_ = subEntity;

      subIndex_ = 0;
      numSubIndices_ = 1;
      subSubEntity_ = subEntity_;
        
      dofNumber_ = 0;
      numDofs_ = lagrangePointSet_->numDofs( codimension, subSubEntity_ );
        
      if( beginIterator ) {
        codim_ = codimension;
        assertDof();
      } else
        codim_ = dimension+1;
    }
  
  public:
    inline SubEntityLagrangePointIterator ()
    : refElementContainer_()
    {
      lagrangePointSet_ = NULL;
      subEntity_ = 0;

      codim_ = dimension + 1;
      subIndex_ = 0;
      dofNumber_ = 0;
    }
    
    inline SubEntityLagrangePointIterator ( const ThisType &other )
    : refElementContainer_()
    {
      lagrangePointSet_ = other.lagrangePointSet_;
      subEntity_ = other.subEntity_;

      codim_ = other.codim_;
      
      subIndex_ = other.subIndex_;
      numSubIndices_ = other.numSubIndices_;
      subSubEntity_ = other.subSubEntity_;
      
      dofNumber_ = other.dofNumber_;
      numDofs_ = other.numDofs_;
    }

    inline ThisType& operator= ( const ThisType &other )
    {
      lagrangePointSet_ = other.lagrangePointSet_;
      subEntity_ = other.subEntity_;

      codim_ = other.codim_;
      
      subIndex_ = other.subIndex_;
      numSubIndices_ = other.numSubIndices_;
      subSubEntity_ = other.subSubEntity_;
      
      dofNumber_ = other.dofNumber_;
      numDofs_ = other.numDofs_;
    }

    inline unsigned int operator* () const
    {
      assert( lagrangePointSet_ != NULL );
      assert( codim_ <= dimension );

      return lagrangePointSet_->entityDofNumber
               ( codim_, subSubEntity_, dofNumber_ );
    }

    inline ThisType& operator++ ()
    {
      if( codim_ <= dimension ) {
        ++dofNumber_;
        assertDof();
      }
      return *this;
    }

    inline bool operator== ( const ThisType& other ) const
    {
      if( (other.codim_ != codim_)
          || (other.subIndex_ != subIndex_)
          || (other.dofNumber_ != dofNumber_) )
        return false;

      return (other.lagrangePointSet_ == lagrangePointSet_)
             && (other.subEntity_ == subEntity_);
    }

    inline bool operator!= ( const ThisType& other ) const
    {
      return !(this->operator==( other ));
    }
 
    
    inline static ThisType begin( const LagrangePointSetType &lagrangePointSet,
                                  unsigned int subEntity )
    {
      return ThisType( lagrangePointSet, subEntity, true );
    }

    inline static ThisType end( const LagrangePointSetType &lagrangePointSet,
                                unsigned int subEntity )
    {
      return ThisType( lagrangePointSet, subEntity, false );
    }

  private:
    inline void assertDof ()
    {
      assert( lagrangePointSet_ != NULL );
        
      while( dofNumber_ >= numDofs_ ) {
        const ReferenceElementType &refElement
          = refElementContainer_( lagrangePointSet_->geometryType() );
          
        dofNumber_ = 0;
        ++subIndex_;
        while( subIndex_ >= numSubIndices_ ) {
          subIndex_ = 0;
          if( ++codim_ > dimension )
            return;
          numSubIndices_ = refElement.size( subEntity_, codimension, codim_ );
        }
        subSubEntity_ = refElement.subEntity( subEntity_, codimension,
                                              subIndex_, codim_ );
        numDofs_ = lagrangePointSet_->numDofs( codim_, subSubEntity_ );
      }
    }
  };
 


  template< class FieldImp,
            unsigned int dim,
            unsigned int polOrder >
  class SubEntityLagrangePointIterator< FieldImp, dim, 0, polOrder >
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = dim };
    enum { codimension = 0 };

    enum { polynomialOrder = polOrder };

    typedef FieldVector< FieldType, dimension > pointType;

  private:
    typedef SubEntityLagrangePointIterator< FieldType,
                                            dimension,
                                            codimension,
                                            polynomialOrder >
      ThisType;

    typedef LagrangePointSetInterface< FieldType, dimension, polynomialOrder >
      LagrangePointSetType;

  private:
    const LagrangePointSetType *lagrangePointSet_;
    unsigned int index_, numDofs_;

  private:
    inline SubEntityLagrangePointIterator
      ( const LagrangePointSetType &lagrangePointSet,
        const unsigned int subEntity,
        const bool beginIterator )
    {
      lagrangePointSet_ = &lagrangePointSet;
      assert( subEntity == 0 );

      numDofs_ = lagrangePointSet_->size();
      index_ = (beginIterator ? 0 : numDofs_);
    }
  
  public:
    inline SubEntityLagrangePointIterator ()
    {
      lagrangePointSet_ = NULL;
      index_ = 0;
      numDofs_ = 0;
    }
    
    inline SubEntityLagrangePointIterator ( const ThisType &other )
    {
      lagrangePointSet_ = other.lagrangePointSet_;

      index_ = other.index_;
      numDofs_ = other.numDofs_;
    }

    inline ThisType& operator= ( const ThisType &other )
    {
      lagrangePointSet_ = other.lagrangePointSet_;

      index_ = other.index_;
      numDofs_ = other.numDofs_;
    }

    inline unsigned int operator* () const
    {
      assert( lagrangePointSet_ != NULL );
      assert( index_ < numDofs_ );
      return index_;
    }

    inline ThisType& operator++ ()
    {
      if( index_ < numDofs_ )
        ++index_;
      return *this;
    }

    inline bool operator== ( const ThisType& other ) const
    {
      return (other.index_ == index_)
             && (other.lagrangePointSet_ == lagrangePointSet_);
    }

    inline bool operator!= ( const ThisType& other ) const
    {
      return !(this->operator==( other ));
    }
 
    
    inline static ThisType begin( const LagrangePointSetType &lagrangePointSet,
                                  unsigned int subEntity )
    {
      return ThisType( lagrangePointSet, subEntity, true );
    }

    inline static ThisType end( const LagrangePointSetType &lagrangePointSet,
                                unsigned int subEntity )
    {
      return ThisType( lagrangePointSet, subEntity, false );
    }
  };


  
  template< class FieldImp, unsigned int dim, unsigned int polOrder >
  class LagrangePointSetInterface
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = dim };

    enum { polynomialOrder = polOrder };

    typedef FieldVector< FieldType, dimension > PointType;

    template< unsigned int codim >
    struct Codim
    {
      //! type of iterator over DoF numbers in a subentity
      typedef SubEntityLagrangePointIterator< FieldType, 
                                              dimension,
                                              codim,
                                              polynomialOrder >
        SubEntityIteratorType;
    };
   
  private:
    typedef LagrangePointSetInterface< FieldType, dimension, polynomialOrder >
      ThisType;

  public:
    //! constructor
    LagrangePointSetInterface ()
    {
    }

  private:
    // Disallow the copy constructor
    LagrangePointSetInterface ( const ThisType &pointSet )
    {
    }

  public:
    //! destructor
    virtual ~LagrangePointSetInterface ()
    {
    }

    //! obtain a Lagrange point
    virtual const PointType& operator[] ( unsigned int index ) const = 0;

    virtual void dofSubEntity ( unsigned int index,
                                unsigned int &codim,
                                unsigned int &subEntity ) const = 0;
 
    virtual void dofSubEntity ( unsigned int index,
                                unsigned int &codim,
                                unsigned int &subEntity,
                                unsigned int &dofNumber ) const = 0;
 
    virtual unsigned int entityDofNumber ( unsigned int codim,
                                           unsigned int subEntity,
                                           unsigned int dofNumber ) const = 0;

    virtual const GeometryType geometryType () const = 0;

    virtual unsigned int maxDofs ( unsigned int codim ) const = 0;
    
    virtual unsigned int numDofs ( unsigned int codim,
                                   unsigned int subEntity ) const = 0;

    //! quadrature-like interface to obtain a Lagrange point
    inline const PointType& point( unsigned int index ) const
    {
      return this->operator[]( index );
    }
   
    //! get number of Lagrange points
    virtual unsigned int size () const = 0;

  public:
    template< unsigned int codim >
    inline typename Codim< codim > :: SubEntityIteratorType
      beginSubEntity ( unsigned int subEntity ) const
    {
      return Codim< codim > :: SubEntityIteratorType :: begin( *this, subEntity );
    }

    template< unsigned int codim >
    inline typename Codim< codim > :: SubEntityIteratorType
      endSubEntity ( unsigned int subEntity ) const
    {
      return Codim< codim > :: SubEntityIteratorType :: end( *this, subEntity );
    }
  };



  template< class FieldImp,
            GeometryType :: BasicType type,
            unsigned int dim,
            unsigned int polOrder >
  class LagrangePointSet
  : public LagrangePointSetInterface< FieldImp, dim, polOrder >
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = dim };

    enum { polynomialOrder = polOrder };
    
    typedef LagrangePoint< type, dimension, polynomialOrder >
      LagrangePointType;

    enum { numLagrangePoints = LagrangePointType :: numLagrangePoints };
    
  private:
    typedef LagrangePointSet< FieldType, type, dimension, polynomialOrder >
      ThisType;
    typedef LagrangePointSetInterface< FieldType, dimension, polynomialOrder >
      BaseType;

  public:
    typedef typename BaseType :: PointType PointType;

  private:
     PointType points_[ numLagrangePoints ];
     unsigned int codim_[ numLagrangePoints ];
     unsigned int subEntity_[ numLagrangePoints ];
     unsigned int dofNumber_[ numLagrangePoints ];

  public:
    LagrangePointSet ()
    : BaseType()
    {
      for( unsigned int i = 0; i < numLagrangePoints; ++i ) {
        LagrangePointType pt( i );
        pt.local( points_[ i ] );
        pt.dofSubEntity( codim_[ i ],
                         subEntity_[ i ],
                         dofNumber_[ i ] );
      }
    }

  private:
    // Disallow the copy constructor
    LagrangePointSet ( const ThisType &pointSet )
    {
    }

  public:
    virtual ~LagrangePointSet ()
    {
    }

    virtual const PointType& operator[] ( unsigned int index ) const
    {
      assert( index < numLagrangePoints );
      return points_[ index ];
    }

    virtual void dofSubEntity ( unsigned int index,
                                unsigned int &codim,
                                unsigned int &subEntity ) const
    {
      assert( index < numLagrangePoints );
      codim = codim_[ index ];
      subEntity = subEntity_[ index ];
    }
 
    virtual void dofSubEntity ( unsigned int index,
                                unsigned int &codim,
                                unsigned int &subEntity,
                                unsigned int &dofNumber ) const
    {
      assert( index < numLagrangePoints );
      codim = codim_[ index ];
      subEntity = subEntity_[ index ];
      dofNumber = dofNumber_[ index ];
    }
    
    virtual unsigned int entityDofNumber ( unsigned int codim,
                                           unsigned int subEntity,
                                           unsigned int dofNumber ) const
    {
      return LagrangePointType :: entityDofNumber
               ( codim, subEntity, dofNumber );
    }

    virtual const GeometryType geometryType () const
    {
      return GeometryType( type, dimension );
    }

    virtual unsigned int maxDofs ( unsigned int codim ) const
    {
      return LagrangePointType :: maxDofs( codim );
    }

    virtual unsigned int numDofs ( unsigned int codim,
                                   unsigned int subEntity ) const
    {
      return LagrangePointType :: numDofs( codim, subEntity );
    }

    virtual unsigned int size () const
    {
      return numLagrangePoints;
    }
  };



  template< class FieldImp, unsigned int dim, unsigned int polOrder >
  class LagrangePointSetFactory
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = dim };

    enum { polynomialOrder = polOrder };

    typedef LagrangePointSetInterface< FieldType, dimension, polynomialOrder >
      LagrangePointSetType;

  private:
    typedef LagrangePointSetFactory< FieldType, dimension, polynomialOrder >
      ThisType;

  public:
    inline static LagrangePointSetType*
       createObject ( const GeometryType type )
    {
      const GeometryType :: BasicType basicType = type.basicType();
      
      switch( basicType ) {
      case GeometryType :: simplex:
        return new LagrangePointSet< FieldType,
                                     GeometryType :: simplex, dimension,
                                     polynomialOrder >();

      case GeometryType :: cube:
        return new LagrangePointSet< FieldType,
                                     GeometryType :: cube, dimension,
                                     polynomialOrder >();

      default:
        DUNE_THROW( NotImplemented, "No such geometry type implemented." );
      }
    }

    inline static void deleteObject( LagrangePointSetType *pointSet )
    {
      delete pointSet;
    }
  };



  template< class FieldImp, unsigned int polOrder >
  class LagrangePointSetFactory< FieldImp, 3, polOrder >
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = 3 };

    enum { polynomialOrder = polOrder };

    typedef LagrangePointSetInterface< FieldType, dimension, polynomialOrder >
      LagrangePointSetType;

  private:
    typedef LagrangePointSetFactory< FieldType, dimension, polynomialOrder >
      ThisType;

  public:
    inline static LagrangePointSetType*
      createObject ( const GeometryType &type )
    {
      const GeometryType :: BasicType basicType = type.basicType();
      
      switch( basicType ) {
      case GeometryType :: simplex:
        return new LagrangePointSet< FieldType,
                                     GeometryType :: simplex, dimension,
                                     polynomialOrder >();

      case GeometryType :: cube:
        return new LagrangePointSet< FieldType,
                                     GeometryType :: cube, dimension,
                                     polynomialOrder >();

       
      case GeometryType :: pyramid:
        return new LagrangePointSet< FieldType,
                                     GeometryType :: pyramid, dimension,
                                     polynomialOrder >();

      
      case GeometryType :: prism:
        return new LagrangePointSet< FieldType,
                                     GeometryType :: prism, dimension,
                                     polynomialOrder >();

      default:
        DUNE_THROW( NotImplemented, "No such geometry type implemented." );
      }
    }

    inline static void deleteObject( LagrangePointSetType *pointSet )
    {
      delete pointSet;
    }
  };

}

#endif
