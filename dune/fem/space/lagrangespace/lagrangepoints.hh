#ifndef DUNE_LAGRANGESPACE_LAGRANGEPOINTS_HH
#define DUNE_LAGRANGESPACE_LAGRANGEPOINTS_HH

#if HAVE_DUNE_GEOMETRY
#include <dune/geometry/genericreferenceelements.hh>
#else
#include <dune/grid/common/genericreferenceelements.hh>
#endif
#include <dune/fem/quadrature/cachingpointlist.hh>

#include "genericgeometry.hh"
#include "genericlagrangepoints.hh"

namespace Dune
{

  /** @ingroup LagrangeDiscreteFunctionSpace 
   * \brief A single lagrange point
   *
   * A single Lagrange point which has methods
   * to determine on which subentity a given
   * lagrange point is located. 
   *
   **/
  template< unsigned int topologyId, unsigned int dim, unsigned int polOrder >
  class LagrangePoint
  : public GenericLagrangePoint< typename GeometryWrapper< topologyId, dim >::GenericGeometryType, polOrder >
  {
    typedef LagrangePoint< topologyId, dim, polOrder > ThisType;
    typedef GenericLagrangePoint< typename GeometryWrapper< topologyId, dim >::GenericGeometryType, polOrder > BaseType;

  public:
    static const unsigned int dimension = BaseType::dimension;

    typedef typename BaseType::DofCoordinateType DofCoordinateType;

    static const unsigned int polynomialOrder = BaseType::polynomialOrder;
    static const unsigned int numLagrangePoints = BaseType::numLagrangePoints;

  public:
    template< unsigned int codim >
    struct Codim
    {
      static unsigned int maxDofs ()
      {
        return BaseType::template Codim< codim >::maxDofs();
      }
    };
    
  public:
    LagrangePoint ( unsigned int index )
    : BaseType( index )
    {}

    LagrangePoint ( const BaseType &point )
    : BaseType( point )
    {}

    void dofSubEntity ( unsigned int &codim, unsigned int &subEntity )
    {
      BaseType::dofSubEntity( codim, subEntity );
    }

    void dofSubEntity ( unsigned int &codim, unsigned int &subEntity, unsigned int &dofNumber )
    {
      BaseType::dofSubEntity( codim, subEntity, dofNumber );
    }

    static unsigned int
    entityDofNumber ( unsigned int codim, unsigned int subEntity, unsigned int dof )
    {
      return BaseType::entityDofNumber( codim, subEntity, dof );
    }
  };

  template< unsigned int dim, unsigned int maxPolOrder >
  class LagrangePointInterface 
  {
    typedef LagrangePointInterface< dim, maxPolOrder > ThisType;
  protected:
    LagrangePointInterface() {}

  public:
    static const unsigned int dimension = dim;

    static const unsigned int maxPolynomialOrder = maxPolOrder;

    //! destructor 
    virtual ~LagrangePointInterface() {}

    virtual unsigned int entityDofNumber ( unsigned int codim,
                                           unsigned int subEntity,
                                           unsigned int dofNumber ) const = 0;
    
    virtual GeometryType geometryType () const = 0;

    /** \brief obtain the maximal number of DoFs in one entity of a codimension
     *
     *  \param[in]  codim  codimension, the information is desired for
     *
     *  \returns maximal number of DoFs for one entity in the codimension
     */
    virtual unsigned int maxDofs ( unsigned int codim ) const = 0;

    inline static int maxOrder ()
    {
      return maxPolynomialOrder;
    }

    /** \brief obtain the number of DoFs on one entity
     *
     *  \param[in]  codim      codimension of the entity
     *  \param[in]  subEntity  number of the subentity (of the given codimension)
     *
     *  \returns the number of DoFs associated with the specified entity
     */
    virtual unsigned int numDofs ( unsigned int codim, unsigned int subEntity ) const = 0;

    /** \brief obtain the total number of DoFs in a codimension
     *
     *  \param[in]  codim      codimension the information is desired for
     *
     *  \returns the number of DoFs associated with the codimension
     */
    virtual unsigned int numDofs ( unsigned int codim ) const = 0;

    virtual int order () const
    {
      return maxPolynomialOrder;
    }
  };

  template< unsigned int topologyId, unsigned int dim, unsigned int maxPolOrder, int polOrder  >
  class LagrangePointImplementation : 
    public LagrangePointInterface< dim, maxPolOrder >
  {
    typedef LagrangePoint< topologyId, dim, polOrder > LagrangePointType;
  public:
    LagrangePointImplementation() {}

    virtual ~LagrangePointImplementation() {}

    virtual unsigned int
    entityDofNumber ( unsigned int codim, unsigned int subEntity, unsigned int dofNumber ) const
    {
      return LagrangePointType::entityDofNumber( codim, subEntity, dofNumber );
    }

    virtual GeometryType geometryType () const
    {
      return GeometryType( topologyId, dim );
    }

    /** \copydoc Dune::LagrangePointListInterface::maxDofs
     */
    virtual unsigned int maxDofs ( unsigned int codim ) const
    {
      return LagrangePointType::maxDofs( codim );
    }

    /** \copydoc Dune::LagrangePointListInterface::numDofs(unsigned int,unsigned int)
     */
    virtual unsigned int
    numDofs ( unsigned int codim, unsigned int subEntity ) const
    {
      return LagrangePointType::numDofs( codim, subEntity );
    }

    /** \copydoc Dune::LagrangePointListInterface::numDofs(unsigned int)
     */
    virtual unsigned int numDofs ( unsigned int codim ) const
    {
      return LagrangePointType::numDofs( codim );
    }

    virtual int order () const { return polOrder; }
  };

  /** @ingroup LagrangeDiscreteFunctionSpace 
   * \brief Set of lagrange points
   *
   * Interface class for a set of lagrange points.
   * An instance of the lagrange points can be
   * obtained from the 
   * \ref Dune::LagrangeDiscreteFunctionSpace "lagrange function space".
   * The set can be wrapped in a quadrature. 
   *
   **/
  template< class FieldImp, unsigned int dim, unsigned int maxPolOrder >
  class LagrangePointListInterface
  : public IntegrationPointListImp< FieldImp, dim >
  {
  public:
    //! field type of points
    typedef FieldImp FieldType;

    //! dimension of points
    enum { dimension = dim };
    
    //! polynomial order of corresponding base functions
    enum { maxPolynomialOrder = maxPolOrder };

    //! type of points
    typedef FieldVector< FieldType, dimension > CoordinateType;

    struct DofInfo
    {
      unsigned int codim;
      unsigned int subEntity;
      unsigned int dofNumber;
    };

  private:
    typedef LagrangePointListInterface< FieldType, dimension, maxPolynomialOrder >
      ThisType;
    typedef IntegrationPointListImp< FieldType, dimension > BaseType;

    typedef LagrangePointInterface< dim, maxPolynomialOrder >
      LagrangePointInterfaceType;
  private:
    std :: vector< DofInfo > dofInfos_;
    const LagrangePointInterfaceType* lagrangePointImpl_;

    const LagrangePointInterfaceType& lagrangePointImpl() const
    {
      assert( lagrangePointImpl_ ) ;
      return *lagrangePointImpl_;
    }

  public:
    LagrangePointListInterface( const size_t id )
    : BaseType( id ),
      dofInfos_(),
      lagrangePointImpl_( 0 )
    {}

    ~LagrangePointListInterface() 
    {
      delete lagrangePointImpl_;
    }
  private:
    // prohibit copy construction  
    LagrangePointListInterface ( const ThisType &other );

  public:
    void setLagrangePointImpl( const LagrangePointInterfaceType* lpImpl ) 
    {
      assert( lagrangePointImpl_ == 0 );
      lagrangePointImpl_ = lpImpl ;
    }

    inline const DofInfo& dofInfo( unsigned int index ) const
    {
      return dofInfos_[ index ];
    }
    
    inline void dofSubEntity( unsigned int index,
                              unsigned int &codim, 
                              unsigned int &subEntity,
                              unsigned int &dofNumber ) const
    {
      const DofInfo &dofInfo = dofInfos_[ index ];
      codim = dofInfo.codim;
      subEntity = dofInfo.subEntity;
      dofNumber = dofInfo.dofNumber;
    }
    
    inline unsigned int entityDofNumber ( unsigned int codim,
                                          unsigned int subEntity,
                                          unsigned int dofNumber ) const 
    {
      return lagrangePointImpl().entityDofNumber( codim, subEntity, dofNumber ); 
    }
    
    inline GeometryType geometryType () const 
    {
      return lagrangePointImpl().geometryType();
    }

    /** \brief obtain the maximal number of DoFs in one entity of a codimension
     *
     *  \param[in]  codim  codimension, the information is desired for
     *
     *  \returns maximal number of DoFs for one entity in the codimension
     */
    inline unsigned int maxDofs ( unsigned int codim ) const 
    {
      return lagrangePointImpl().maxDofs( codim );
    }

    inline static int maxOrder ()
    {
      return LagrangePointInterfaceType :: maxOrder();
    }

    /** \brief obtain the number of DoFs on one entity
     *
     *  \param[in]  codim      codimension of the entity
     *  \param[in]  subEntity  number of the subentity (of the given codimension)
     *
     *  \returns the number of DoFs associated with the specified entity
     */
    inline unsigned int numDofs ( unsigned int codim, unsigned int subEntity ) const 
    {
      return lagrangePointImpl().numDofs( codim, subEntity );
    }

    /** \brief obtain the total number of DoFs in a codimension
     *
     *  \param[in]  codim      codimension the information is desired for
     *
     *  \returns the number of DoFs associated with the codimension
     */
    inline unsigned int numDofs ( unsigned int codim ) const 
    {
      return lagrangePointImpl().numDofs( codim ); 
    }

    inline int order () const
    {
      return lagrangePointImpl().order();
    }

  protected:
    inline void addDofInfo( const DofInfo &dofInfo )
    {
      dofInfos_.push_back( dofInfo );
    }
  };


  template< class FieldImp, unsigned int topologyId, unsigned int dim, unsigned int maxPolOrder >
  class LagrangePointListImplementation
  : public LagrangePointListInterface< FieldImp, dim, maxPolOrder >
  {
    typedef LagrangePointListImplementation< FieldImp, topologyId, dim, maxPolOrder > ThisType;
    typedef LagrangePointListInterface< FieldImp, dim, maxPolOrder > BaseType;

  public:
    //! field type of points
    typedef FieldImp FieldType;

    //! dimension of points
    enum { dimension = dim };

    //! polynomial order of corresponding base functions
    static const int maxPolynomialOrder = maxPolOrder ;

    //! type of points
    typedef FieldVector< FieldType, dimension > CoordinateType;

  private:
    template <int pOrd>
    struct CreateLagrangePoint  
    {
      typedef LagrangePoint< topologyId, dimension, pOrd > LagrangePointType;
      enum { numLagrangePoints = LagrangePointType::numLagrangePoints };

      static void apply( ThisType& lp, const int order ) 
      {
        // if order is not equal to pOrd, do nothing 
        if( order != pOrd ) return ;

        for( unsigned int i = 0; i < numLagrangePoints; ++i )
        {
          LagrangePointType pt( i );
          
          CoordinateType local;
          pt.local( local );
          lp.addIntegrationPoint( local );
          
          typename BaseType :: DofInfo dofInfo;
          pt.dofSubEntity( dofInfo.codim, dofInfo.subEntity, dofInfo.dofNumber );
          lp.addDofInfo( dofInfo );
        }

        typedef LagrangePointImplementation< topologyId, dim, maxPolynomialOrder, pOrd >
            LagrangePointImplementationType;
        lp.setLagrangePointImpl( new LagrangePointImplementationType() );
      }
    };
  public:
    LagrangePointListImplementation ( const size_t id )
    : BaseType( id )
    {
      ForLoop< CreateLagrangePoint, 1, maxPolynomialOrder > :: apply( *this, maxPolynomialOrder );
    }
    
    LagrangePointListImplementation ( const GeometryType &geo, const int order, const size_t id )
    : BaseType( id )
    {
      ForLoop< CreateLagrangePoint, 1, maxPolynomialOrder > :: apply( *this, order );

      // assert this after lagrangePointImpl has been created since 
      // this->geometry() uses this class
      assert( order <= maxPolynomialOrder );
      assert( geo == this->geometryType() );
    }

  private:
    LagrangePointListImplementation ( const ThisType &other );
  };


  template< class FieldImp, unsigned int topologyId, unsigned int dim, unsigned int maxPolOrder >
  const int LagrangePointListImplementation< FieldImp, topologyId, dim, maxPolOrder >::maxPolynomialOrder;




  template< class GridPartImp, unsigned int polOrder >
  struct LagrangePointSetTraits
  {
    //! type of grid partition
    typedef GridPartImp GridPartType;

    //! polynomial order of corresponding base functions
    enum { polynomialOrder = polOrder };

    //! type of grid
    typedef typename GridPartType :: GridType GridType;
    
    //! field type of coordinates
    typedef typename GridType :: ctype FieldType;

    //! dimension of points
    enum { dimension = GridType :: dimension };

    //! codimension of point set
    enum { codimension = 0 };

    //! default defines for used point lists
    template< typename ct, int dim >
    struct PointListTraits
    {
      typedef LagrangePointListImplementation< ct, 0, 0, polynomialOrder >
        PointQuadratureType;

      typedef LagrangePointListImplementation< ct, 0, 1, polynomialOrder >
        LineQuadratureType;
      
      typedef LagrangePointListImplementation< ct, 0, dimension, polynomialOrder >
        SimplexQuadratureType;
 
      typedef LagrangePointListImplementation< ct, (1 << dimension)-1, dimension, polynomialOrder >
        CubeQuadratureType;
      
      typedef LagrangePointListImplementation< ct, (1 << (dimension-1)), dimension, polynomialOrder >
        PrismQuadratureType;
      
      typedef LagrangePointListImplementation< ct, (1 << (dimension-1))-1, dimension, polynomialOrder >
        PyramidQuadratureType;
      
      //! type of integration point list implemementation 
      typedef LagrangePointListInterface< ct, dim, polynomialOrder > IntegrationPointListType; 
    }; 

    //! type of used integration point list 
    typedef IntegrationPointList< FieldType, dimension, PointListTraits > IntegrationPointListType;

    //! type of global coordinate 
    typedef typename IntegrationPointListType :: CoordinateType CoordinateType;
  };



  template< class GridPartImp, unsigned int polOrder >
  class LagrangePointSet;

  

  template< class GridPartImp,
            unsigned int codim,
            unsigned int polOrder >
  class SubEntityLagrangePointIterator
  {
  public:
    typedef GridPartImp GridPartType;

    typedef typename GridPartType :: GridType GridType;

    typedef typename GridType :: ctype FieldType;

    enum { dimension = GridType :: dimension };
    enum { codimension = codim };

    enum { polynomialOrder = polOrder };

    typedef FieldVector< FieldType, dimension > pointType;

  private:
    typedef SubEntityLagrangePointIterator< GridPartType,
                                            codimension,
                                            polynomialOrder >
      ThisType;

    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;

    typedef GenericReferenceElementContainer< FieldType, dimension >
      ReferenceElementContainerType;
    typedef typename ReferenceElementContainerType :: value_type
      ReferenceElementType;

  private:
    const ReferenceElementContainerType *refElementContainer_;
    
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
    : refElementContainer_( &ReferenceElementContainerType :: instance () )
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
    : refElementContainer_( &ReferenceElementContainerType :: instance () )
    {
      lagrangePointSet_ = NULL;
      subEntity_ = 0;

      codim_ = dimension + 1;
      subIndex_ = 0;
      dofNumber_ = 0;
    }
    
    inline SubEntityLagrangePointIterator ( const ThisType &other )
    : refElementContainer_( &ReferenceElementContainerType :: instance () )
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
    const ReferenceElementType & referenceElement ( const GeometryType& geo) const
    {
      assert( refElementContainer_ );
      return (*refElementContainer_)( geo );
    }

    inline void assertDof ()
    {
      assert( lagrangePointSet_ != NULL );
        
      while( dofNumber_ >= numDofs_ ) {
        const ReferenceElementType &refElement = referenceElement( lagrangePointSet_->geometryType() );
          
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
 


  template< class GridPartImp,
            unsigned int polOrder >
  class SubEntityLagrangePointIterator< GridPartImp, 0, polOrder >
  {
  public:
    typedef GridPartImp GridPartType;

    typedef typename GridPartType :: GridType GridType;

    typedef typename GridType :: ctype FieldType;

    enum { dimension = GridType :: dimension };
    enum { codimension = 0 };

    enum { polynomialOrder = polOrder };

    typedef FieldVector< FieldType, dimension > pointType;

  private:
    typedef SubEntityLagrangePointIterator< GridPartType,
                                            codimension,
                                            polynomialOrder >
      ThisType;

    typedef LagrangePointSet< GridPartType, polynomialOrder >
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



  template< class GridPartImp, unsigned int maxPolOrder >
  class LagrangePointSet
  : public CachingPointList< GridPartImp, 0, 
                             LagrangePointSetTraits< GridPartImp, maxPolOrder > >
  {
  public:
    typedef LagrangePointSetTraits< GridPartImp, maxPolOrder > Traits;

    typedef typename Traits :: GridPartType GridPartType;

    typedef typename Traits :: FieldType FieldType;

    enum { dimension = Traits :: dimension };
    
    enum { polynomialOrder = Traits :: polynomialOrder };

    typedef typename Traits :: CoordinateType CoordinateType;
    typedef typename Traits :: CoordinateType PointType;

    template< unsigned int codim >
    struct Codim
    {
      //! type of iterator over DoF numbers in a subentity
      typedef SubEntityLagrangePointIterator< GridPartType, 
                                              codim,
                                              polynomialOrder >
        SubEntityIteratorType;
    };
   
  private:
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      ThisType;
    typedef CachingPointList< GridPartType, 0, Traits > BaseType;

    typedef typename BaseType :: IntegrationPointListType
                              :: IntegrationPointListType
      LagrangePointListType;

  public:
    typedef typename LagrangePointListType :: DofInfo DofInfo;

  private:
    const LagrangePointListType &lagrangePointList_;

  public:
    //! constructor
    inline LagrangePointSet ( const GeometryType &geometry, const int order )
    : BaseType( geometry, order ),
      lagrangePointList_( this->quadImp().ipList() )
    {
    }

    //! copy constructor
    inline LagrangePointSet ( const ThisType &other )
    : BaseType( other ),
      lagrangePointList_( this->quadImp().ipList() )
    {
    }

  private:
    // kill the assignment operator
    inline ThisType& operator=( const ThisType &other )
    {
      assert( false );
      abort();
    }

  public:
    inline const DofInfo& dofInfo( unsigned int index ) const
    {
      return lagrangePointList_.dofInfo( index );
    }

    inline void dofSubEntity ( unsigned int index,
                               unsigned int &codim,
                               unsigned int &subEntity ) const
    {
      unsigned int dofNumber;
      lagrangePointList_.dofSubEntity
               ( index, codim, subEntity, dofNumber );
    }
 
    inline void dofSubEntity ( unsigned int index,
                               unsigned int &codim,
                               unsigned int &subEntity,
                               unsigned int &dofNumber ) const
    {
      lagrangePointList_.dofSubEntity
             ( index, codim, subEntity, dofNumber );
    }
 
    inline unsigned int entityDofNumber ( unsigned int codim,
                                          unsigned int subEntity,
                                          unsigned int dofNumber ) const
    {
      return lagrangePointList_.entityDofNumber( codim, subEntity, dofNumber );
    }

    inline unsigned int maxDofs ( unsigned int codim ) const
    {
      return lagrangePointList_.maxDofs( codim );
    }
    
    inline unsigned int numDofs ( unsigned int codim,
                                  unsigned int subEntity ) const
    {
      return lagrangePointList_.numDofs( codim, subEntity );
    }

    inline unsigned int numDofs ( unsigned int codim ) const
    {
      return lagrangePointList_.numDofs( codim );
    }

    //! get number of Lagrange points
    inline unsigned int size () const
    {
      return this->nop();
    }

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

} // end namespace dune 
#endif
