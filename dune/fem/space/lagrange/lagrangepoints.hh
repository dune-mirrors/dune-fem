#ifndef DUNE_FEM_SPACE_LAGRANGE_LAGRANGEPOINTS_HH
#define DUNE_FEM_SPACE_LAGRANGE_LAGRANGEPOINTS_HH

// dune-common includes
#include <dune/fem/common/forloop.hh>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>

// dune-fem includes
#include <dune/fem/quadrature/cachingpointlist.hh>
#include <dune/fem/space/mapper/localkey.hh>

// local includes
#include "genericgeometry.hh"
#include "genericlagrangepoints.hh"


namespace Dune
{

  namespace Fem
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
    : public GenericLagrangePoint< typename GeometryWrapper< topologyId, dim >::ImplType, polOrder >
    {
      typedef LagrangePoint< topologyId, dim, polOrder > ThisType;
      typedef GenericLagrangePoint< typename GeometryWrapper< topologyId, dim >::ImplType, polOrder > BaseType;

    public:
      static const unsigned int dimension = BaseType::dimension;

      typedef typename BaseType::DofCoordinateType DofCoordinateType;

      static const unsigned int polynomialOrder = BaseType::polynomialOrder;
      static const unsigned int numLagrangePoints = BaseType::numLagrangePoints;

      template< unsigned int codim >
      struct Codim
      {
        static unsigned int maxDofs ()
        {
          return BaseType::template Codim< codim >::maxDofs();
        }
      };

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
      LagrangePointInterface() = default;

    public:
      static const unsigned int dimension = dim;

      static const unsigned int maxPolynomialOrder = maxPolOrder;

      //! destructor
      virtual ~LagrangePointInterface() = default;

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

      static int maxOrder ()
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
      LagrangePointImplementation() = default;

      virtual ~LagrangePointImplementation() = default;

      virtual unsigned int
      entityDofNumber ( unsigned int codim, unsigned int subEntity, unsigned int dofNumber ) const
      {
        return LagrangePointType::entityDofNumber( codim, subEntity, dofNumber );
      }

      virtual GeometryType geometryType () const
      {
        return GeometryType( topologyId, dim );
      }

      /** \copydoc Dune::Fem::LagrangePointListInterface::maxDofs
       */
      virtual unsigned int maxDofs ( unsigned int codim ) const
      {
        return LagrangePointType::maxDofs( codim );
      }

      /** \copydoc Dune::Fem::LagrangePointListInterface::numDofs(unsigned int,unsigned int)
       */
      virtual unsigned int
      numDofs ( unsigned int codim, unsigned int subEntity ) const
      {
        return LagrangePointType::numDofs( codim, subEntity );
      }

      /** \copydoc Dune::Fem::LagrangePointListInterface::numDofs(unsigned int)
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
     * \ref Dune::Fem::LagrangeDiscreteFunctionSpace "lagrange function space".
     * The set can be wrapped in a quadrature.
     *
     **/
    template< class FieldImp, int dim, unsigned int maxPolOrder >
    class LagrangePointListInterface
    : public IntegrationPointListImp< FieldImp, dim >
    {
      typedef LagrangePointListInterface< FieldImp, dim, maxPolOrder > ThisType;
      typedef IntegrationPointListImp< FieldImp, dim > BaseType;

    public:
      //! field type of points
      typedef FieldImp FieldType;

      //! dimension of points
      static const int dimension = dim;

      //! polynomial order of corresponding base functions
      static const unsigned int maxPolynomialOrder = maxPolOrder;

      //! type of points
      typedef FieldVector< FieldType, dimension > CoordinateType;

    private:
      typedef LagrangePointInterface< dim, maxPolynomialOrder >
        LagrangePointInterfaceType;

      const LagrangePointInterfaceType &lagrangePointImpl () const
      {
        assert( lagrangePointImpl_ ) ;
        return *lagrangePointImpl_;
      }

    public:
      explicit LagrangePointListInterface ( const size_t id )
      : BaseType( id ),
        dofInfos_(),
        lagrangePointImpl_( 0 )
      {}

      ~LagrangePointListInterface ()
      {
        delete lagrangePointImpl_;
      }

      LagrangePointListInterface ( const ThisType& ) = delete;

      void setLagrangePointImpl ( const LagrangePointInterfaceType* lpImpl )
      {
        assert( lagrangePointImpl_ == 0 );
        lagrangePointImpl_ = lpImpl ;
      }

      const LocalKey &dofInfo ( unsigned int index ) const
      {
        return dofInfos_[ index ];
      }

      void dofSubEntity ( unsigned int index,
                          unsigned int &codim,
                          unsigned int &subEntity,
                          unsigned int &dofNumber ) const
      {
        const LocalKey &dofInfo = this->dofInfo( index );
        codim = dofInfo.codim();
        subEntity = dofInfo.subEntity();
        dofNumber = dofInfo.index();
      }

      unsigned int entityDofNumber ( unsigned int codim,
                                     unsigned int subEntity,
                                     unsigned int dofNumber ) const
      {
        return lagrangePointImpl().entityDofNumber( codim, subEntity, dofNumber );
      }

      GeometryType geometryType () const
      {
        return lagrangePointImpl().geometryType();
      }

      /** \brief obtain the maximal number of DoFs in one entity of a codimension
       *
       *  \param[in]  codim  codimension, the information is desired for
       *
       *  \returns maximal number of DoFs for one entity in the codimension
       */
      unsigned int maxDofs ( unsigned int codim ) const
      {
        return lagrangePointImpl().maxDofs( codim );
      }

      static int maxOrder ()
      {
        return LagrangePointInterfaceType::maxOrder();
      }

      /** \brief obtain the number of DoFs on one entity
       *
       *  \param[in]  codim      codimension of the entity
       *  \param[in]  subEntity  number of the subentity (of the given codimension)
       *
       *  \returns the number of DoFs associated with the specified entity
       */
      unsigned int numDofs ( unsigned int codim, unsigned int subEntity ) const
      {
        return lagrangePointImpl().numDofs( codim, subEntity );
      }

      /** \brief obtain the total number of DoFs in a codimension
       *
       *  \param[in]  codim      codimension the information is desired for
       *
       *  \returns the number of DoFs associated with the codimension
       */
      unsigned int numDofs ( unsigned int codim ) const
      {
        return lagrangePointImpl().numDofs( codim );
      }

      int order () const
      {
        return lagrangePointImpl().order();
      }

    protected:
      void addDofInfo ( const LocalKey &dofInfo )
      {
        dofInfos_.push_back( dofInfo );
      }

    private:
      std::vector< LocalKey > dofInfos_;
      const LagrangePointInterfaceType* lagrangePointImpl_;
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

            unsigned int codim, subEntity, dofNumber;
            pt.dofSubEntity( codim, subEntity, dofNumber );
            lp.addDofInfo( LocalKey( subEntity, codim, dofNumber ) );
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
        Fem::ForLoop< CreateLagrangePoint, 1, maxPolynomialOrder > :: apply( *this, maxPolynomialOrder );
      }

      LagrangePointListImplementation ( const GeometryType &geo, const int order, const size_t id )
      : BaseType( id )
      {
        Fem::ForLoop< CreateLagrangePoint, 1, maxPolynomialOrder > :: apply( *this, order );

        // assert this after lagrangePointImpl has been created since
        // this->geometry() uses this class
        assert( order <= maxPolynomialOrder );
        assert( geo == this->geometryType() );
      }

      LagrangePointListImplementation ( const ThisType& ) = delete;
    };


    template< class FieldImp, unsigned int topologyId, unsigned int dim, unsigned int maxPolOrder >
    const int LagrangePointListImplementation< FieldImp, topologyId, dim, maxPolOrder >::maxPolynomialOrder;




    template< class Field, int dim, unsigned int maxPolOrder >
    struct LagrangePointSetTraits
    {
      //! field type of coordinates
      typedef Field FieldType;

      //! polynomial order of corresponding base functions
      static const unsigned int maxPolynomialOrder = maxPolOrder;

      //! dimension of points
      static const int dimension = dim;

      //! codimension of point set
      static const int codimension = 0;

      //! default defines for used point lists
      template< typename ct, int quaddim >
      struct PointListTraits
      {
        typedef LagrangePointListImplementation< ct, 0, 0, maxPolynomialOrder >
          PointQuadratureType;

        typedef LagrangePointListImplementation< ct, 0, 1, maxPolynomialOrder >
          LineQuadratureType;

        typedef LagrangePointListImplementation< ct, 0, dimension, maxPolynomialOrder >
          SimplexQuadratureType;

        typedef LagrangePointListImplementation< ct, (1 << dimension)-1, dimension, maxPolynomialOrder >
          CubeQuadratureType;

        typedef LagrangePointListImplementation< ct, (1 << (dimension-1)), dimension, maxPolynomialOrder >
          PrismQuadratureType;

        typedef LagrangePointListImplementation< ct, (1 << (dimension-1))-1, dimension, maxPolynomialOrder >
          PyramidQuadratureType;

        //! type of integration point list implemementation
        typedef LagrangePointListInterface< ct, quaddim, maxPolynomialOrder > IntegrationPointListType;

        typedef int QuadratureKeyType ;
      };

      //! type of used integration point list
      typedef Fem::IntegrationPointList< FieldType, dimension, PointListTraits > IntegrationPointListType;

      //! type of global coordinate
      typedef typename IntegrationPointListType::CoordinateType CoordinateType;
    };



    template< class GridPart, unsigned int maxPolOrder >
    class LagrangePointSet;



    // SubEntityLagrangePointIterator
    // ------------------------------

    template< class GridPart, int codim, unsigned int polOrder >
    class SubEntityLagrangePointIterator
    {
      typedef SubEntityLagrangePointIterator< GridPart, codim, polOrder > ThisType;

    public:
      typedef GridPart GridPartType;

      typedef typename GridPartType::ctype FieldType;

      static const int dimension = GridPartType::dimension;
      static const int codimension = codim;

      static const unsigned int polynomialOrder = polOrder;

      typedef FieldVector< FieldType, dimension > pointType;

      typedef LagrangePointSet< GridPartType, polynomialOrder > LagrangePointSetType;

    private:
      typedef Dune::ReferenceElement< FieldType, dimension > ReferenceElementType;
      typedef Dune::ReferenceElements< FieldType, dimension > ReferenceElementsType;

      SubEntityLagrangePointIterator ( const LagrangePointSetType &lagrangePointSet,
                                       const unsigned int subEntity,
                                       const bool beginIterator )
      : lagrangePointSet_( &lagrangePointSet ),
        refElement_( &ReferenceElementsType::general( lagrangePointSet_->geometryType() ) ),
        subEntity_( subEntity ),
        codim_( beginIterator ? codimension : dimension+1 ),
        subIndex_( 0 ),
        numSubIndices_( 1 ),
        subSubEntity_( subEntity_ ),
        dofNumber_( 0 ),
        numDofs_( lagrangePointSet_->numDofs( codimension, subSubEntity_ ) )
      {
        if( beginIterator )
          assertDof();
      }

    public:
      SubEntityLagrangePointIterator ()
      : lagrangePointSet_( 0 ),
        refElement_( 0 ),
        codim_( dimension+1 ),
        subEntity_( 0 ),
        subIndex_( 0 ),
        dofNumber_( 0 )
      {}

      unsigned int operator* () const
      {
        assert( lagrangePointSet_ );
        assert( codim_ <= dimension );

        return lagrangePointSet_->entityDofNumber( codim_, subSubEntity_, dofNumber_ );
      }

      ThisType &operator++ ()
      {
        assert( codim_ <= dimension );
        ++dofNumber_;
        assertDof();
        return *this;
      }

      bool operator== ( const ThisType& other ) const
      {
        if( (other.codim_ != codim_)
            || (other.subIndex_ != subIndex_)
            || (other.dofNumber_ != dofNumber_) )
          return false;

        return (other.lagrangePointSet_ == lagrangePointSet_)
               && (other.subEntity_ == subEntity_);
      }

      bool operator!= ( const ThisType& other ) const
      {
        return !(*this == other);
      }

      static ThisType begin ( const LagrangePointSetType &lagrangePointSet,
                              unsigned int subEntity )
      {
        return ThisType( lagrangePointSet, subEntity, true );
      }

      static ThisType end ( const LagrangePointSetType &lagrangePointSet,
                            unsigned int subEntity )
      {
        return ThisType( lagrangePointSet, subEntity, false );
      }

    private:
      void assertDof ()
      {
        assert( lagrangePointSet_ );
        assert( refElement_ );

        while( dofNumber_ >= numDofs_ )
        {
          const ReferenceElementType &refElement = *refElement_;

          dofNumber_ = 0;
          ++subIndex_;
          while( subIndex_ >= numSubIndices_ )
          {
            subIndex_ = 0;
            if( ++codim_ > dimension )
              return;
            numSubIndices_ = refElement.size( subEntity_, codimension, codim_ );
          }
          subSubEntity_ = refElement.subEntity( subEntity_, codimension, subIndex_, codim_ );
          numDofs_ = lagrangePointSet_->numDofs( codim_, subSubEntity_ );
        }
      }

      const LagrangePointSetType *lagrangePointSet_;
      const ReferenceElementType *refElement_;
      unsigned int subEntity_;

      int codim_;
      unsigned int subIndex_, numSubIndices_;
      unsigned int subSubEntity_;
      unsigned int dofNumber_, numDofs_;
    };



    // SubEntityLagrangePointIterator for codimension 0
    // ------------------------------------------------

    template< class GridPart, unsigned int polOrder >
    class SubEntityLagrangePointIterator< GridPart, 0, polOrder >
    {
      typedef SubEntityLagrangePointIterator< GridPart, 0, polOrder > ThisType;

    public:
      typedef GridPart GridPartType;

      typedef typename GridPartType::ctype FieldType;

      static const int dimension = GridPartType::dimension;
      static const int codimension = 0;

      static const unsigned int polynomialOrder = polOrder;

      typedef FieldVector< FieldType, dimension > pointType;

      typedef LagrangePointSet< GridPartType, polynomialOrder > LagrangePointSetType;

    private:
      SubEntityLagrangePointIterator ( const LagrangePointSetType &lagrangePointSet,
                                       const unsigned int subEntity,
                                       const bool beginIterator )
      : lagrangePointSet_( &lagrangePointSet ),
        numDofs_( lagrangePointSet_->size() ),
        index_( beginIterator ? 0 : numDofs_ )
      {
        assert( subEntity == 0 );
      }

    public:
      SubEntityLagrangePointIterator ()
      : lagrangePointSet_( 0 ),
        numDofs_( 0 ),
        index_( 0 )
      {}

      unsigned int operator* () const
      {
        assert( lagrangePointSet_ );
        assert( index_ < numDofs_ );
        return index_;
      }

      ThisType& operator++ ()
      {
        assert( index_ < numDofs_ );
        ++index_;
        return *this;
      }

      bool operator== ( const ThisType& other ) const
      {
        return (other.index_ == index_)
               && (other.lagrangePointSet_ == lagrangePointSet_);
      }

      bool operator!= ( const ThisType& other ) const
      {
        return !(*this == other);
      }


      static ThisType begin ( const LagrangePointSetType &lagrangePointSet,
                              unsigned int subEntity )
      {
        return ThisType( lagrangePointSet, subEntity, true );
      }

      static ThisType end ( const LagrangePointSetType &lagrangePointSet,
                            unsigned int subEntity )
      {
        return ThisType( lagrangePointSet, subEntity, false );
      }

    private:
      const LagrangePointSetType *lagrangePointSet_;
      unsigned int numDofs_, index_;
    };



    template< class GridPart, unsigned int maxPolOrder >
    class LagrangePointSet
    : public CachingPointList< GridPart, 0, LagrangePointSetTraits< typename GridPart::ctype, GridPart::dimension, maxPolOrder > >
    {
      typedef LagrangePointSet< GridPart, maxPolOrder > ThisType;
      typedef CachingPointList< GridPart, 0, LagrangePointSetTraits< typename GridPart::ctype, GridPart::dimension, maxPolOrder > > BaseType;

    public:
      typedef LagrangePointSetTraits< typename GridPart::ctype, GridPart::dimension, maxPolOrder > Traits;

      typedef GridPart GridPartType;

      typedef typename GridPartType::ctype FieldType;

      static const int dimension = BaseType::dimension;

      enum { maxPolynomialOrder = Traits::maxPolynomialOrder };

      typedef typename BaseType::CoordinateType CoordinateType;
      typedef typename Traits::CoordinateType PointType;

      template< unsigned int codim >
      struct Codim
      {
        //! type of iterator over DoF numbers in a subentity
        typedef SubEntityLagrangePointIterator< GridPartType, codim, maxPolynomialOrder >
          SubEntityIteratorType;
      };

    private:
      typedef typename BaseType::IntegrationPointListType::IntegrationPointListType
        LagrangePointListType;

    public:
      //! constructor
      LagrangePointSet ( const GeometryType &geometry, const int polynomialOrder )
      : BaseType( geometry, polynomialOrder ), lagrangePointList_( this->quadImp().ipList() )
      {}

      //! copy constructor
      LagrangePointSet ( const ThisType &other )
      : BaseType( other ), lagrangePointList_( this->quadImp().ipList() )
      {}

      ThisType& operator=( const ThisType& ) = delete;

      const LocalKey &localKey ( unsigned int index ) const
      {
        return lagrangePointList_.dofInfo( index );
      }

      void dofSubEntity ( unsigned int index, unsigned int &codim,
                          unsigned int &subEntity ) const
      {
        unsigned int dofNumber;
        lagrangePointList_.dofSubEntity( index, codim, subEntity, dofNumber );
      }

      void dofSubEntity ( unsigned int index, unsigned int &codim,
                          unsigned int &subEntity, unsigned int &dofNumber ) const
      {
        lagrangePointList_.dofSubEntity( index, codim, subEntity, dofNumber );
      }

      unsigned int entityDofNumber ( unsigned int codim, unsigned int subEntity,
                                     unsigned int dofNumber ) const
      {
        return lagrangePointList_.entityDofNumber( codim, subEntity, dofNumber );
      }

      unsigned int maxDofs ( unsigned int codim ) const
      {
        return lagrangePointList_.maxDofs( codim );
      }

      unsigned int numDofs ( unsigned int codim, unsigned int subEntity ) const
      {
        return lagrangePointList_.numDofs( codim, subEntity );
      }

      unsigned int numDofs ( unsigned int codim ) const
      {
        return lagrangePointList_.numDofs( codim );
      }

      //! get number of Lagrange points
      std::size_t size () const
      {
        return this->nop();
      }

      template< unsigned int codim >
      typename Codim< codim >::SubEntityIteratorType
        beginSubEntity ( unsigned int subEntity ) const
      {
        return Codim< codim >::SubEntityIteratorType::begin( *this, subEntity );
      }

      template< unsigned int codim >
      typename Codim< codim >::SubEntityIteratorType
        endSubEntity ( unsigned int subEntity ) const
      {
        return Codim< codim >::SubEntityIteratorType::end( *this, subEntity );
      }

    private:
      const LagrangePointListType &lagrangePointList_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_LAGRANGEPOINTS_HH
