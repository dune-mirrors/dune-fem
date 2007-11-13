#ifndef DUNE_FEM_TEMPORARYLOCALFUNCTION_HH
#define DUNE_FEM_TEMPORARYLOCALFUNCTION_HH

#include <dune/fem/storage/array.hh>
#include <dune/fem/function/common/localfunction.hh>

namespace Dune
{

  /** \ingroup LocalFunction
   *  \class TemporaryLocalFunction
   *  \brief A temporary function carrying values for one entity
   *
   *  A TemporaryLocalFunction is a LocalFunction which is not associated with
   *  any DiscreteFunction. It can be used when generating discrete functions
   *  to temporarily store values for one entity.
   *
   *  \note Local DoF numbers correspond directly to array indices. Hence it
   *  may be more cache efficient to generate a TemporaryLocalFunction and then
   *  do only one update step on the discrete function's LocalFunction.
   *
   *  \param DiscreteFunctionSpaceImp type of the discrete function space, the
   *                                  local function shall belong to
   */
  template< class DiscreteFunctionSpaceImp,
            template< class > class ArrayAllocatorImp = DefaultArrayOverAllocator >
  class TemporaryLocalFunction
  : public LocalFunctionDefault
    < DiscreteFunctionSpaceImp,
      TemporaryLocalFunction< DiscreteFunctionSpaceImp, ArrayAllocatorImp > >
  {
  public:
    //! type of the discrete function space
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  private:
    typedef TemporaryLocalFunction< DiscreteFunctionSpaceType, ArrayAllocatorImp >
      ThisType;
    typedef LocalFunctionDefault< DiscreteFunctionSpaceType, ThisType > BaseType;

  public:
    //! type of domain vectors
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
    //! type of range vectors
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

    //! type of the jacobian
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;
    
    //! field type of domain vectors
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType;
    //! field type of range vectors
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;

    enum { DimDomain = DiscreteFunctionSpaceType :: DimDomain };
    enum { DimRange = DiscreteFunctionSpaceType :: DimRange };

    //! type of the base function set
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

  private:
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;
    typedef typename GridType :: template Codim< 0 > :: Entity EntityCodim0Type;

    typedef DynamicArray< RangeFieldType, ArrayAllocatorImp > DofArrayType;

  protected:
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;
    const EntityCodim0Type *entity_;

    BaseFunctionSetType baseFunctionSet_;

    DofArrayType dofs_;

  protected:
    using BaseType :: evaluate;
    using BaseType :: jacobian;

  public:
    /** \brief constructor creating a local function without binding it to an 
     *         entity
     *
     *  Creates the local function without initializing the fields depending on
     *  the current entity.
     *
     *  \note Before using the local function it must be initilized by
     *  \code
     *  localFunction.init( entity );
     *  \endcode
     *
     *  \param[in] dfSpace discrete function space the local function shall
     *                     belong to
     */
    inline explicit TemporaryLocalFunction
      ( const DiscreteFunctionSpaceType &dfSpace )
    : discreteFunctionSpace_( dfSpace ),
      entity_( 0 ),
      baseFunctionSet_(),
      dofs_()
    {
    }
    
    /** \brief constructor creating a local function and binding it to an
     *         entity
     *
     *  Creates the local function and initilizes the fields depending on the
     *  current entity. It is not necessary, though allowed, to call init
     *  before using the discrete function.
     *
     *  \note The degrees of freedom are not initialized by this function.
     *  
     *  \param[in] dfSpace discrete function space the local function shall
     *                     belong to
     *  \param[in] entity  entity for initialize the local function to
     */
    inline TemporaryLocalFunction ( const DiscreteFunctionSpaceType &dfSpace,
                                    const EntityCodim0Type &entity )
    : discreteFunctionSpace_( dfSpace ),
      entity_( &entity ),
      baseFunctionSet_( discreteFunctionSpace_.baseFunctionSet( entity ) ),
      dofs_( baseFunctionSet_.numBaseFunctions() )
    {
    }

    /** \brief copy constructor
     *
     *  Creates the local function as a copy of the specified one. It is bound
     *  to an entity if and only if the copied local function is bound to an
     *  entity.
     *
     *  \note The degrees of freedom are always copied, even if the copied
     *        local function is not bound to an entity.
     * 
     *  \param[in]  other  TemporaryLocalFunction to copy
     */
    inline TemporaryLocalFunction ( const ThisType &other )
    : discreteFunctionSpace_( other.discreteFunctionSpace_ ),
      entity_( other.entity_ ),
      baseFunctionSet_( other.baseFunctionSet_ ),
      dofs_( other.dofs_ )
    {
    }

    /** \brief access a local DoF
     *
     *  Grants read access to a local DoF.
     *
     *  \param[in] index local number of the DoF to access
     *
     *  \returns constant reference to the DoF
     */
    inline const RangeFieldType &operator[] ( const int index ) const
    {
      return dofs_[ index ];
    }

    /** \brief access a local DoF
     *
     *  Grants read and write access to a local DoF.
     *
     *  \param[in] index local number of the DoF to access
     *
     *  \returns reference to the DoF
     */
    inline RangeFieldType &operator[] ( const int index )
    {
      return dofs_[ index ];
    }

    /** obtain the base function set
     *
     *  The local function caches the base function set of the entity. This
     *  provides faster access than obtaining it from the space (which has to
     *  look up the base function set for the entity type first).
     *
     *  \returns constant reference to the base function set
     */
    inline const BaseFunctionSetType &baseFunctionSet () const
    {
      assert( entity_ != NULL );
      return baseFunctionSet_;
    }

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const DomainType &x,RangeType &ret) const */
    inline void evaluate ( const DomainType &x,
                           RangeType &ret ) const
    {
      assert( entity_ != NULL );
      
      ret = 0;
      const int numDofs = dofs_.size();
      for( int i = 0; i < numDofs; ++i )
      {
        RangeType phi;
        baseFunctionSet_.evaluate( i, x, phi );
        ret.axpy( dofs_[ i ], phi );
      }
    }

    /** \copydoc Dune::LocalFunctionInterface::evaluate(const QuadratureType &quadrature,const int quadPoint,RangeType &ret) const */
    template< class QuadratureType >
    inline void evaluate ( const QuadratureType &quadrature,
                           const int quadPoint,
                           RangeType &ret ) const
    {
      assert( entity_ != NULL );
      
      ret = 0;
      const int numDofs = dofs_.size();
      for( int i = 0; i < numDofs; ++i )
      {
        RangeType phi;
        baseFunctionSet_.evaluate( i, quadrature, quadPoint, phi );
        ret.axpy( dofs_[ i ], phi );
      }
    }

    /** \brief initialize the local function for an entity
     *
     *  Binds the local function to an entity.
     *
     *  \note A local function must be initialized to an entity before it can
     *        be used.
     *        
     *  \note This function can be called multiple times to use the local
     *        function for more than one entity.
     *
     *  \param[in] entity entity to bind the local function to
     */
    inline void init ( const EntityCodim0Type &entity )
    {
      entity_ = &entity;
      baseFunctionSet_ = discreteFunctionSpace_.baseFunctionSet( entity );
      dofs_.resize( baseFunctionSet_.numBaseFunctions() );
    }

    /** \copydoc Dune::LocalFunctionInterface::jacobian(const DomainType &x,JacobianRangeType &ret) const */
    inline void jacobian ( const DomainType &x, JacobianRangeType &ret ) const
    {
      assert( entity_ != NULL );
      
      typedef typename EntityCodim0Type :: Geometry GeometryType;
      typedef FieldMatrix< typename GeometryType :: ctype,
                           GeometryType :: mydimension,
                           GeometryType :: mydimension > GeometryJacobianType;
      
      const GeometryType &geometry = entity_->geometry();
      const GeometryJacobianType &inv = geometry.jacobianInverseTransposed( x );
      
      ret = 0;
      const int numDofs = dofs_.size();
      for( int i = 0; i < numDofs; ++i )
      {
        JacobianRangeType tmp;
        baseFunctionSet_.jacobian( i, x, tmp );
        ret.axpy( dofs_[ i ], tmp );
      }

      for( int i = 0; i < DimRange; ++i )
        ret[ i ] = FMatrixHelp :: mult( inv, ret[ i ] );
    }

    /** \copydoc Dune::LocalFunctionInterface::jacobian(const QuadratureType &quadrature,const int quadPoint,JacobianRangeType &ret) const */
    template< class QuadratureType >
    inline void jacobian ( const QuadratureType &quadrature,
                           const int quadPoint,
                           JacobianRangeType &ret ) const
    {
      assert( entity_ != NULL );
      
      typedef typename EntityCodim0Type :: Geometry GeometryType;
      typedef FieldMatrix< typename GeometryType :: ctype,
                           GeometryType :: mydimension,
                           GeometryType :: mydimension > GeometryJacobianType;
      
      const GeometryType &geometry = entity_->geometry();
      const GeometryJacobianType &inv
        = geometry.jacobianInverseTransposed( quadrature.point( quadPoint ) );
      
      ret = 0;
      const int numDofs = dofs_.size();
      for( int i = 0; i < numDofs; ++i )
      {
        JacobianRangeType tmp;
        baseFunctionSet_.jacobian( i, quadrature, quadPoint, tmp );
        ret.axpy( dofs_[ i ], tmp );
      }

      for( int i = 0; i < DimRange; ++i )
        ret[ i ] = FMatrixHelp :: mult( inv, ret[ i ] );
    }


    /** \copydoc Dune::LocalFunctionInterface::numDofs */
    inline int numDofs () const
    {
      return dofs_.size();
    }
  };


  
  template< class DiscreteFunctionSpaceImp,
            template< class > class ArrayAllocatorImp = DefaultArrayAllocator >
  class TemporaryLocalFunctionFactory
  {
  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  private:
    typedef TemporaryLocalFunctionFactory
      < DiscreteFunctionSpaceType, ArrayAllocatorImp >
      ThisType;

  public:
    typedef TemporaryLocalFunction< DiscreteFunctionSpaceType, ArrayAllocatorImp >
      ObjectType;

  protected:
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;

  public:
    inline explicit
    TemporaryLocalFunctionFactory ( const DiscreteFunctionSpaceType &dfSpace )
    : discreteFunctionSpace_( dfSpace )
    {
    }

    inline ObjectType *newObject () const
    {
      return new ObjectType( discreteFunctionSpace_ );
    }
  };

}

#endif
