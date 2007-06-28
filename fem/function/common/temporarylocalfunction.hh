#ifndef DUNE_FEM_TEMPORARYLOCALFUNCTION_HH
#define DUNE_FEM_TEMPORARYLOCALFUNCTION_HH

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
  template< class DiscreteFunctionSpaceImp >
  class TemporaryLocalFunction
  : public LocalFunctionDefault
    < DiscreteFunctionSpaceImp,
      TemporaryLocalFunction< DiscreteFunctionSpaceImp >
    >
  {
  public:
    //! type of the discrete function space
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  private:
    typedef TemporaryLocalFunction< DiscreteFunctionSpaceType > ThisType;
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

    //! type of the base function set
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

  private:
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename GridPartType :: EntityCodim0Type EntityCodim0Type;

  protected:
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;
    const EntityCodim0Type *entity_;

    const BaseFunctionSetType *baseFunctionSet_;
    int numBaseFunctions_;

    RangeFieldType *dofs_;

  public:
    /** \brief constructor creating a local function without binding it to an 
     *         entity
     *
     *  Creates the local function without initializing the fields depending on
     *  the current entity.
     *
     *  \note Before using the local function it must be initilized by
     *        \code
     *          localFunction.init( entity );
     *        \endcode
     *
     *  \param[in] dfSpace discrete function space the local function shall
     *                     belong to
     */
    inline TemporaryLocalFunction ( const DiscreteFunctionSpaceType &dfSpace )
    : discreteFunctionSpace_( dfSpace ),
      entity_( NULL ),
      baseFunctionSet_( NULL ),
      numBaseFunctions_( 0 ),
      dofs_( NULL );
    {
    }
    
     /** \brief constructor creating a local function without binding it to an 
     *         entity
     *
     *  Creates the local function without initializing the fields depending on
     *  the current entity.
     *  
     *  \param[in] dfSpace discrete function space the local function shall
     *                     belong to
     *  \param[in] entity  entity for initialize the local function to
     */
    inline TemporaryLocalFunction ( const DiscreteFunctionSpaceType &dfSpace,
                                    const EntityCodim0Type &entity )
    : discreteFunctionSpace_( dfSpace ),
      entity_( &entity ),
      baseFunctionSet_( &(discreteFunctionSpace_( entity )) ),
      numBaseFunctions_( baseFunctionSet_->numBaseFunctions() ),
      dofs_( new RangeFieldType[ numBaseFunctions_ ] )
    {
      assert( dofs_ != NULL );
    }

    inline ~TemporaryLocalFunction ()
    {
      if( dofs_ != NULL )
        delete[] dofs_;
    }

    /** \brief access a local DoF
     *
     *  Grants read access to a local DoF.
     *
     *  \param[in] index local number of the DoF to access
     *
     *  \returns constant reference to the DoF
     */
    inline const RangeFieldType &operator[] ( int index ) const
    {
      assert( (index >= 0) && (index < numBaseFunctions_) );
      return fields_[ index ];
    }

    /** \brief access a local DoF
     *
     *  Grants read and write access to a local DoF.
     *
     *  \param[in] index local number of the DoF to access
     *
     *  \returns reference to the DoF
     */
    inline RangeFieldType &operator[] ( int index )
    {
      assert( (index >= 0) && (index < numBaseFunctions_) );
      return fields_[ index ];
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
      assert( entity != NULL );
      return *baseFunctionSet_;
    }

    /** evaluate the function in local coordinate x
     *
     *  Evaluate the local function in a point x, given in local coordinates
     *  (i.e. coordinates within the reference element).
     *
     *  \param[in]  x   local coordinate to evaluate the function in
     *  \param[out] phi variable to hold the function value
     */
    inline void evaluate ( const DomainType &x, RangeType &phi )
    {
      assert( entity != NULL );
      
      phi = 0;
      for( int i = 0; i < numBaseFunctions_; ++i )
      {
        RangeType psi;
        baseFunctionSet_.evaluate( i, x, psi );
        phi.axpy( fields_[ i ], psi );
      }
    }

    /** evaluate the function in a quadrature point
     *
     *  Evaluate the local function in a quadrature point. If the discrete
     *  function space uses a caching storage, this is more efficient than
     *  \code
     *    localFunction.evaluate( quadrature.point( point ), phi );
     *  \endcode
     *
     *  \param[in]  quadrature quadrature to use
     *  \param[in]  point      index of the quadrature point to evaluate the
     *                         function in
     *  \param[out] phi        variable to hold the function value
     */
    template< class QuadratureType >
    inline void evaluate ( QuadratureType &quadrature, int point, RangeType &phi )
    {
      assert( entity != NULL );
      
      phi = 0;
      for( int i = 0; i < numBaseFunctions_; ++i )
      {
        RangeType psi;
        baseFunctionSet_.evaluate( i, quadrature, point, psi );
        phi.axpy( fields_[ i ], psi );
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
      const int oldNumDofs = numBaseFunctions_;

      entity_ = &entity;
      baseFunctionSet_ = &(discreteFunctionSpace_( entity ));
      numBaseFunctions_( baseFunctionSet_->numBaseFunctions() );
      if( numBaseFunctions_ != oldNumDofs ) {
        if( dofs_ != NULL )
          delete[] dofs_;
        dofs_ = new RangeFieldType[ numBaseFunctions_ ];
      }
      assert( dofs_ != NULL );
    }

    /** evaluate Jacobian in local coordinate x
     *
     *  Evaluate the Jacobian of the local function in a point x, given in
     *  local coordinates (i.e. coordinates within the reference element).
     *
     *  \param[in]  x    local coordinate to evaluate the Jacobian in
     *  \param[out] grad variable to hold the Jacobian
     */
    inline void jacobian ( const DomainType &x, JacobianRangeType &grad )
    {
      assert( entity != NULL );
      
      typedef typename EntityCodim0Type :: Geometry GeometryType;
      typedef FieldMatrix< typename GeometryType :: ctype,
                           GeometryType :: mydimension,
                           GeometryType :: mydimension > GeometryJacobianType;
      
      const GeometryType &geometry = entity->geometry();
      const GeometryJacobianType &inv = geometry.jacobianInverseTransposed( x );
      
      grad = 0;
      for( int i = 0; i < numBaseFunctions_; ++i )
      {
        JacobianRangeType tmp;
        baseFunctionSet_.jacobian( i, x, tmp );
        grad.axpy( fields_[ i ], psi );
      }

      for( int i = 0; i < DimRange; ++i )
        grad[ i ] = FMatrixHelp :: mult( inv, grad[ i ] );
    }

    /** evaluate Jacobian in a quadrature point
     *
     *  Evaluate the Jacobian of the local function in a quadrature point. If
     *  the discrete function space uses a caching storage, this is more
     *  efficient than
     *  \code
     *    localFunction.jacobian( quadrature.point( point ), phi );
     *  \endcode
     *
     *  \param[in]  quadrature quadrature to use
     *  \param[in]  point      index of the quadrature point to evaluate the
     *                         Jacobian in
     *  \param[out] phi        variable to hold the Jacobian
     */
    template< class QuadratureType >
    inline void jacobian ( QuadratureType &quadrature, int point, JacobianRangeType &grad )
    {
      assert( entity != NULL );
      
      typedef typename EntityCodim0Type :: Geometry GeometryType;
      typedef FieldMatrix< typename GeometryType :: ctype,
                           GeometryType :: mydimension,
                           GeometryType :: mydimension > GeometryJacobianType;
      
      const GeometryType &geometry = entity->geometry();
      const GeometryJacobianType &inv = geometry.jacobianInverseTransposed( x );
      
      grad = 0;
      for( int i = 0; i < numBaseFunctions_; ++i )
      {
        JacobianRangeType tmp;
        baseFunctionSet_.jacobian( i, quadrature, point, tmp );
        grad.axpy( fields_[ i ], psi );
      }

      for( int i = 0; i < DimRange; ++i )
        grad[ i ] = FMatrixHelp :: mult( inv, grad[ i ] );
    }


    /** obtain the number of local DoFs
     *
     *  Obtain the number of local DoFs of this local function. The value is
     *  identical to the number of base functions on the entity.
     *
     *  \returns number of local DoFs
     */
    inline int numDofs () const
    {
      return numBaseFunctions_;
    }
  };
  
}

#endif
