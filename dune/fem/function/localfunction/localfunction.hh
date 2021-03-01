#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_LOCALFUNCTION_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_LOCALFUNCTION_HH

#include <utility>
#include <type_traits>
#include <tuple>

#include <dune/common/deprecated.hh>
#include <dune/common/fvector.hh>


namespace Dune
{
  namespace Fem {
    // forward declaration for DenseMatVecTraits
    template< class BasisFunctionSet, class LocalDofVector >
    class LocalFunction;
  }


  // DenseMatVecTraits for LocalFunction for Discrete Functions
  // ----------------------------------------------------------

  template< class BasisFunctionSet, class LocalDofVector >
  struct DenseMatVecTraits< Fem::LocalFunction< BasisFunctionSet, LocalDofVector > >
  {
    typedef Fem::LocalFunction< BasisFunctionSet, LocalDofVector >  derived_type;
    typedef typename LocalDofVector::value_type                     value_type;
    //typedef LocalDofVector                                          container_type;
    typedef typename LocalDofVector::size_type                      size_type;
  };



  // FieldTraits for LocalContribution for Discrete Functions
  // --------------------------------------------------------

  template< class BasisFunctionSet, class LocalDofVector >
  struct FieldTraits< Fem::LocalFunction< BasisFunctionSet, LocalDofVector > >
  {
    typedef typename FieldTraits< typename LocalDofVector::value_type >::field_type field_type;
    typedef typename FieldTraits< typename LocalDofVector::value_type >::real_type  real_type;
  };


  namespace Fem
  {

    /** \addtogroup LocalFunction
     *
     * On every element from a discrete function the local funtion can be
     * accessed. With the local function one has access to the dof and on the
     * other hand to the basis function set of this actual element. Therefore
     * this is called a local function.
     *
     * \remarks The interface for using a LocalFunction is defined by the class
     *          LocalFunction.
     *
     * \{
     */



    // LocalFunction
    // -------------

    /** \class LocalFunction
     *  \brief interface for local functions
     *
     *  Local functions are used to represend a discrete function on one entity.
     *  The LocalFunctionInterface defines the functionality that can be expected
     *  from such a local function.
     */
    template< class BasisFunctionSet, class LocalDofVector >
    class LocalFunction
      : public Dune::DenseVector< LocalFunction< BasisFunctionSet, LocalDofVector > >
    {
      typedef LocalFunction< BasisFunctionSet, LocalDofVector > ThisType;
      typedef Dune::DenseVector< ThisType > BaseType;
    public:

      //! type of basis function set
      typedef BasisFunctionSet BasisFunctionSetType;

      /** \brief  type of local Dof Vector  */
      typedef LocalDofVector LocalDofVectorType;

      //! type of DoF use with the discrete function
      typedef typename LocalDofVectorType::value_type DofType;

      //! type of index
      typedef typename LocalDofVectorType::size_type SizeType;

      //! type of the entity, the local function lives on is given by the space
      typedef typename BasisFunctionSetType::EntityType EntityType;

      //! type of functionspace
      typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;

      //! field type of the domain
      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
      //! field type of the range
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      //! type of domain vectors, i.e., type of coordinates
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! type of range vectors, i.e., type of function values
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! type of the Jacobian, i.e., type of evaluated Jacobian matrix
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! type of the Hessian
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      //! type of local coordinates
      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;

      //! dimension of the domain
      static const int dimDomain = FunctionSpaceType::dimDomain;
      //! dimension of the range
      static const int dimRange = FunctionSpaceType::dimRange;

      //! default constructor, calls default ctor of BasisFunctionSetType and LocalDofVectorType
      LocalFunction () {}

      //! ctor taking a basisFunctionSet, calling default ctor for LocalDofVectorType, and resize
      explicit LocalFunction ( const BasisFunctionSetType &basisFunctionSet )
      : basisFunctionSet_( basisFunctionSet )
      {
        localDofVector_.resize( basisFunctionSet.size() );
      }

      //! ctor taking a localDofVector, calling default ctor for BasisFunctionSetType
      explicit LocalFunction ( const LocalDofVectorType &localDofVector ) : localDofVector_( localDofVector ) {}

      //! copy given agruments
      LocalFunction ( const BasisFunctionSetType &basisFunctionSet, const LocalDofVector &localDofVector )
      : basisFunctionSet_( basisFunctionSet ),
        localDofVector_( localDofVector )
      {
        localDofVector_.resize( basisFunctionSet.size() );
      }

      //! half move ctor
      explicit LocalFunction ( LocalDofVectorType &&localDofVector ) : localDofVector_( localDofVector ) {}

      //! half move ctor
      LocalFunction ( const BasisFunctionSetType &basisFunctionSet, LocalDofVector &&localDofVector )
      : basisFunctionSet_( basisFunctionSet ),
        localDofVector_( localDofVector )
      {
        localDofVector_.resize( basisFunctionSet.size() );
      }

      //! move constructor
      LocalFunction ( ThisType && other )
      : basisFunctionSet_( std::move( other.basisFunctionSet_ ) ),
        localDofVector_( std::move( other.localDofVector_ ) )
      {}

      //! copy constructor
      LocalFunction ( const ThisType & other )
      : basisFunctionSet_( other.basisFunctionSet_ ),
        localDofVector_( other.localDofVector_ )
      {}

      /** \brief access to local dofs (read-only)
       *
       *  \param[in]  num  local dof number
       *  \return reference to dof
       */
      const DofType &operator[] ( SizeType num ) const { return localDofVector_[ num ]; }

      /** \brief access to local dofs (read-write)
       *
       *  \param[in]  num  local DoF number
       *  \return reference to DoF
       */
      DofType &operator[] ( SizeType num ) { return localDofVector_[ num ]; }

      using BaseType :: operator +=;
      using BaseType :: operator -=;
      using BaseType :: operator *=;
      using BaseType :: operator /=;

      /** \brief assign all DoFs of this local function
       *
       *  \param[in]  lf  local function to assign DoFs from
       */
      template< class T >
      void assign ( const LocalFunction< BasisFunctionSet, T > &other )
      {
        localDofVector() = other.localDofVector();
      }

      /** \brief set all DoFs to zero */
      void clear ()
      {
        std::fill( localDofVector().begin(), localDofVector().end(), DofType( 0 ) );
      }

#if 0
      /** \brief add a multiple of another local function to this one
       *
       *  \note The local function to add may be any implementation of a
       *        LocalFunction.
       *
       *  \param[in]  s   scalar factor to scale lf with
       *  \param[in]  lf  local function to add
       *
       *  \returns a reference to this local function (i.e., *this)
       */
      template< class T >
      ThisType &axpy ( const RangeFieldType s, const LocalFunction< BasisFunctionSet, T > &other )
      {
        localDofVector().axpy( s, other.localDofVector() );
        return *this;
      }
#endif
      using BaseType :: axpy;

      /** \brief axpy operation for local function
       *
       *  Denoting the DoFs of the local function by \f$u_i\f$ and the basis
       *  functions by \f$\varphi_i\f$, this function performs the following
       *  operation:
       *  \f[
       *  u_i = u_i + factor \cdot \varphi_i( x )
       *  \f]
       *
       *  \param[in]  x       point to evaluate basis functions in
       *  \param[in]  factor  axpy factor
       */
      template< class PointType >
      void axpy ( const PointType &x, const RangeType &factor )
      {
        basisFunctionSet().axpy( x, factor, localDofVector() );
      }

      /** \brief axpy operation for local function
       *
       *  Denoting the DoFs of the local function by \f$u_i\f$ and the basis
       *  functions by \f$\varphi_i\f$, this function performs the following
       *  operation:
       *  \f[
       *  u_i = u_i + factor \cdot \nabla\varphi_i( x )
       *  \f]
       *
       *  \param[in]  x       point to evaluate jacobian of basis functions in
       *  \param[in]  factor  axpy factor
       */
      template< class PointType >
      void axpy ( const PointType &x, const JacobianRangeType &factor)
      {
        basisFunctionSet().axpy( x, factor, localDofVector() );
      }
      template< class PointType >
      void axpy ( const PointType &x, const HessianRangeType &factor)
      {
        basisFunctionSet().axpy( x, factor, localDofVector() );
      }

      /** \brief axpy operation for local function
       *
       *  Denoting the DoFs of the local function by \f$u_i\f$ and the basis
       *  functions by \f$\varphi_i\f$, this function performs the following
       *  operation:
       *  \f[
       *  u_i = u_i + factor1 \cdot \varphi_i( x ) + factor2 \cdot \nabla\varphi_i( x )
       *  \f]
       *
       *  \param[in]  x        point to evaluate basis functions in
       *  \param[in]  factor1  axpy factor for \f$\varphi( x )\f$
       *  \param[in]  factor2  axpy factor for \f$\nabla\varphi( x )\f$
       */
      template< class PointType >
      void axpy ( const PointType &x, const RangeType &factor1, const JacobianRangeType &factor2 )
      {
        basisFunctionSet().axpy( x, factor1, factor2, localDofVector() );
      }

      /** \brief obtain the order of this local function
       *
       *  The order of a local function refers to the polynomial order required
       *  to integrate it exactly.
       *
       *  \note It is not completely clear what this value should be, e.g., for
       *        bilinear basis functions.
       *
       *  \returns order of the local function
       */
      int order () const { return basisFunctionSet().order(); }

      /** \brief obtain the basis function set for this local function
       *
       *  \returns reference to the basis function set
       */
      const BasisFunctionSetType &basisFunctionSet () const { return basisFunctionSet_; }

      /** \brief obtain the entity, this local function lives on
       *
       *  \returns reference to the entity
       */
      const EntityType &entity () const { return basisFunctionSet().entity(); }


      /** \brief evaluate the local function
       *
       *  \param[in]   x    evaluation point in local coordinates
       *  \param[out]  ret  value of the function in the given point
       */
      template< class PointType >
      void evaluate ( const PointType &x, RangeType &ret ) const
      {
        basisFunctionSet().evaluateAll( x, localDofVector(), ret );
      }

      /** \brief evaluate Jacobian of the local function
       *
       *  \note Though the Jacobian is evaluated on the reference element, the
       *        return value is the Jacobian with respect to the actual entity.
       *
       *  \param[in]   x    evaluation point in local coordinates
       *  \param[out]  ret  Jacobian of the function in the evaluation point
       */
      template< class PointType >
      void jacobian ( const PointType &x, JacobianRangeType &ret ) const
      {
        basisFunctionSet().jacobianAll( x, localDofVector(), ret );
      }

      /** \brief evaluate Hessian of the local function
       *
       *  \note Though the Hessian is evaluated on the reference element, the
       *        return value is the Hessian with respect to the actual entity.
       *
       *  \param[in]   x        evaluation point in local coordinates
       *  \param[out]  ret  Hessian of the function in the evaluation point
       */
      template< class PointType >
      void hessian ( const PointType &x, HessianRangeType &ret ) const
      {
        basisFunctionSet().hessianAll( x, localDofVector(), ret );
      }

      /** \brief obtain the number of local DoFs
       *
       *  Obtain the number of local DoFs of this local function. The value is
       *  identical to the number of basis functons on the entity.
       *
       *  \returns number of local DoFs
       */
      int numDofs () const { return localDofVector().size(); }

      /** \brief obtain the number of local DoFs
       *
       *  Obtain the number of local DoFs of this local function. The value is
       *  identical to the number of basis functons on the entity.
       *
       *  \returns number of local DoFs
       */
      SizeType size () const { return localDofVector().size(); }

      /** \brief evaluate all basisfunctions for all quadrature points, multiply with the given factor and
                 add the result to the local coefficients  */
      template< class QuadratureType, class ... Vectors >
      void axpyQuadrature ( const QuadratureType &quad, const Vectors& ... values )
      {
        static_assert( sizeof...( Vectors ) > 0, "evaluateQuadrature needs to be called with at least one vector." );
        std::ignore =
          std::make_tuple( ( basisFunctionSet().axpy( quad, values, localDofVector() ), 1 )...);
      }

      /** \brief evaluate all basisfunctions for all quadrature points, multiply with the given factor and
                 add the result to the local coefficients */
      template< class QuadratureType, class RangeVectorType, class JacobianRangeVectorType >
      void axpyQuadrature ( const QuadratureType &quad,
                            const RangeVectorType& rangeVector,
                            const JacobianRangeVectorType& jacobianVector )
      {
        basisFunctionSet().axpy( quad, rangeVector, jacobianVector, localDofVector() );
      }

      /** \brief evaluate all basisfunctions for all quadrature points and store the results in the result vector */
      template< class QuadratureType, class ... Vectors >
      void evaluateQuadrature( const QuadratureType &quad, Vectors& ... vec ) const
      {
        static_assert( sizeof...( Vectors ) > 0, "evaluateQuadrature needs to be called with at least one vector." );
        std::ignore =
          std::make_tuple( ( evaluateQuadrature( quad, vec, typename std::decay< decltype( vec[ 0 ] ) >::type() ), 1 )
            ... );
      }

      /** \brief evaluate all Jacobians for all basis functions for all quadrature points and
       *         store the results in the result vector */
      template< class QuadratureType, class ... Vectors >
      void jacobianQuadrature( const QuadratureType &quad, Vectors& ... vec ) const
      {
        static_assert( sizeof...( Vectors ) > 0, "evaluateQuadrature needs to be called with at least one vector." );
        std::ignore =
          std::make_tuple( ( evaluateQuadrature( quad, vec, typename std::decay< decltype( vec[ 0 ] ) >::type() ), 1 )
            ... );
      }

      /** \brief evaluate all hessians of all basis functions for all quadrature points and store the results in the result vector */
      template< class QuadratureType, class ... Vectors >
      void hessianQuadrature( const QuadratureType &quad, Vectors& ... vec ) const
      {
        // make sure vec size if large enough
        static_assert( sizeof...( Vectors ) > 0, "evaluateQuadrature needs to be called with at least one vector." );
        std::ignore =
           std::make_tuple( ( evaluateQuadrature( quad, vec, typename std::decay< decltype( vec[ 0 ] ) >::type() ), 1 )
           ... );
      }

      /** \brief return const reference to local Dof Vector  */
      const LocalDofVectorType &localDofVector () const { return localDofVector_; }

      /** \brief return mutable reference to local Dof Vector  */
      LocalDofVectorType &localDofVector () { return localDofVector_; }


      /** \brief Returns true if local function if bind or init was previously called.
       */
      bool valid () const
      {
        return basisFunctionSet_.valid();
      }

    protected:
      /** \brief initialize the local function for an entity
       *
          Binds the local function to an basisFunctionSet and entity.

          \note Must be overloaded on the derived implementation class.

          \param[in] entity to bind the local function to
       */
      void init ( const EntityType& entity )
      {
        DUNE_THROW(NotImplemented,"LocalFunction::init( entity ) must be overloaded on derived class!");
      }

      /** \brief initialize the local function for an entity
       *
          Binds the local function to an basisFunctionSet and entity.

          \note Must be overloaded on the derived implementation class.

          \param[in] entity to bind the local function to
       */
      void bind ( const EntityType& entity )
      {
        DUNE_THROW(NotImplemented,"LocalFunction::bind( entity ) must be overloaded on derived class!");
      }

      /** \brief clears the local function by removing the basisFunctionSet
       *
       */
      void unbind ()
      {
        // basically sets entity pointer inside basis function set to nullptr
        basisFunctionSet_ = BasisFunctionSetType();
      }

      /** \brief initialize the local function for an basisFunctionSet
       *
          Binds the local function to an basisFunctionSet and entity.

          \note A local function must be initialized to an entity before it can
                be used.

          \note This function can be called multiple times to use the local
                function for more than one entity.

          \param[in] basisFunctionSet to bind the local function to
       */
      void init ( const BasisFunctionSetType &basisFunctionSet )
      {
        basisFunctionSet_ = basisFunctionSet;
        localDofVector_.resize( basisFunctionSet.size() );
      }

      // evaluate local function and store results in vector of RangeTypes
      // this method only helps to identify the correct method on
      // the basis function set
      template< class QuadratureType, class VectorType >
      void evaluateQuadrature( const QuadratureType &quad, VectorType &result, const RangeType & ) const
      {
        basisFunctionSet().evaluateAll( quad, localDofVector(), result );
      }

      // evaluate jacobian of local function and store result in vector of
      // JacobianRangeTypes, this method only helps to identify the correct method on
      // the basis function set
      template< class QuadratureType, class VectorType >
      void evaluateQuadrature( const QuadratureType &quad, VectorType &result, const JacobianRangeType & ) const
      {
        basisFunctionSet().jacobianAll( quad, localDofVector(), result );
      }

      // evaluate jacobian of local function and store result in vector of
      // JacobianRangeTypes, this method only helps to identify the correct method on
      // the basis function set
      template< class QuadratureType, class VectorType >
      void evaluateQuadrature( const QuadratureType &quad, VectorType &result, const HessianRangeType & ) const
      {
        basisFunctionSet().hessianAll( quad, localDofVector(), result );
      }

      BasisFunctionSetType basisFunctionSet_;
      LocalDofVectorType localDofVector_;
    };

    /** \} */

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_LOCALFUNCTION_HH
