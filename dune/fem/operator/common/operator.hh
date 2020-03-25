#ifndef DUNE_FEM_OPERATOR_HH
#define DUNE_FEM_OPERATOR_HH

#include <cassert>
#include <typeindex>
#include <cassert>

#include <dune/common/exceptions.hh>

#include <dune/fem/version.hh>
#include <dune/common/exceptions.hh>
#include <dune/fem/operator/common/mapping.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class Operator
     *  \brief abstract operator
     *
     *  Operators map a discrete function onto another discrete function.
     *  Their interface is described by the abstract class Operator.
     *
     *  \tparam  DomainFunction  type of discrete function for the domain
     *  \tparam  RangeFunction   type of discrete function for the range
     *                           (defaults to DomainFunction)
     *
     *  \interfaceclass
     */
    template< class DomainFunction, class RangeFunction = DomainFunction >
    struct Operator
    {
      /** \brief type of discrete function in the operator's domain */
      typedef DomainFunction DomainFunctionType;
      /** \brief type of discrete function in the operator's range */
      typedef RangeFunction RangeFunctionType;

      /** \brief field type of the operator's domain */
      typedef typename DomainFunction::RangeFieldType DomainFieldType;
      /** \brief field type of the operator's range */
      typedef typename RangeFunction::RangeFieldType RangeFieldType;

      virtual ~Operator ()
      {}

      /** \brief application operator
       *
       *  \param[in]   u  argument discrete function
       *  \param[out]  w  destination discrete function
       *
       *  \note This method has to be implemented by all derived classes.
       */
      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const = 0;

      /** \brief finalization of operator
       *
       *  \note The default implementation is empty.
       */
      virtual void finalize () {}
    };

    /** \class LinearOperator
     *  \brief abstract affine-linear operator
     *
     *  Operators map a discrete function onto another discrete
     *  function. Their interface is described by the abstract class
     *  Operator. Implementation should derive from LinearOperator to
     *  indicate that they model an affine linear operator of the form
     *
     *  @f[
     *  u\mapsto A\,u + b
     *  @f]
     *
     *  with a linear Operator @f$A@f$ and an affine translation @f$b@f$.
     *
     *  \tparam  DomainFunction  type of discrete function for the domain
     *  \tparam  RangeFunction   type of discrete function for the range
     *                           (defaults to DomainFunction)
     *
     *  \interfaceclass
     */
    template< class DomainFunction, class RangeFunction = DomainFunction >
    struct LinearOperator
      : public virtual Operator<DomainFunction, RangeFunction>
    {
      /**Return @c true if the Operator is symmetric. */
      virtual bool symmetric() const
      {
	      return false;
      }
      /**Return @c true if the Operator is positive definite. */
      virtual bool positiveDefinite() const
      {
	      return false;
      }
    };

    /** \class AssembledOperator
     *  \brief abstract matrix operator
     *
     *  Operators map a discrete function onto another discrete
     *  function. Their interface is described by the abstract class
     *  Operator. Implementation should derive from AssembledOperator to
     *  indicate that they model an affine linear operator of the form
     *
     *  @f[
     *  u\mapsto A\,u
     *  @f]
     *
     *  with a matrix @f$A@f$. Jacobians of LinearOperator classes,
     *  for instance, could be modelled as matrices.
     *
     *  \tparam  DomainFunction  type of discrete function for the domain
     *  \tparam  RangeFunction   type of discrete function for the range
     *                           (defaults to DomainFunction)
     *
     *  \interfaceclass
     */
    template< class DomainFunction, class RangeFunction = DomainFunction >
    class AssembledOperator
      : public virtual LinearOperator<DomainFunction, RangeFunction>
    {
    public:
      /** \brief commit intermediate states of linear operator assembly */
      virtual void flushAssembly() {}

      /** \brief Initiate the assemble of values using the LocalContribution concept
       *  \tparam AssembleOperation the specific operation (Add, Set, ...)
       */
      template< class AssembleOperation >
      void beginAssemble ()
      {
        const std::type_index id( typeid( AssembleOperation ) );
        if( assembleOperation_ != id )
        {
          if( assembleOperation_ != std::type_index( typeid( void ) ) )
            DUNE_THROW( InvalidStateException, "Another assemble operation in progress" );
          assembleOperation_ = id;
          assert( assembleCount_ == 0 );
          AssembleOperation::begin( *this );
        }
        ++assembleCount_;
      }

      /** \brief Finalize the assemble of values using the LocalContribution concept
       *  \tparam AssembleOperation the specific operation (Add, Set, ...)
       */
      template< class AssembleOperation >
      void endAssemble ()
      {
        const std::type_index id( typeid( AssembleOperation ) );
        if( assembleOperation_ != id )
          DUNE_THROW( InvalidStateException, "Assemble operation not in progress" );
        assert( assembleCount_ > 0 );
        if( --assembleCount_ == 0 )
        {
          AssembleOperation::end( *this );
          assembleOperation_ = std::type_index( typeid( void ) );
        }
      }

    protected:
      std::type_index assembleOperation_ = std::type_index( typeid( void ) );
      std::size_t assembleCount_ = 0;
    };


    // Is*Operator
    // -----------

    namespace Impl
    {

      template< class T >
      using DomainFunctionType_t = typename T::DomainFunctionType;

      template< class T >
      using RangeFunctionType_t = typename T::RangeFunctionType;

      template< class T >
      using IsOperatorImpl = std::is_base_of< ::Dune::Fem::Operator< DomainFunctionType_t< T >, RangeFunctionType_t< T > >, T >;

      template< class T >
      using IsLinearOperatorImpl = std::is_base_of< ::Dune::Fem::LinearOperator< DomainFunctionType_t< T >, RangeFunctionType_t< T > >, T >;

      template< class T >
      using IsAssembledOperatorImpl = std::is_base_of< ::Dune::Fem::AssembledOperator< DomainFunctionType_t< T >, RangeFunctionType_t< T > >, T >;

    } // namespace Impl

    template< class T >
    using IsOperator = Impl::IsOperatorImpl< std::decay_t< T > >;

    template< class T >
    using IsLinearOperator = Impl::IsOperatorImpl< std::decay_t< T > >;

    template< class T >
    using IsAssembledOperator = Impl::IsAssembledOperatorImpl< std::decay_t< T > >;

  } // namespace Fem



  /** @addtogroup OperatorCommon
      Operators are mappings from function spaces into function spaces.

      \remarks
      The most general interface for a mapping is defined by Mapping. From
      Mapping the Function is derived and also Operator is derived.
      Operator defines the interface for operations on discrete functions
      while Function is an interface for functions like sinus.

      @{
   */

  /** \brief An abstract operator
   Interface class for Operators. Operators are applied to Functions and
   the result is a Function again.

   \interfaceclass
  */
  template <typename DFieldType, typename RFieldType,
            typename DType , typename RType>
  class Operator
  : public Fem::Mapping < DFieldType, RFieldType, DType, RType >,
    public virtual Fem::Operator< DType, RType >
  {
    typedef Fem::Operator< DType, RType > BaseType;

  protected:
    //! \brief type of mapping base class
    typedef Fem::Mapping <DFieldType,RFieldType,DType,RType> MappingType;

  public:
    //- remember template parameters for derived classes
    typedef DType DomainType;
    typedef RType  RangeType;
    typedef DFieldType DomainFieldType;
    typedef RFieldType RangeFieldType;

    using BaseType::operator();
    using BaseType::finalize;

  protected:
    /** \brief The method apply calls the application operator. The method
        has to be implemented here, because this method called when a mapping list
        is evaluated.
        \param[in] arg argument
        \param[out] dest destination
    */
    virtual void apply (const DomainType& arg, RangeType& dest) const
    {
      this->operator() (arg, dest);
    }
  }; // class Operator

  ///@}

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_HH
