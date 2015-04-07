#ifndef DUNE_FEM_FUNCTION_COMMON_INSTATIONARY_HH
#define DUNE_FEM_FUNCTION_COMMON_INSTATIONARY_HH

#include <functional>
#include <type_traits>
#include <utility>

#include <dune/fem/function/common/function.hh>

namespace Dune
{

  namespace Fem
  {

    // BasicInstationaryFunction
    // -------------------------

    /** \brief basic wrapper class (still a CRTP) for instationary functions
     *
     *  This class implements two methods
     *  \code
     *    double setTime(double);
     *    double time() const;
     *  \endcode
     *  for wrapping some time dependent function and making it a
     *  Dune::Fem::Function that may only depend on spatial variables.
     *
     *  \ingroup Functions
     */
    template< class FunctionSpace, class Function >
    class BasicInstationaryFunction
      : public Dune::Fem::Function< FunctionSpace, Function >
    {
    public:
      /** \name Constructon
       *  \{
       */

      explicit BasicInstationaryFunction ( double time )
        : time_( time )
      {}

      /** \} */

      /** \name Set time
       *  \{
       */

      /** \brief set time to give value
       *
       *  \param[in]  time  time to be used
       *
       *  \returns set time
       */
      double setTime ( double time ) { return (time_ = time); }

      /** \brief return set time */
      double time () const { return time_; }

      /** \} */

    private:
      double time_;
    };



#ifndef DOXYGEN

    namespace __InstationaryFunction
    {

      // HoldCopy
      // --------

      template< class Function >
      struct HoldCopy
      {
        explicit HoldCopy ( Function function )
          : function_( std::move( function ) )
        {}

        const Function &get () const noexcept { return function_; }

      private:
        Function function_;
      };



      // HoldReference
      // -------------

      template< class Function >
      struct HoldReference
      {
        explicit HoldReference ( const Function &function )
          : function_( function )
        {}

        const Function &get () const noexcept { return function_.get(); }

      private:
        std::reference_wrapper< const Function > function_;
      };

    } // namespace __InstationaryFunction

#endif // #ifndef DOXYGEN



    // InstationaryFunction
    // --------------------

    /** \brief implementation of a Dune::Fem::Function taking an instationary
     *         function
     *
     *  It is assumed that all evaluation methods are present on the parameter
     *  function and have a second parameter for the time:
     *  \code
     *    void Function::evaluate(const DomainType &x, double time, RangeType &value) const;
     *    void Function::jacobian(const DomainType &x, double time, JacobianRangeType &value) const;
     *    void Function::hessian(const DomainType &x, double time, HessianRangeType &value) const;
     *  \endcode
     *
     *  Users may prescribe how the parameter function is stored by providing a
     *  second template parameter, the storage policy. The policy is class that
     *  must be constructible from a function object or reference and that has
     *  a single method:
     *  \code
     *    const Function &Policy::get() const;
     *  \endcode
     *  The default policy is to copy the function parameter. The free-standing
     *  method
     *  \code
     *    Dune::Fem::instationaryFunction
     *  \endcode
     *  may be used to conveniently create a new instance of InstationaryFunction. Use
     *  \code
     *    auto g = instationaryFunction( std::cref( f ), 0. );
     *  \endcode
     *  to create an instationary function that holds a reference to \c f
     *  instead of a copy.
     *
     *  \tparam  Function  an instationary function
     *  \tparam  StoragePolicy  storage policy
     *
     *  \ingroup Functions
     */
    template< class Function,
              template< class > class StoragePolicy = __InstationaryFunction::HoldCopy >
    class InstationaryFunction
      : public BasicInstationaryFunction< typename Function::FunctionSpaceType, InstationaryFunction< Function, StoragePolicy > >,
        private StoragePolicy< Function >
    {
      typedef BasicInstationaryFunction< typename Function::FunctionSpaceType, InstationaryFunction< Function, StoragePolicy > > BaseType;
      typedef StoragePolicy< Function > StoragePolicyType;

    public:
      /** \copydoc Dune::Fem::Function::DomainType */
      typedef typename BaseType::DomainType DomainType;

      /** \name Constructon
       *  \{
       */

      InstationaryFunction ( const Function &function, double time )
        : BaseType( time ),
          StoragePolicyType( function )
      {}

      InstationaryFunction ( Function &&function, double time )
        : BaseType( time ),
          StoragePolicyType( std::move( function ) )
      {}

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::Function::evaluate */
      void evaluate ( const DomainType &x, typename BaseType::RangeType &value ) const
      {
        this->get().evaluate( x, this->time(), value );
      }

      /** \copydoc Dune::Fem::Function::jacobian */
      void jacobian ( const DomainType &x, typename BaseType::JacobianRangeType &jacobian ) const
      {
        this->get().jacobian( x, this->time(), jacobian );
      }

      /** \copydoc Dune::Fem::Function::hessian */
      void hessian ( const DomainType &x, typename BaseType::HessianRangeType &hessian ) const
      {
        this->get().hessian( x, this->time(), hessian );
      }

      /** \} */
    };



    // instationaryFunction
    // --------------------

    template< class Function >
    InstationaryFunction< Function, __InstationaryFunction::HoldCopy >
    instationaryFunction ( Function function, double time )
    {
      typedef InstationaryFunction< Function, __InstationaryFunction::HoldCopy > InstationaryFunctionType;
      return InstationaryFunctionType( std::move( function ), time );
    }

    template< class Function >
    InstationaryFunction< typename std::remove_const< Function >::type, __InstationaryFunction::HoldReference >
    instationaryFunction ( std::reference_wrapper< Function > function, double time )
    {
      typedef InstationaryFunction< typename std::remove_const< Function >::type, __InstationaryFunction::HoldReference > InstationaryFunctionType;
      return InstationaryFunctionType( function.get(), time );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_COMMON_INSTATIONARY_HH
