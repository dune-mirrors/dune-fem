#ifndef DUNE_FEM_PASSS_DISCRETEFUNCTIONTUPLE_HH
#define DUNE_FEM_PASSS_DISCRETEFUNCTIONTUPLE_HH

#include <dune/common/typetraits.hh>

#include "localfunctionselector.hh"
#include "tupletypetraits.hh"

namespace
{

  // LocalFunctionEvaluator
  // ----------------------

  /*
   * \brief Choose local function for given discrete function.
   *        Use LocalFunctionSelector to define the type of 
   *        the local function.
   */
  template< class DiscreteFunction >
  struct LocalFunctionEvaluator
  {
    typedef typename Dune::Fem::LocalFunctionSelector< 
        typename Dune::TypeTraits< DiscreteFunction >::ReferredType
      >::Type Type;

    static Type apply ( const DiscreteFunction &discreteFunction )
    {
      return Type( discreteFunction );
    }
  };

} // namespace



namespace Dune
{
  
  namespace Fem
  {
    
    // DiscreteFunctionTuple
    // ---------------------

    template< class DiscreteFunctions >
    class DiscreteFunctionTuple
    {
      typedef DiscreteFunctionTuple< DiscreteFunctions > ThisType;
  
    public:
      //! \brief forward template argument
      typedef DiscreteFunctions DiscreteFunctionTupleType;

      //! \brief type of local function tuple
      typedef typename Dune::ForEachType< LocalFunctionEvaluator, DiscreteFunctionTupleType >::Type LocalFunctionTupleType;

      DiscreteFunctionTuple ( typename Dune::ReferenceTuple< DiscreteFunctions > discreteFunctions )
      : discreteFunctions_( discreteFunctions )
      {}

      //! \brief return tuple of discrete function references
      typename Dune::ReferenceTuple< DiscreteFunctions > discreteFunctions ()
      {
        return discreteFunctions_;
      }

      //! \brief return tuple of discrete function references
      const typename Dune::ReferenceTuple< DiscreteFunctions > discreteFunctions () const
      {
        return discreteFunctions_;
      }

    private:
      typename Dune::ReferenceTuple< DiscreteFunctions > discreteFunctions_;
    };

   } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASSS_DISCRETEFUNCTIONTUPLE_HH
