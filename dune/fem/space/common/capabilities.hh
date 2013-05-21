#ifndef DUNE_FEM_SPACE_COMMON_CAPABILITIES_HH
#define DUNE_FEM_SPACE_COMMON_CAPABILITIES_HH


namespace Dune
{
  
  namespace Fem
  {

    namespace Capabilities
    {

      /** \class hasFixedPolynomialOrder
       *
       *  \brief specialize with \b true if polynomial order does
       *         not depend on the grid (part) entity
       */
      template< class DiscreteFunctionSpace >
      struct hasFixedPolynomialOrder 
      {
        static const bool v = false;
      };



      /** \class hasStaticPolynomialOrder
       *
       *  \brief specialize with \b true if polynomial order fixed
       *         and compile time static
       */
      template< class DiscreteFunctionSpace >
      struct hasStaticPolynomialOrder 
      {
        static const bool v = false;
        static const int order = -1;
      };



      /** \class isContinuous 
       *
       *  \brief specialize with \b true if space is always continuous
       */
      template< class DiscreteFunctionSpace >
      struct isContinuous 
      {
        static const bool v = false;
      };



      /** \class isLocalized 
       *
       *  \brief specialize with \b true if the space is localized, *
       *  i.e., the basis function set is based on a shape function set.
       *
       *  We require, that a localized space has a method
\code
  ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity );
\endcode
       */
      template< class DiscreteFunctionSpace >
      struct isLocalized 
      {
        static const bool v = false;
      };



      /** \class isParallel
       *
       *  \brief specialize with \b true if space can be used in parallel
       */
      template< class DiscreteFunctionSpace >
      struct isParallel 
      {
        static const bool v = false;
      };



      /** \class isAdaptive
       *
       *  \brief specialize with \b true if space can be used with
       *         AdaptiveDiscreteFunction
       */
      template< class DiscreteFunctionSpace >
      struct isAdaptive
      {
        static const bool v = false;
      };



      /** \class threadSafe
       *
       *  \brief specialize with \b true if the space implementation 
       *         is thread safe
       */
      template< class DiscreteFunctionSpace >
      struct threadSafe 
      {
        static const bool v = false;
      };



      /** \class viewThreadSafe 
       *
       * \brief specialize with \b true if the space implementation is 
       *        thread safe, while it is not modified
       *
       */
      template< class DiscreteFunctionSpace >
      struct viewThreadSafe 
      {
        static const bool v = false;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune 

#endif // #ifndef DUNE_FEM_SPACE_COMMON_CAPABILITIES_HH
