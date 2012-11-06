#ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_CAPABILITIES_HH
#define DUNE_FEM_SPACE_SHAPEFUNCTIONSET_CAPABILITIES_HH

/**
  @file
  @author Christoph Gersbacher
  @brief Shape function set capabilities
*/

namespace Dune
{

  namespace Fem
  {

    namespace ShapeFunctionSetCapabilities
    {

      /** \brief specialize with 'true' for if the shape function set
       *         has static size (default = false, size = undefined)
       */
      template< class ShapeFunctionSet >
      struct hasStaticSize
      {
        static const bool v = false;
        static const unsigned int size = ~0u;
      };

    } // namespace ShapeFunctionSetCapabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_CAPABILITIES_HH
