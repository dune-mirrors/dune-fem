#ifndef DUNE_FEM_SPACE_MAPPER_CAPABILITIES_HH
#define DUNE_FEM_SPACE_MAPPER_CAPABILITIES_HH


namespace Dune
{

  namespace Fem
  {

    namespace Capabilities
    {

      /** \class isAdaptiveDofeMapper
       *
       *  \brief specialize with \b true if the mapper supports adaptivity
       *
       **/

      template< class Mapper >
      struct isAdaptiveDofMapper
      {
        static const bool v = false;
      };


      // const specialization
      // --------------------

      template< class Mapper >
      struct isAdaptiveDofMapper< const Mapper >
      {
        static const bool v = isAdaptiveDofMapper< Mapper >::v;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_MAPPER_CAPABILITIES_HH
