#ifndef DUNE_FEM_COMMON_LOCALCONTRIBUTION_HH
#define DUNE_FEM_COMMON_LOCALCONTRIBUTION_HH

namespace Dune
{

  namespace Fem
  {

    // Forward Declarations
    // --------------------

    template< class DiscreteFunction, template< class > class AssemblyOperation, class = void >
    class LocalContribution;



    namespace Assembly
    {

      namespace Global
      {

        // Forward Declarations
        // --------------------

        template< class DiscreteFunction, class = void >
        struct AddBase;

        template< class DiscreteFunction, class = void >
        struct SetBase;



        // Add
        // ---

        template< class DiscreteFunction >
        struct Add
          : public AddBase< DiscreteFunction >
        {};



        // Set
        // ---

        template< class DiscreteFunction >
        struct Set
          : public SetBase< DiscreteFunction >
        {};

      } // namespace Global



      // Forward Declarations
      // --------------------

      template< class DiscreteFunction, class = void >
      struct AddBase;

      template< class DiscreteFunction, class = void >
      struct AddScaledBase;

      template< class DiscreteFunction, class = void >
      struct SetBase;

      template< class DiscreteFunction, class = void >
      struct SetSelectedBase;


      // Add
      // ---

      template< class DiscreteFunction >
      struct Add
        : public AddBase< DiscreteFunction >
      {
        using AddBase< DiscreteFunction >::AddBase;
      };



      // AddScaled
      // ---------

      template< class DiscreteFunction >
      struct AddScaled
        : public AddScaledBase< DiscreteFunction >
      {
        using AddScaledBase< DiscreteFunction >::AddScaledBase;
      };



      // Set
      // ---

      template< class DiscreteFunction >
      struct Set
        : public SetBase< DiscreteFunction >
      {
        using SetBase< DiscreteFunction >::SetBase;
      };


      // SetSelected
      // ---

      template< class DiscreteFunction >
      struct SetSelected
        : public SetSelectedBase< DiscreteFunction >
      {
        using SetSelectedBase< DiscreteFunction >::SetSelectedBase;
      };

    } // namespace Assembly



    // AddLocalContribution
    // --------------------

    template< class DiscreteFunction >
    using AddLocalContribution = LocalContribution< DiscreteFunction, Assembly::Add >;



    // AddScaledLocalContribution
    // --------------------------

    template< class DiscreteFunction >
    using AddScaledLocalContribution = LocalContribution< DiscreteFunction, Assembly::AddScaled >;



    // SetLocalContribution
    // --------------------

    template< class DiscreteFunction >
    using SetLocalContribution = LocalContribution< DiscreteFunction, Assembly::Set >;


    // SetSelectedLocalContribution
    // ----------------------------

    template< class DiscreteFunction >
    using SetSelectedLocalContribution = LocalContribution< DiscreteFunction, Assembly::SetSelected >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_LOCALCONTRIBUTION_HH
