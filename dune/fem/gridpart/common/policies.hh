#ifndef DUNE_FEM_GRIDPART_COMMON_POLICIES_HH
#define DUNE_FEM_GRIDPART_COMMON_POLICIES_HH

#warning "This header should not be needed anymore. Remove it from the include list!"


#include <type_traits>
#include <utility>

#include <dune/grid/common/gridview.hh>

#include <dune/fem/gridpart/common/gridpart2gridview.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------

    template< class Traits >
    class GridPartPolicies;



#ifndef DOXYGEN

    namespace __GridPartPolicies
    {

      // hasGridView
      // -----------

      template< class Traits >
      std::true_type __hasGridView ( const Traits &, const typename Traits::GridViewType * = nullptr );

      std::false_type __hasGridView ( ... );

      template< class Traits >
      struct hasGridView
      {
        static const bool v = decltype( __hasGridView( std::declval< Traits >() ) )::value;
      };



      // HasGridView
      // -----------

      template< class Traits, bool hasGridView = hasGridView< Traits >::v >
      class HasGridView;

      template< class Traits >
      class HasGridView< Traits, true >
      {
      public:
        typedef typename Traits::GridViewType GridViewType;
      };

      template< class Traits >
      class HasGridView< Traits, false >
      {
        typedef typename Traits::GridPartType GridPartType;

      public:
        typedef Dune::GridView< GridPart2GridViewTraits< GridPartType > > GridViewType;

      protected:
        HasGridView () = default;

      public:
        operator GridPart2GridViewImpl< GridPartType > () const
        {
          return GridPart2GridViewImpl< GridPartType >( impl() );
        }

      private:
        const GridPartType &impl () const
        {
          return static_cast< const GridPartType & >( *this );
        }
      };

    } // namespace __GridPartPolicies

#endif // #ifndef DOXYGEN



    // GridPartPolicies
    // ----------------

    template< class Traits >
    class GridPartPolicies
      : public __GridPartPolicies::HasGridView< Traits >
    {};

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_POLICIES_HH
