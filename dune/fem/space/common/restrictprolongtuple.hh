#ifndef DUNE_FEM_SPACE_COMMON_RESTRICTPROLONGTUPLE_HH
#define DUNE_FEM_SPACE_COMMON_RESTRICTPROLONGTUPLE_HH

#include <tuple>
#include <utility>

#include <dune/common/deprecated.hh>
#include <dune/common/forloop.hh>
#include <dune/common/tupleutility.hh>

#include <dune/fem/space/common/restrictprolonginterface.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------

    template< class... RestrictProlongInterfaces >
    class RestrictProlongTuple;
    template< class... DiscreteFunctions >
    class RestrictProlongDefaultTuple;



    // RestrictProlongTuple
    // --------------------

    /** \addtogroup RestrictProlongInterface
     *  \{
     **/

    /** \class RestrictProlongTuple
     *
     *  \brief combine a variadic number of Dune::Fem::RestrictProlongInterface
     *         instances into a single object again derived from
     *         Dune::Fem::RestrictProlongInterface
     *
     *  \tparam  Head  a Dune::Fem::RestrictProlongInterface type
     *  \tparam  Tail  additional Dune::Fem::RestrictProlongInterface types
     */
    template< class Head, class... Tail >
    class RestrictProlongTuple< Head, Tail... >
    : public Dune::Fem::RestrictProlongInterface< RestrictProlongTraits< RestrictProlongTuple< Head, Tail... >, typename Head::DomainFieldType > >
    {
      typedef Dune::Fem::RestrictProlongInterface< RestrictProlongTraits< RestrictProlongTuple< Head, Tail... >, typename Head::DomainFieldType > > BaseType;

      template< int i > struct AddToList;
      template< int i > struct AddToLoadBalancer;
      template< int i > struct ProlongLocal;
      template< int i > struct RestrictLocal;
      template< int i > struct SetFatherChildWeight;

    public:
      /** \copydoc Dune::Fem::RestrictProlongInterface::DomainFieldType */
      typedef typename BaseType::DomainFieldType DomainFieldType;

      /** \name Construction
       *  \{
       */

      explicit RestrictProlongTuple ( Head &&head, Tail &&... tail )
        : tuple_( std::forward< Head >( head ), std::forward< Tail >( tail )... )
      {}

      explicit RestrictProlongTuple ( std::tuple< Head, Tail... > &&tuple )
        : tuple_( tuple )
      {}

      /** \} */

      /** \name Interface methods
       *  \{
       */

      /** \copydoc Dune::Fem::RestrictProlongInterface::setFatherChildWeight */
      void setFatherChildWeight ( const DomainFieldType &weight ) const
      {
        Dune::ForLoop< SetFatherChildWeight, 0, sizeof...( Tail ) >::apply( weight, tuple_ );
      }

      /** \copydoc Dune::Fem::RestrictProlongInterface::restrictLocal */
      template< class Entity >
      void restrictLocal ( const Entity &father, const Entity &child, bool initialize ) const
      {
        Dune::ForLoop< RestrictLocal, 0, sizeof...( Tail ) >::apply( father, child, initialize, tuple_ );
      }

      /** \copydoc Dune::Fem::RestrictProlongInterface::restrictLocal */
      template< class Entity, class LocalGeometry >
      void restrictLocal ( const Entity &father, const Entity &child,
                           const LocalGeometry &geometryInFather, bool initialize ) const
      {
        Dune::ForLoop< RestrictLocal, 0, sizeof...( Tail ) >::apply( father, child, geometryInFather, initialize, tuple_ );
      }

      /** \copydoc Dune::Fem::RestrictProlongInterface::prolongLocal */
      template< class Entity >
      void prolongLocal ( const Entity &father, const Entity &child, bool initialize ) const
      {
        Dune::ForLoop< ProlongLocal, 0, sizeof...( Tail ) >::apply( father, child, initialize, tuple_ );
      }

      /** \copydoc Dune::Fem::RestrictProlongInterface::prolongLocal */
      template< class Entity, class LocalGeometry >
      void prolongLocal ( const Entity &father, const Entity &child,
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {
        Dune::ForLoop< ProlongLocal, 0, sizeof...( Tail ) >::apply( father, child, geometryInFather, initialize, tuple_ );
      }

      /** \copydoc Dune::Fem::RestrictProlongInterface::addToList */
      template< class Communicator >
      void addToList ( Communicator &comm )
      {
        Dune::ForLoop< AddToList, 0, sizeof...( Tail ) >::apply( comm, tuple_ );
      }

      /** \copydoc Dune::Fem::RestrictProlongInterface::addToLoadBalancer */
      template< class LoadBalancer >
      void addToLoadBalancer ( LoadBalancer &loadBalancer )
      {
        Dune::ForLoop< AddToLoadBalancer, 0, sizeof...( Tail ) >::apply( loadBalancer, tuple_ );
      }

      /** \} */

    private:
      std::tuple< Head, Tail... > tuple_;
    };



    // RestrictProlongTuple< Head, Tail... >::AddToList
    // ------------------------------------------------

    template< class Head, class... Tail >
    template< int i >
    struct RestrictProlongTuple< Head, Tail... >::AddToList
    {
      template< class Communicator >
      static void apply ( Communicator &comm, std::tuple< Head, Tail... > &tuple )
      {
        std::get< i >( tuple ).addToList( comm );
      }
    };



    // RestrictProlongTuple< Head, Tail... >::AddToLoadBalancer
    // --------------------------------------------------------

    template< class Head, class... Tail >
    template< int i >
    struct RestrictProlongTuple< Head, Tail... >::AddToLoadBalancer
    {
      template< class LoadBalancer >
      static void apply ( LoadBalancer &loadBalancer, std::tuple< Head, Tail... > &tuple )
      {
        std::get< i >( tuple ).addToLoadBalancer( loadBalancer );
      }
    };



    // RestrictProlongTuple< Head, Tail... >::ProlongLocal
    // ---------------------------------------------------

    template< class Head, class... Tail >
    template< int i >
    struct RestrictProlongTuple< Head, Tail... >::ProlongLocal
    {
      template< class Entity >
      static void apply ( const Entity &father, const Entity &child, bool initialize,
                          const std::tuple< Head, Tail... > &tuple )
      {
        std::get< i >( tuple ).prolongLocal( father, child, initialize );
      }

      template< class Entity, class LocalGeometry >
      static void apply ( const Entity &father, const Entity &child, const LocalGeometry &geometryInFather, bool initialize,
                          const std::tuple< Head, Tail... > &tuple )
      {
        std::get< i >( tuple ).prolongLocal( father, child, geometryInFather, initialize );
      }
    };



    // RestrictProlongTuple< Head, Tail... >::RestrictLocal
    // ----------------------------------------------------

    template< class Head, class... Tail >
    template< int i >
    struct RestrictProlongTuple< Head, Tail... >::RestrictLocal
    {
      template< class Entity >
      static void apply ( const Entity &father, const Entity &child, bool initialize,
                          const std::tuple< Head, Tail... > &tuple )
      {
        std::get< i >( tuple ).restrictLocal( father, child, initialize );
      }

      template< class Entity, class LocalGeometry >
      static void apply ( const Entity &father, const Entity &child, const LocalGeometry &geometryInFather, bool initialize,
                          const std::tuple< Head, Tail... > &tuple )
      {
        std::get< i >( tuple ).restrictLocal( father, child, geometryInFather, initialize );
      }
    };



    // RestrictProlongTuple< Head, Tail... >::SetFatherChildWeight
    // -----------------------------------------------------------

    template< class Head, class... Tail >
    template< int i >
    struct RestrictProlongTuple< Head, Tail... >::SetFatherChildWeight
    {
      static void apply ( const DomainFieldType &weight, const std::tuple< Head, Tail... > &tuple )
      {
        std::get< i >( tuple ).setFatherChildWeight( weight );
      }
    };



    // RestrictProlongDefaultTuple
    // ---------------------------

    /** \class RestrictProlongDefaultTuple
     *
     *  \brief conveniently set up a tuple of Dune::Fem::RestrictProlongDefault
     *         restriction/prolongation objects created from a variadic list of
     *         discrete functions
     *
     *  \tparam  DiscreteFunctions  a variadic list of discrete function types
     */
    template< class... DiscreteFunctions >
    class RestrictProlongDefaultTuple
    : public RestrictProlongTuple< RestrictProlongDefault< DiscreteFunctions >... >
    {
      typedef RestrictProlongTuple< RestrictProlongDefault< DiscreteFunctions >... > BaseType;

      template< class DiscreteFunction >
      struct Operation
      {
        typedef typename std::decay< DiscreteFunction >::type DiscreteFunctionType;
        typedef RestrictProlongDefault< DiscreteFunctionType > Type;

        static Type apply ( DiscreteFunctionType &discreteFunction )
        {
          return Type( discreteFunction );
        }
      };

    public:
      explicit RestrictProlongDefaultTuple ( DiscreteFunctions &... discreteFunctions )
        : BaseType( RestrictProlongDefault< DiscreteFunctions >( discreteFunctions )... )
      {}

      explicit RestrictProlongDefaultTuple ( std::tuple< DiscreteFunctions &... > tuple )
        : BaseType( Dune::transformTuple< Operation >( tuple ) )
      {}
    };



    // makeRestrictProlongDefault
    // --------------------------

    /** \fn makeRestrictProlongDefault
     *
     *  \brief conveniently set up a tuple of Dune::Fem::RestrictProlongDefault
     *         restriction/prolongation objects created from a variadic list of
     *         discrete functions
     *
     *  \param[in]  discreteFunctions  a variadic list of discrete function types
     */
    template< class... DiscreteFunctions >
    static inline RestrictProlongDefaultTuple< DiscreteFunctions... >
    DUNE_DEPRECATED makeRestrictProlongDefault ( DiscreteFunctions &... discreteFunctions )
    {
      return RestrictProlongDefaultTuple< DiscreteFunctions... >( discreteFunctions... );
    }

    template< class... DiscreteFunctions >
    static inline RestrictProlongDefaultTuple< DiscreteFunctions... >
    DUNE_DEPRECATED makeRestrictProlongDefault ( std::tuple< DiscreteFunctions &... > tuple )
    {
      return RestrictProlongDefaultTuple< DiscreteFunctions ... >( tuple );
    }

    /** \} */

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMMON_RESTRICTPROLONGTUPLE_HH
