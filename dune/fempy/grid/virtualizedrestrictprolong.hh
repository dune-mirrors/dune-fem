#ifndef DUNE_FEMPY_GRID_VIRTUALIZEDRESTRICTPROLONG_HH
#define DUNE_FEMPY_GRID_VIRTUALIZEDRESTRICTPROLONG_HH

#include <functional>
#include <memory>
#include <type_traits>
#include <utility>

#include <dune/fem/space/common/loadbalancer.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>

#include <dune/fempy/grid/discretefunctionmanager.hh>

namespace Dune
{

  namespace FemPy
  {

    template< class T >
    using DefaultUniquePtr = std::unique_ptr< T >;



    // BasicVirtualizedRestrictProlong
    // -------------------------------

    template< class Grid, template< class > class Holder >
    struct BasicVirtualizedRestrictProlong
    {
      typedef typename Grid::ctype ctype;

      typedef typename Grid::template Codim< 0 >::Entity Element;

    protected:
      struct Interface
      {
        virtual ~Interface () = default;
        virtual Interface *clone () const = 0;

        virtual void initialize () = 0;
        virtual void finalize   () = 0;

        virtual void setFatherChildWeight ( const ctype &weight ) const = 0;

        virtual void restrictLocal ( const Element &father, const Element &child, bool initialize ) const = 0;
        virtual void restrictFinalize ( const Element &father ) const = 0;
        virtual void prolongLocal ( const Element &father, const Element &child, bool initialize ) const = 0;

        virtual void addToList ( Fem::CommunicationManagerList &commList ) = 0;
        virtual void removeFromList ( Fem::CommunicationManagerList &commList ) = 0;

        virtual void addToLoadBalancer ( Fem::LoadBalancer< Grid > &lb ) = 0;
        virtual void addToLoadBalancer ( DiscreteFunctionManager< Grid > &lb ) = 0;
      };

      template< class Impl >
      struct Implementation final
        : public Interface
      {
        template< class... Args >
        Implementation ( Args &&... args )
          : impl_( std::forward< Args >( args )... )
        {}

        virtual Interface *clone () const override { return new Implementation( *this ); }

        virtual void initialize () override { impl().initialize(); }
        virtual void finalize   () override { impl().finalize(); }

        virtual void setFatherChildWeight ( const ctype &weight ) const override { impl().setFatherChildWeight( weight ); }

        virtual void restrictLocal ( const Element &father, const Element &child, bool initialize ) const override { impl().restrictLocal( father, child, initialize ); }
        virtual void restrictFinalize ( const Element &father ) const override { impl().restrictFinalize( father ); }
        virtual void prolongLocal ( const Element &father, const Element &child, bool initialize ) const override { impl().prolongLocal( father, child, initialize ); }

        virtual void addToList ( Fem::CommunicationManagerList &commList ) override { impl().addToList( commList ); }
        virtual void removeFromList ( Fem::CommunicationManagerList &commList ) override { impl().removeFromList( commList ); }

        virtual void addToLoadBalancer ( Fem::LoadBalancer< Grid > &lb ) override { impl().addToLoadBalancer( lb ); }
        virtual void addToLoadBalancer ( DiscreteFunctionManager< Grid > &lb ) override { impl().addToLoadBalancer( lb ); }

      private:
        const auto &impl () const { using std::cref; return cref( impl_ ).get(); }
        auto &impl () { using std::ref; return ref( impl_ ).get(); }

        Impl impl_;
      };

    protected:
      BasicVirtualizedRestrictProlong ( Interface *impl ) : impl_( impl ) {}

    public:
      template< class DF, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DF >::value, int > = 0 >
      explicit BasicVirtualizedRestrictProlong ( DF &df )
        : impl_( new Implementation< Fem::RestrictProlongDefault< DF > >( df ) )
      {}

      void initialize () { impl_->initialize(); }
      void finalize   () { impl_->finalize(); }

      void setFatherChildWeight ( const ctype &weight ) const { impl_->setFatherChildWeight( weight ); }

      void restrictLocal ( const Element &father, const Element &child, bool initialize ) const { impl_->restrictLocal( father, child, initialize ); }
      void restrictFinalize ( const Element &father ) const { impl_->restrictFinalize( father ); }
      void prolongLocal ( const Element &father, const Element &child, bool initialize ) const { impl_->prolongLocal( father, child, initialize ); }

      void addToList ( Fem::CommunicationManagerList &commList ) { impl_->addToList( commList ); }
      void removeFromList ( Fem::CommunicationManagerList &commList ) { impl_->removeFromList( commList ); }

      void addToLoadBalancer ( Fem::LoadBalancer< Grid > &lb ) { impl_->addToLoadBalancer( lb ); }
      void addToLoadBalancer ( DiscreteFunctionManager< Grid > &lb ) { impl_->addToLoadBalancer( lb ); }

    protected:
      Holder< Interface > impl_;
    };



    // VirtualizedRestrictProlong
    // --------------------------

    template< class Grid, bool shared = false >
    class VirtualizedRestrictProlong;

    template< class Grid >
    class VirtualizedRestrictProlong< Grid, false >
      : public BasicVirtualizedRestrictProlong< Grid, DefaultUniquePtr >
    {
      typedef BasicVirtualizedRestrictProlong< Grid, DefaultUniquePtr > Base;

    public:
      template< class DF, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DF >::value, int > = 0 >
      explicit VirtualizedRestrictProlong ( DF &df )
        : Base( df )
      {}

      VirtualizedRestrictProlong ( const VirtualizedRestrictProlong &other ) : Base( other.impl_->clone() ) {}
      VirtualizedRestrictProlong ( VirtualizedRestrictProlong && ) = default;

      VirtualizedRestrictProlong &operator= ( const VirtualizedRestrictProlong &other ) { Base::impl_.reset( other.impl_->clone() ); return *this; }
      VirtualizedRestrictProlong &operator= ( VirtualizedRestrictProlong && ) = default;
    };

    template< class Grid >
    class VirtualizedRestrictProlong< Grid, true >
      : public BasicVirtualizedRestrictProlong< Grid, std::shared_ptr >
    {
      typedef BasicVirtualizedRestrictProlong< Grid, std::shared_ptr > Base;

    public:
      template< class DF, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DF >::value, int > = 0 >
      explicit VirtualizedRestrictProlong ( DF &df )
        : Base( df )
      {}
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRID_VIRTUALIZEDRESTRICTPROLONG_HH
