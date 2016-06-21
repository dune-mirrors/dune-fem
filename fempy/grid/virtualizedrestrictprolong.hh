#ifndef DUNE_FEMPY_GRID_VIRTUALIZEDRESTRICTPROLONG_HH
#define DUNE_FEMPY_GRID_VIRTUALIZEDRESTRICTPROLONG_HH

#include <functional>
#include <memory>
#include <type_traits>
#include <utility>

#include <dune/fem/space/common/restrictprolonginterface.hh>

namespace Dune
{

  namespace FemPy
  {

    // VirtualizedRestrictProlong
    // --------------------------

    template< class Grid >
    struct VirtualizedRestrictProlong
    {
      typedef typename Grid::ctype ctype;

      typedef typename Grid::template Codim< 0 >::Entity Element;
      typedef typename Grid::template Codim< 0 >::LocalGeometry LocalGeometry;

    private:
      struct Interface
      {
        virtual ~Interface () = default;

        virtual void setFatherChildWeight ( const ctype &weight ) const = 0;

        virtual void restrictLocal ( const Element &father, const Element &child, bool initialize ) const = 0;
        virtual void restrictLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) const = 0;

        virtual void prolongLocal ( const Element &father, const Element &child, bool initialize ) const = 0;
        virtual void prolongLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) const = 0;

        virtual void addToList ( Fem::CommunicationManagerList &commList ) = 0;
        virtual void removeFromList ( Fem::CommunicationManagerList &commList ) = 0;

        virtual void addToLoadBalancer ( Fem::LoadBalancer< Grid > &lb ) = 0;
      };

      template< class Impl >
      struct Implementation final
      {
        template< class... Args >
        Implementation ( Args &&... args )
          : impl_( std::forward< Args >( args )... )
        {}

        virtual void setFatherChildWeight ( const ctype &weight ) const override { impl().setFatherChildWeight( weight ); }

        virtual void restrictLocal ( const Element &father, const Element &child, bool initialize ) const override { impl().restrictLocal( father, child, initialize ); }
        virtual void restrictLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) const override { impl().restrictLocal( father, child, geometryInFather, initialize ); }

        virtual void prolongLocal ( const Element &father, const Element &child, bool initialize ) const override { impl().prolongLocal( father, child, initialize ); }
        virtual void prolongLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) const override { impl().prolongLocal( father, child, geometryInFather, initialize ); }

        virtual void addToList ( Fem::CommunicationManagerList &commList ) override { impl().addToList( commList ); }
        virtual void removeFromList ( Fem::CommunicationManagerList &commList ) override { impl().removeFromList( commList ); }

        virtual void addToLoadBalancer ( Fem::LoadBalancer< Grid > &lb ) override { impl().addToLoadBalancer( lb ); }

      private:
        const auto &impl () const { using std::cref; return cref( impl_ ).get(); }
        auto &impl () { using std::ref; return ref( impl_ ).get(); }

        Impl impl_;
      };

    public:
      template< class DF, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DF >::value, int > = 0 >
      explicit VirtualizedRestrictProlong ( DF &df )
        : impl_( new Implementation< Fem::RestrictProlongDefault< DF > >( df ) )
      {}

      void setFatherChildWeight ( const ctype &weight ) const { impl_->setFatherChildWeight( weight ); }

      void restrictLocal ( const Element &father, const Element &child, bool initialize ) const { impl_->restrictLocal( father, child, initialize ); }
      void restrictLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) const { impl_->restrictLocal( father, child, geometryInFather, initialize ); }

      void prolongLocal ( const Element &father, const Element &child, bool initialize ) const { impl_->prolongLocal( father, child, initialize ); }
      void prolongLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) const { impl_->prolongLocal( father, child, geometryInFather, initialize ); }

      void addToList ( Fem::CommunicationManagerList &commList ) { impl_->addToList( commList ); }
      void removeFromList ( Fem::CommunicationManagerList &commList ) { impl_->removeFromList( commList ); }

      void addToLoadBalancer ( Fem::LoadBalancer< Grid > &lb ) { impl_->addToLoadBalancer( lb ); }

    private:
      std::unique_ptr< Interface > impl_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRID_VIRTUALIZEDRESTRICTPROLONG_HH
