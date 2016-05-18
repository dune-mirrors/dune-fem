#ifndef DUNE_FEMPY_GRID_RESTRICTPROLONG_HH
#define DUNE_FEMPY_GRID_RESTRICTPROLONG_HH

#include <cstddef>

#include <memory>
#include <utility>

#include <dune/fem/space/common/restrictprolonginterface.hh>

namespace Dune
{

  namespace FemPy
  {

    // RestrictProlongInterface
    // ------------------------

    template< class Grid >
    struct RestrictProlongInterface
    {
      typedef typename Grid::ctype ctype;

      typedef typename Grid::template Codim< 0 >::Entity Element;
      typedef typename Grid::template Codim< 0 >::LocalGeometry LocalGeometry;

      virtual ~RestrictProlongInterface () = default;

      virtual void setFatherChildWeight ( const ctype &weight ) = 0;

      virtual void restrictLocal ( const Element &father, const Element &child, bool initialize ) = 0;
      virtual void restrictLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) = 0;

      virtual void prolongLocal ( const Element &father, const Element &child, bool initialize ) = 0;
      virtual void prolongLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) = 0;

      virtual void addToList ( Fem::CommunicationManagerList &commList ) = 0;
      virtual void removeFromList ( Fem::CommunicationManagerList &commList ) = 0;
    };



    // RestrictProlongImpl
    // -------------------

    template< class Grid, class Impl >
    class RestrictProlongImpl final
      : public RestrictProlongInterface< Grid >
    {
      typedef RestrictProlongInterface< Grid > Base;

    public:
      typedef typename Base::ctype ctype;

      typedef typename Base::Element Element;
      typedef typename Base::LocalGeometry LocalGeometry;

      template< class ... Args >
      RestrictProlongImpl ( Args && ... args )
        : impl_( std::forward< Args >( args ) ... )
      {}

      virtual void setFatherChildWeight ( const ctype &weight ) override
      {
        impl_.setFatherChildWeight( weight );
      }

      virtual void restrictLocal ( const Element &father, const Element &child, bool initialize ) override
      {
        impl_.restrictLocal( father, child, initialize );
      }

      virtual void restrictLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) override
      {
        impl_.restrictLocal( father, child, geometryInFather, initialize );
      }

      virtual void prolongLocal ( const Element &father, const Element &child, bool initialize ) override
      {
        impl_.prolongLocal( father, child, initialize );
      }

      virtual void prolongLocal ( const Element &father, const Element &child, const LocalGeometry &geometryInFather, bool initialize ) override
      {
        impl_.prolongLocal( father, child, geometryInFather, initialize );
      }

      virtual void addToList ( Fem::CommunicationManagerList &commList ) override
      {
        //  impl_.addToList( commList );
      }

      virtual void removeFromList ( Fem::CommunicationManagerList &commList ) override
      {
        // impl_.removeFromList( commList );
      }

    private:
      Impl impl_;
    };



    // makeDefaultRestrictProlong
    // --------------------------

    template< class DiscreteFunction >
    inline static std::shared_ptr< RestrictProlongInterface< typename DiscreteFunction::GridPartType::GridType > >
    makeDefaultRestrictProlong ( DiscreteFunction &discreteFunction )
    {
      typedef RestrictProlongImpl< typename DiscreteFunction::GridPartType::GridType, Fem::RestrictProlongDefault< DiscreteFunction > > RP;
      return std::make_shared< RP >( discreteFunction );
    }



    // AdaptiveDofVector
    // -----------------

    template< class Grid, class D = double >
    struct AdaptiveDofVector
    {
      typedef D DofType;

      typedef typename Grid::template Codim< 0 >::Entity Element;

      typedef RestrictProlongInterface< Grid > RestrictProlong;

      explicit AdaptiveDofVector ( std::shared_ptr< RestrictProlong > restrictProlong )
        : restrictProlong_( std::move( restrictProlong ) )
      {}

      virtual ~AdaptiveDofVector () = default;

      virtual void enableDofCompression () = 0;

      virtual std::size_t numLocalDofs ( const Element &element ) const = 0;

      virtual DofType *getLocalDofs ( const Element &element, DofType *localDofs ) const = 0;
      virtual const DofType *setLocalDofs ( const Element &element, const DofType *localDofs ) = 0;

      RestrictProlong &restrictProlong () { assert( restrictProlong_ ); return *restrictProlong_; }

    protected:
      void resetRestrictProlong ( std::shared_ptr< RestrictProlong > restrictProlong )
      {
        restrictProlong_ = std::move( restrictProlong );
      }

    private:
      std::shared_ptr< RestrictProlong > restrictProlong_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_GRID_RESTRICTPROLONG_HH
