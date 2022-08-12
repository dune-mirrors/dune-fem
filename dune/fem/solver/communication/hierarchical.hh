#ifndef DUNE_FEM_SOLVER_COMMUNICATION_HIERARCHICAL_HH
#define DUNE_FEM_SOLVER_COMMUNICATION_HIERARCHICAL_HH

#include <memory>
#include <type_traits>
#include <utility>

#include <dune/common/ftraits.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/solvercategory.hh>

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/function/hierarchical/dofvector.hh>
#include <dune/fem/space/common/commoperations.hh>

namespace Dune
{

  namespace Fem
  {

    namespace ISTL
    {

      // HierarchicalCommunicationVector
      // -------------------------------

      template< class DiscreteFunctionSpace, class DofContainer >
      class HierarchicalCommunicationVector
      {
        typedef HierarchicalCommunicationVector< DiscreteFunctionSpace, DofContainer > ThisType;

      public:
        struct DofVector
        {
          typedef typename Dune::Fem::Impl::BlockIndicesFor< DofContainer >::Type BlockIndices;
          static constexpr std::size_t blockSize = Dune::Hybrid::size( BlockIndices() );

          typedef HierarchicalDofBlock< const DofContainer > ConstDofBlockType;
          typedef HierarchicalDofBlock< DofContainer > DofBlockType;

          explicit DofVector ( DofContainer &data ) : data_( data ) {}

          ConstDofBlockType operator[] ( std::size_t i ) const { return ConstDofBlockType( data_, i ); }
          DofBlockType operator[] ( std::size_t i ) { return DofBlockType( data_, i ); }

        private:
          DofContainer &data_;
        };

        typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

        typedef typename DiscreteFunctionSpaceType::LocalBlockIndices BlockIndices;

        static constexpr std::size_t blockSize = Hybrid::size( BlockIndices() );

        typedef typename DofContainer::field_type DofType;

        template< class Operation >
        struct CommDataHandle
        {
          typedef typename DiscreteFunctionSpaceType::template CommDataHandle< ThisType, Operation >::Type Type;
        };

        HierarchicalCommunicationVector ( const DiscreteFunctionSpace &dfSpace, DofContainer &dofContainer )
          : dfSpace_( dfSpace ), dofVector_( dofContainer )
        {}

        template< class Operation >
        typename CommDataHandle< Operation >::Type dataHandle ( const Operation &operation )
        {
          return space().createDataHandle( *this, operation );
        }

        const DofVector &dofVector () const { return dofVector_; }
        DofVector &dofVector () { return dofVector_; }

        const DiscreteFunctionSpaceType &space () const { return dfSpace_; }

      private:
        const DiscreteFunctionSpace &dfSpace_;
        DofVector dofVector_;
      };



      // HierarchicalCommunication
      // -------------------------

      template< class DiscreteFunctionSpace >
      class HierarchicalCommunication
      {
        typedef HierarchicalCommunication< DiscreteFunctionSpace > ThisType;

        typedef typename DiscreteFunctionSpace::AuxiliaryDofsType AuxiliaryDofsType;

      public:
        typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

        explicit HierarchicalCommunication ( const DiscreteFunctionSpaceType &dfSpace,  Dune::SolverCategory::Category solverCategory = Dune::SolverCategory::sequential )
          : dfSpace_( dfSpace ), solverCategory_( solverCategory )
        {}

        const typename DiscreteFunctionSpace::GridPartType::CommunicationType &communicator () const { return dfSpace_.gridPart().comm(); }

        template< class T >
        void copyOwnerToAll ( const T &x, T &y ) const
        {
          y = x;
          project( y );
          HierarchicalCommunicationVector< DiscreteFunctionSpaceType, T > z( dfSpace_, y );
          dfSpace_.communicator().exchange( z, DFCommunicationOperation::Add() );
        }

        template< class T >
        void project ( T &x ) const
        {
          project( dfSpace_.auxiliaryDofs(), x );
        }

        template< class T, class F >
        void dot ( const T &x, const T &y, F &scp ) const
        {
          dot( dfSpace_.auxiliaryDofs(), x, y, scp );
          scp = communicator().sum( scp );
        }

        template< class T >
        typename Dune::FieldTraits< typename T::field_type >::real_type norm ( const T &x ) const
        {
          using std::sqrt;
          typename Dune::FieldTraits< typename T::field_type >::real_type norm2( 0 );
          dot( x, x, norm2 );
          return sqrt( norm2 );
        }

        Dune::SolverCategory::Category getSolverCategory () const { return solverCategory_; }

      private:
        template< class... V >
        static void project ( const AuxiliaryDofsType &auxiliaryDofs, MultiTypeBlockVector< V... > &x )
        {
          Dune::Hybrid::forEach( std::index_sequence_for< V... >(), [ &auxiliaryDofs, &x ] ( auto &&i ) { ThisType::project( auxiliaryDofs, x[ i ] ); } );
        }

        template< class B, class A >
        static void project ( const AuxiliaryDofsType &auxiliaryDofs, BlockVector< B, A > &x )
        {
          typedef typename B::field_type field_type;
          for( int i : auxiliaryDofs )
            x[ i ] = field_type( 0 );
        }

        template< class... V, class F >
        static void dot ( const AuxiliaryDofsType &auxiliaryDofs, const MultiTypeBlockVector< V... > &x, const MultiTypeBlockVector< V... > &y, F &scp )
        {
          Dune::Hybrid::forEach( std::index_sequence_for< V... >(), [ &auxiliaryDofs, &x, &y, &scp ] ( auto &&i ) { ThisType::dot( auxiliaryDofs, x[ i ], y[ i ], scp ); } );
        }

        template< class B, class A, class F >
        static void dot ( const AuxiliaryDofsType &auxiliaryDofs, const BlockVector< B, A > &x, const BlockVector< B, A > &y, F &scp )
        {
          const int numAuxiliarys = auxiliaryDofs.size();
          for( int auxiliary = 0, i = 0; auxiliary < numAuxiliarys; ++auxiliary, ++i )
          {
            const int nextAuxiliary = auxiliaryDofs[ auxiliary ];
            for( ; i < nextAuxiliary; ++i )
              scp += x[ i ] * y[ i ];
          }
        }

        const DiscreteFunctionSpaceType &dfSpace_;
        Dune::SolverCategory::Category solverCategory_;
      };



      // buildCommunication
      // ------------------

      template< class DiscreteFunctionSpace >
      void buildCommunication ( const DiscreteFunctionSpace &dfSpace,
                                Dune::SolverCategory::Category solverCategory,
                                std::shared_ptr< HierarchicalCommunication< DiscreteFunctionSpace > > &communication )
      {
        communication.reset( new HierarchicalCommunication< DiscreteFunctionSpace >( dfSpace, solverCategory ) );
      }



      // SupportsAMG for HierarchicalCommunication
      // -----------------------------------------

      template< class T >
      struct SupportsAMG;

      template< class DiscreteFunctionSpace >
      struct SupportsAMG< HierarchicalCommunication< DiscreteFunctionSpace > >
        : public std::false_type
      {};

    } // namespace ISTL

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SOLVER_COMMUNICATION_HIERARCHICAL_HH
