#ifndef DUNE_FEM_SOLVER_COMMUNICATION_FEM_HH
#define DUNE_FEM_SOLVER_COMMUNICATION_FEM_HH

#include <memory>
#include <type_traits>
#include <utility>

#include <dune/common/ftraits.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/istl/solvercategory.hh>

#include <dune/fem/space/common/commoperations.hh>

namespace Dune
{

  namespace Fem
  {

    namespace ISTL
    {

      // FemCommunicationVector
      // ----------------------

      template< class DiscreteFunctionSpace, class DofVector >
      class FemCommunicationVector
      {
        typedef FemCommunicationVector< DiscreteFunctionSpace, DofVector > ThisType;

      public:
        typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

        typedef typename DiscreteFunctionSpaceType::LocalBlockIndices BlockIndices;

        static constexpr std::size_t blockSize = Hybrid::size( BlockIndices() );

        typedef typename DofVector::field_type DofType;

        template< class Operation >
        struct CommDataHandle
        {
          typedef typename DiscreteFunctionSpaceType::template CommDataHandle< ThisType, Operation >::Type Type;
        };

        FemCommunicationVector ( const DiscreteFunctionSpace &dfSpace, DofVector &dofVector )
          : dfSpace_( dfSpace ), dofVector_( dofVector )
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
        DofVector &dofVector_;
      };



      // FemCommunication
      // ----------------

      template< class DiscreteFunctionSpace >
      class FemCommunication
      {
        typedef FemCommunication< DiscreteFunctionSpace > ThisType;

      public:
        typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

        explicit FemCommunication ( const DiscreteFunctionSpaceType &dfSpace,  Dune::SolverCategory::Category solverCategory = Dune::SolverCategory::sequential )
          : dfSpace_( dfSpace ), solverCategory_( solverCategory )
        {}

        const typename DiscreteFunctionSpace::GridPartType::CommunicationType &communicator () const { return dfSpace_.gridPart().comm(); }

        template< class T >
        void copyOwnerToAll ( const T &x, T &y ) const
        {
          y = x;
          project( y );
          FemCommunicationVector< DiscreteFunctionSpaceType, T > z( dfSpace_, y );
          dfSpace_.communicator().exchange( z, DFCommunicationOperation::Add() );
        }

        template< class T >
        void project ( T &x ) const
        {
          typedef typename T::field_type field_type;

          // clear auxiliary DoFs
          for( int i : dfSpace_.auxiliaryDofs() )
            x[ i ] = field_type( 0 );
        }

        template< class T, class F >
        void dot ( const T &x, const T &y, F &scp ) const
        {
          const auto &auxiliaryDofs = dfSpace_.auxiliaryDofs();

          const int numAuxiliarys = auxiliaryDofs.size();
          for( int auxiliary = 0, i = 0; auxiliary < numAuxiliarys; ++auxiliary, ++i )
          {
            const int nextAuxiliary = auxiliaryDofs[ auxiliary ];
            for( ; i < nextAuxiliary; ++i )
              scp += x[ i ] * y[ i ];
          }

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
        const DiscreteFunctionSpaceType &dfSpace_;
        Dune::SolverCategory::Category solverCategory_;
      };



      // buildCommunication
      // ------------------

      template< class DiscreteFunctionSpace >
      void buildCommunication ( const DiscreteFunctionSpace &dfSpace,
                                Dune::SolverCategory::Category solverCategory,
                                std::shared_ptr< FemCommunication< DiscreteFunctionSpace > > &communication )
      {
        communication.reset( new FemCommunication< DiscreteFunctionSpace >( dfSpace, solverCategory ) );
      }



      // SupportsAMG for FemCommunication
      // --------------------------------

      template< class T >
      struct SupportsAMG;

      template< class DiscreteFunctionSpace >
      struct SupportsAMG< FemCommunication< DiscreteFunctionSpace > >
        : public std::false_type
      {};

    } // namespace ISTL

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SOLVER_COMMUNICATION_FEM_HH
