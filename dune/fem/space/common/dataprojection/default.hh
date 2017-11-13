#ifndef DUNE_FEM_HPDG_SPACE_COMMON_DATAPROJECTION_DEFAULT_HH
#define DUNE_FEM_HPDG_SPACE_COMMON_DATAPROJECTION_DEFAULT_HH

#include <cassert>
#include <cstddef>

#include <functional>
#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/fem/function/localfunction/localfunction.hh>

#include "dataprojection.hh"

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // DefaultDataProjection
      // ---------------------

      /** \brief Local \f$L^2(\Omega)\f$-projection for the restriction and prolongation of discrete functions
       *
       *  \tparam  DiscreteFunction  type of the discrete function
       *
       *  \ingroup DiscreteFunctionSpace_RestrictProlong
       */
      template< class DiscreteFunction >
      class DefaultDataProjection
      : public DataProjection< typename DiscreteFunction::DiscreteFunctionSpaceType, DefaultDataProjection< DiscreteFunction > >
      {
        using ThisType = DefaultDataProjection< DiscreteFunction >;
        using BaseType = DataProjection< typename DiscreteFunction::DiscreteFunctionSpaceType, DefaultDataProjection< DiscreteFunction > >;

      public:
        /** \copydoc Dune::Fem::hpDG::DataProjection::DiscreteFunctionSpaceType */
        using DiscreteFunctionSpaceType = typename BaseType::DiscreteFunctionSpaceType;
        /** \copydoc Dune::Fem::hpDG::DataProjection::BasisFunctionSetType */
        using BasisFunctionSetType = typename BaseType::BasisFunctionSetType;
        /** \copydoc Dune::Fem::hpDG::DataProjection::EntityType */
        using EntityType = typename BaseType::EntityType;

      private:
        using RangeFieldType = typename DiscreteFunction::RangeFieldType;
        using LocalDofVectorType = Dune::DynamicVector< RangeFieldType >;
        using LocalFunctionType = Dune::Fem::LocalFunction< BasisFunctionSetType, LocalDofVectorType >;

        static const std::size_t localBlockSize = DiscreteFunctionSpaceType::localBlockSize;

      public:
        /** \name Construction
         *  \{
         */

        explicit DefaultDataProjection ( DiscreteFunction &discreteFunction )
          : discreteFunction_( discreteFunction )
        {
          discreteFunction_.get().enableDofCompression();
        }

        /** \} */

#ifndef DOXYGEN

        DefaultDataProjection ( const ThisType & ) = delete;

        DefaultDataProjection ( ThisType && ) = default;

        ThisType &operator= ( const ThisType & ) = delete;

        ThisType &operator= ( ThisType && ) = default;

#endif // #ifndef DOXYGEN

        /** \copydoc Dune::Fem::hpDG::DataProjection::operator() */
        void operator() ( const EntityType &entity,
                          const BasisFunctionSetType &prior,
                          const BasisFunctionSetType &present,
                          const std::vector< std::size_t > &origin,
                          const std::vector< std::size_t > &destination )
        {
          LocalFunctionType localFunction( prior );
          read( origin, localFunction.localDofVector() );

          assert( present.size() == space().basisFunctionSet( entity ).size() );
          LocalDofVectorType localDofVector( present.size() );

          const auto interpolation = space().interpolation( entity );
          interpolation( localFunction, localDofVector );

          write( destination, localDofVector );
        }

        /** \brief transfer of discrete function from old to new space using intermediate storage */
        template <class TemporaryStorage>
        void operator () ( TemporaryStorage& tmp )
        {
          auto& df = discreteFunction();

          // copy dofs to temporary storage, old order still exists
          auto dfit = df.dbegin();
          const auto endtmp = tmp.dend();
          for( auto it = tmp.dbegin(); it != endtmp; ++it, ++dfit )
          {
            assert( dfit != df.dend() );
            *it = *dfit;
          }

          // interpolate to new space, this can be a
          // Lagrange interpolation or a L2 projection, both locally
          const auto end = df.space().end();
          for( auto it = df.space().begin(); it != end; ++it )
          {
            const auto& entity = *it;
            const auto tmpLF = tmp.localFunction( entity );
            auto lf = df.localFunction( entity );

            // if size is the same we can just copy the dof values
            if( tmpLF.size() == lf.size() )
            {
              lf.assign( tmpLF );
            }
            else
            {
              // otherwise a local interpolation is needed
              df.space().interpolation( entity )( tmpLF, lf );
            }
          }
        }

        /** \copydoc Dune::Fem::Adaptive::DataProjection::addToList () */
        template< class Communicator >
        void addToList ( Communicator &comm )
        {
          comm.addToList( discreteFunction() );
        }

      private:
        template< class LocalDofVector >
        void read ( const std::vector< std::size_t > &blocks, LocalDofVector &localDofVector ) const
        {
          assert( localDofVector.size() == localBlockSize*blocks.size() );
          std::size_t index = 0;
          for( auto i : blocks )
          {
            const auto block = discreteFunction().block( i );
            for( std::size_t j = 0; j < localBlockSize; ++j )
              localDofVector[ index++ ] = (*block)[ j ];
          }
          assert( index == localDofVector.size() );
        }

        template< class LocalDofVector >
        void write ( const std::vector< std::size_t > &blocks, const LocalDofVector &localDofVector )
        {
          assert( localDofVector.size() == localBlockSize*blocks.size() );
          std::size_t index = 0;
          for( auto i : blocks )
          {
            auto block = discreteFunction().block( i );
            for( std::size_t j = 0; j < localBlockSize; ++j )
              (*block)[ j ] = localDofVector[ index++ ];
          }
          assert( index == localDofVector.size() );
        }

        DiscreteFunction &discreteFunction () { return discreteFunction_.get(); }

        const DiscreteFunction &discreteFunction () const { return discreteFunction_.get(); }

        const DiscreteFunctionSpaceType &space () const { return discreteFunction().space(); }

        std::reference_wrapper< DiscreteFunction > discreteFunction_;
      };

    } // namespace hpDG

    // forward types to Fem namespace for convenience
    using hpDG::DefaultDataProjection;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_SPACE_COMMON_DATAPROJECTION_DEFAULT_HH
