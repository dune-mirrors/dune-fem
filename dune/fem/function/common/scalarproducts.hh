#ifndef DUNE_FEM_SCALARPRODURCTS_HH
#define DUNE_FEM_SCALARPRODURCTS_HH

#include <iostream>
#include <memory>
#include <set>
#include <map>
#include <limits>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/genericiterator.hh>
#include <dune/common/ftraits.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/scalarproducts.hh>
#endif

#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/slavedofs.hh>
#include <dune/fem/space/common/commindexmap.hh>
#include <dune/fem/function/blockvectorfunction/declaration.hh>
#include <dune/fem/function/blockvectors/defaultblockvectors.hh>

namespace Dune
{

  namespace Fem
  {

  /** @addtogroup Communication Communication
      @{
  **/

#if HAVE_DUNE_ISTL
    template <class DofVector>
    struct ISTLScalarProductSelector
    {
      typedef Dune::FieldVector< typename DofVector::FieldType, DofVector::blockSize > Block;
      typedef Dune::BlockVector< Block > type;
    };

    template <class Block>
    struct ISTLScalarProductSelector< Dune::Fem::ISTLBlockVector< Block > >
      : public Dune::ScalarProduct< typename Dune::Fem::ISTLBlockVector< Block > :: DofContainerType >
    {
      //! define the category
      enum { category=Dune::SolverCategory::sequential };

      typedef typename ISTLBlockVector< Block > :: DofContainerType type;
    };
#endif

    //! Proxy class to evaluate ScalarProduct
    //! holding SlaveDofs which is singleton per space and mapper
    template< class DiscreteFunction >
    class ParallelScalarProduct
#if HAVE_DUNE_ISTL
      : public ISTLScalarProductSelector< typename DiscreteFunction :: DofVectorType >
#endif
    {
    public:
      typedef DiscreteFunction DiscreteFunctionType;

      //! type of the discrete function space
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;

    private:
      typedef ParallelScalarProduct< DiscreteFunctionType > ThisType;

    public:
      //! type of range field
      typedef typename DiscreteFunctionSpaceType :: RangeFieldType  RangeFieldType;

      //! type of used mapper
      typedef typename DiscreteFunctionSpaceType :: BlockMapperType MapperType;

      enum { blockSize = DiscreteFunctionSpaceType :: localBlockSize };

      // type of communication manager object which does communication
      typedef SlaveDofs< DiscreteFunctionSpaceType, MapperType > SlaveDofsType;

      typedef RangeFieldType  field_type;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type real_type;

      //! constructor taking space
      ParallelScalarProduct ( const DiscreteFunctionSpaceType &space )
      : space_( space )
      {}

      const DiscreteFunctionSpaceType &space() const
      {
        return space_;
      }

      //! evaluate scalar product and omit slave nodes
      template < class OtherDiscreteFunctionType >
      RangeFieldType scalarProductDofs ( const DiscreteFunctionType &x, const OtherDiscreteFunctionType &y ) const
      {
        assert(x.space() == y.space());
        assert(x.space() == space());
        return dotProduct( x.dofVector(), y.dofVector() );
      }

      const SlaveDofsType &slaveDofs() const
      {
        return space().slaveDofs();
      }

    protected:
      //! evaluate scalar product on dofVector and omit slave nodes
      template < class DofVector, class OtherDofVector >
      RangeFieldType dotProduct ( const DofVector &x, const OtherDofVector &y ) const
      {
        auto &slaveDofs = space().slaveDofs();

        RangeFieldType scp = 0;

        const int numSlaves = slaveDofs.size();
        for( int slave = 0, i = 0 ; slave < numSlaves; ++slave )
        {
          const int nextSlave = slaveDofs[ slave ];
          for(; i < nextSlave; ++i )
            for( unsigned int j = 0; j < blockSize; ++j )
              scp += x[ i ][ j ] * y[ i ][ j ];

          // skip the slave dof
          ++i;
        }

        // do global sum
        scp = space().gridPart().comm().sum( scp );
        return scp;
      }

#if HAVE_DUNE_ISTL
    protected:
      typedef typename ISTLScalarProductSelector< typename DiscreteFunction :: DofVectorType > :: type BlockVectorType;

    public:
      //! dot product for ISTL solvers
      virtual field_type dot (const BlockVectorType& x,
                              const BlockVectorType& y)
      {
        return dotProduct( x, y );
      }

      //! norm for ISTL solvers
      virtual real_type norm( const BlockVectorType& x )
      {
        return std::abs( std::sqrt( dotProduct( x, x ) ) );
      }

      //! delete slave values (for debugging)
      void deleteNonInterior( BlockVectorType& x) const
      {
#if HAVE_MPI
        // case of ALUGrid and DGSpace or FVSpace
        const bool deleteGhostEntries = (space().gridPart().grid().overlapSize( 0 ) == 0) && !space().continuous();

        // only delete ghost entries
        if( deleteGhostEntries )
        {
          const auto &slaveDofs = space().slaveDofs();

          // don't delete the last since this is the overall Size
          const int slaveSize = slaveDofs.size() - 1;
          for(int slave = 0; slave<slaveSize; ++slave)
            x[ slaveDofs[slave] ] = 0;
        }
#endif
      }
#endif
      const DiscreteFunctionSpaceType &space_;
    };

  //@}

  } // end namespace Fem

} // end namespace Dune
#endif // #ifndef DUNE_FEM_SCALARPRODURCTS_HH
