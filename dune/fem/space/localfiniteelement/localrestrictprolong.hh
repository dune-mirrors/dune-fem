#ifndef DUNE_FEM_LFESPACE_LOCALRESTRICTPROLONG_HH
#define DUNE_FEM_LFESPACE_LOCALRESTRICTPROLONG_HH

#include <dune/common/dynmatrix.hh>
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Impl
    {
      // BasisFunctionWrapper
      // --------------------
      template< class BasisFunctionSet, class LocalGeometry >
      struct BasisFunctionWrapper
      {
        typedef std::remove_reference_t< BasisFunctionSet > BasisFunctionSetType;
        typedef std::remove_reference_t< LocalGeometry > LocalGeometryType;
        struct Traits
        {
          typedef typename BasisFunctionSetType::DomainType DomainType;
          typedef typename BasisFunctionSetType::RangeType RangeType;
        };
        typedef typename BasisFunctionSetType::EntityType EntityType;
        typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;
        typedef typename BasisFunctionSetType::DomainType DomainType;
        typedef typename BasisFunctionSetType::RangeType RangeType;
        typedef typename BasisFunctionSetType::JacobianRangeType JacobianRangeType;
        typedef typename BasisFunctionSetType::HessianRangeType HessianRangeType;

        BasisFunctionWrapper ( const LocalGeometryType &localGeo, const BasisFunctionSetType &bset, int i )
          : localGeo_( localGeo ), bset_( bset ), i_( i ), values_( bset.size() ) {}

        template< class Arg >
        void evaluate ( const Arg &x, typename Traits::RangeType &y ) const
        {
          bset_.evaluateAll( localGeo_.global(x), values_);
          y = values_[i_];
        }

      private:
        const LocalGeometryType &localGeo_;
        const BasisFunctionSetType &bset_;
        const int i_;
        mutable std::vector<RangeType> values_;
      };
    } // namespce Impl

    template< class LFEMap, class FunctionSpace, template< class > class Storage >
    struct DefaultLocalRestrictProlong< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
    {
      typedef LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > DiscreteFunctionSpaceType;
      typedef DefaultLocalRestrictProlong< DiscreteFunctionSpaceType > ThisType;

      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;

      DefaultLocalRestrictProlong (const DiscreteFunctionSpaceType &space)
      : space_( space ), weight_(-1),
        lfFather_(0)
      {}

      /** \brief explicit set volume ratio of son and father
       *
       *  \param[in]  weight  volume of son / volume of father
       *
       *  \note If this ratio is set, it is assume to be constant.
       */
      void setFatherChildWeight ( const DomainFieldType &weight )
      {
        weight_ = weight;
      }

      //! restrict data to father
      template< class LFFather, class LFSon, class LocalGeometry >
      void restrictLocal ( LFFather &lfFather, const LFSon &lfSon,
                           const LocalGeometry &geometryInFather, bool initialize ) const
      {
        const DomainFieldType weight = (weight_ < DomainFieldType( 0 ) ? calcWeight( lfFather.entity(), lfSon.entity() ) : weight_);

        assert( weight > 0.0 );

        const int numDofs = lfFather.numDofs();
        assert( lfFather.numDofs() == lfSon.numDofs() );

        if (initialize)
          lfFather_.resize(numDofs,0);

        std::cout << "initialize=" << initialize;
        for (int r=0;r<numDofs;++r)
          std::cout << "   " << lfFather_[r] << " , " << lfSon[r];

        auto P = calcProlongMatrix(lfFather.entity(), lfSon.entity(), numDofs);
        P.usmv(weight, lfSon.localDofVector(), lfFather_);

        std::cout << "   ->    ";
        for (int r=0;r<numDofs;++r)
          std::cout << "   " << lfFather_[r];
        std::cout << std::endl;
      }

      template <class LFFather>
      void restrictFinalize( LFFather &lfFather ) const
      {
        const int numDofs = lfFather.numDofs();
        assert( numDofs == lfFather_.size() );
        for (int r=0;r<numDofs;++r)
          lfFather[r] = lfFather_[r];
        std::cout << "finalize";
        for (int r=0;r<numDofs;++r)
          std::cout << "   " << lfFather[r];
        std::cout << std::endl << std::endl;
        lfFather_.clear();
      }

      //! prolong data to children
      template< class LFFather, class LFSon, class LocalGeometry >
      void prolongLocal ( const LFFather &lfFather, LFSon &lfSon,
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {
        const int numDofs = lfFather.numDofs();
        assert( lfFather.numDofs() == lfSon.numDofs() );
        auto P = calcProlongMatrix(lfFather.entity(), lfSon.entity(), numDofs);

        // Note: we want lfSon = P^T lfFather
        // but lfSon and lfFather might share dofs (e.g. Lagrange space)
        // and then the necessary lfSon.clear() before computing the matrix
        // vector multiplication would lead to errors so we first need to
        // copy the lfFather dofs:
        DynamicVector<double> lfCopy(numDofs);
        for (int r=0;r<numDofs;++r)
          lfCopy[r] = lfFather[r];
        // y = P^T x
        P.mtv(lfCopy,lfSon.localDofVector());
      }

      //! do discrete functions need a communication after restriction / prolongation?
      bool needCommunication () const { return true; }

    protected:
      template< class Entity >
      static DomainFieldType calcWeight ( const Entity &father, const Entity &son )
      {
        return son.geometry().volume() / father.geometry().volume();
      }
      // note that this method returns
      // P = ( l(phi_r)_c )_rc so it is the transposed of the matrix needed
      template< class Entity >
      DynamicMatrix<double> calcProlongMatrix( const Entity &father, const Entity &son, int numDofs ) const
      {
        DynamicMatrix<double> P(numDofs,numDofs);
        const auto &localBasis = space_.basisFunctionSet(father);
        const auto &localInterpolation = space_.interpolation(son);
        const auto &localGeo = son.geometryInFather();
        for (int r=0;r<numDofs;++r)
        {
          Impl::BasisFunctionWrapper< decltype(localBasis), decltype(localGeo) >
              basisWrapper( localGeo, localBasis, r );
          localInterpolation( basisWrapper, P[r] );
        }
        /*
        for (int r=0;r<numDofs;++r)
        {
          for (int c=0;c<numDofs;++c)
            std::cout << P[c][r] << " ";
          std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
        */
        return P;
      }
      const DiscreteFunctionSpaceType &space_;
      DomainFieldType weight_;
      mutable std::vector<double> lfFather_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LOCALRESTRICTPROLONG_HH
