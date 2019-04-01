#ifndef DUNE_FEM_LFESPACE_LOCALRESTRICTPROLONG_HH
#define DUNE_FEM_LFESPACE_LOCALRESTRICTPROLONG_HH

#include <dune/common/dynmatrix.hh>
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>

namespace Dune
{

  namespace Fem
  {

    template< class LFEMap, class FunctionSpace, template< class > class Storage >
    struct DefaultLocalRestrictProlong< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
    {
      typedef LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > DiscreteFunctionSpaceType;
      typedef DefaultLocalRestrictProlong< DiscreteFunctionSpaceType > ThisType;

      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;

      DefaultLocalRestrictProlong (const DiscreteFunctionSpaceType &space)
      : space_( space ), weight_(-1)
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
        auto P = calcProlongMatrix(lfFather.entity(), lfSon.entity(), numDofs);

        if( initialize )
          for( int i = 0; i < numDofs; ++i )
            lfFather[ i ] = 0;
        // y += alpha A^T x
        P.usmtv(weight,lfSon.localDofVector(),lfFather.localDofVector());
      }

      //! prolong data to children
      template< class LFFather, class LFSon, class LocalGeometry >
      void prolongLocal ( const LFFather &lfFather, LFSon &lfSon,
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {
        const int numDofs = lfFather.numDofs();
        assert( lfFather.numDofs() == lfSon.numDofs() );
        auto P = calcProlongMatrix(lfFather.entity(), lfSon.entity(), numDofs);
        for( int i = 0; i < numDofs; ++i )
          P.mv(lfFather.localDofVector(),lfSon.localDofVector());
      }

      //! do discrete functions need a communication after restriction / prolongation?
      bool needCommunication () const { return true; }

    protected:
      template< class Entity >
      static DomainFieldType calcWeight ( const Entity &father, const Entity &son )
      {
        return son.geometry().volume() / father.geometry().volume();
      }
      template< class Entity >
      DynamicMatrix<double> calcProlongMatrix( const Entity &father, const Entity &son, int numDofs ) const
      {
        // TODO
        DynamicMatrix<double> P(numDofs,numDofs);
        return P;
      }

      const DiscreteFunctionSpaceType &space_;
      DomainFieldType weight_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LOCALRESTRICTPROLONG_HH
