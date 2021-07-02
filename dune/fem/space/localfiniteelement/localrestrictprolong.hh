#ifndef DUNE_FEM_LFESPACE_LOCALRESTRICTPROLONG_HH
#define DUNE_FEM_LFESPACE_LOCALRESTRICTPROLONG_HH

#include <dune/common/dynmatrix.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/space/localfiniteelement/space.hh>
#include <dune/fem/space/common/localinterpolation.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Impl
    {
      template< class LocalGeometry, class LF>
      struct FatherWrapper
      {
        typedef std::remove_reference_t< LF > LocalFunctionType;
        typedef std::remove_reference_t< LocalGeometry > LocalGeometryType;
        struct Traits
        {
          typedef typename LocalFunctionType::DomainType DomainType;
          typedef typename LocalFunctionType::RangeType RangeType;
        };
        typedef typename LocalFunctionType::EntityType EntityType;
        typedef typename LocalFunctionType::FunctionSpaceType FunctionSpaceType;
        typedef typename LocalFunctionType::DomainType DomainType;
        typedef typename LocalFunctionType::RangeType RangeType;
        typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;
        typedef typename LocalFunctionType::HessianRangeType HessianRangeType;

        FatherWrapper ( const LocalGeometryType &localGeo, const LocalFunctionType &lfFather)
          : localGeo_(localGeo), lfFather_(lfFather) {}
        template <class Point>
        void evaluate ( const Point &x, RangeType &y ) const
        {
          lfFather_.evaluate( localGeo_.global(coordinate(x)), y);
        }
        template <class Quadrature, class RangeArray>
        void evaluateQuadrature( const Quadrature& quadrature, RangeArray& values ) const
        {
          const unsigned int nop = quadrature.nop();
          values.resize( nop );
          for( unsigned int qp=0; qp<nop; ++qp)
            evaluate( quadrature[ qp ], values[ qp ]);
        }

        const EntityType& entity () const { return lfFather_.entity(); }
        const unsigned int order () const { return lfFather_.order(); }

      private:
        const LocalGeometryType &localGeo_;
        const LocalFunctionType &lfFather_;
      };
      template< class BasisFunctionSet, class LF>
      struct SonsWrapper
      {
        typedef std::remove_reference_t< LF > LocalFunctionType;
        typedef std::remove_reference_t< BasisFunctionSet > BasisFunctionSetType;
        struct Traits
        {
          typedef typename LocalFunctionType::DomainType DomainType;
          typedef typename LocalFunctionType::RangeType RangeType;
        };
        typedef typename LocalFunctionType::FunctionSpaceType FunctionSpaceType;
        typedef typename LocalFunctionType::DomainType DomainType;
        typedef typename LocalFunctionType::RangeType RangeType;
        typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;
        typedef typename LocalFunctionType::HessianRangeType HessianRangeType;
        typedef typename LocalFunctionType::EntityType EntityType;
        typedef typename EntityType::LocalGeometry LocalGeometryType;

        template <class LFFather>
        SonsWrapper ( const LFFather &father,
          const std::vector< EntityType >& childEntities,
          const std::vector< BasisFunctionSetType >& childBasisSets,
          const std::vector< std::vector<double> >& childDofs )
          : father_(father.entity())
          , order_(father.order())
          , childEntities_(childEntities), childBasisSets_(childBasisSets), childDofs_(childDofs)
        {}
        template <class Point>
        void evaluate ( const Point &x, RangeType &val ) const
        {
          val = RangeType(0);
          RangeType tmp;
          double weight = 0;
          for (unsigned int i=0; i<childEntities_.size();++i)
          {
            const auto &refSon = Dune::ReferenceElements< typename LocalGeometryType::ctype, LocalGeometryType::mydimension >
              ::general( childEntities_[i].type() );
            auto y = childEntities_[i].geometryInFather().local(coordinate(x));
            if( refSon.checkInside( y ) )
            {
              childBasisSets_[i].evaluateAll( y, childDofs_[i], tmp );
              val += tmp;
              weight += 1.;
            }
          }
          assert( weight > 0); // weight==0 would mean that point was found in none of the children
          val /= weight;
        }
        template <class Quadrature, class RangeArray>
        void evaluateQuadrature( const Quadrature& quadrature, RangeArray& values ) const
        {
          const unsigned int nop = quadrature.nop();
          values.resize( nop );
          for( unsigned int qp=0; qp<nop; ++qp)
            evaluate( quadrature[ qp ], values[ qp ]);
        }
        const EntityType& entity () const { return father_; }
        const unsigned int order () const { return order_; }

      private:
        const EntityType &father_;
        unsigned int order_;
        const std::vector< EntityType >& childEntities_;
        const std::vector< BasisFunctionSetType >& childBasisSets_;
        const std::vector< std::vector<double> >& childDofs_;
      };

      // a detailed description is given in MR308
      // https://gitlab.dune-project.org/dune-fem/dune-fem/merge_requests/308
      template< class LFESpace >
      struct DefaultLocalRestrictProlongLFE
      {
        typedef LFESpace DiscreteFunctionSpaceType;
        typedef DefaultLocalRestrictProlongLFE< DiscreteFunctionSpaceType > ThisType;

        typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
        typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
        typedef typename DiscreteFunctionSpaceType::EntityType EntityType;
        typedef typename EntityType::LocalGeometry LocalGeometryType;
        typedef typename EntityType::EntitySeed EntitySeedType;

        DefaultLocalRestrictProlongLFE (const DiscreteFunctionSpaceType &space)
        : space_( space ),
          interpolation_( space_ ),
          childSeeds_(0), childDofs_(0)
        {}

        /** \brief explicit set volume ratio of son and father
         *
         *  \param[in]  weight  volume of son / volume of father
         *
         *  \note If this ratio is set, it is assume to be constant.
         */
        void setFatherChildWeight ( const DomainFieldType &weight )
        {
        }

        //! restrict data to father
        template< class LFFather, class LFSon, class LocalGeometry >
        void restrictLocal ( LFFather &lfFather, const LFSon &lfSon,
                             const LocalGeometry &geometryInFather, bool initialize ) const
        {
          assert( lfFather.numDofs() == lfSon.numDofs() );

          if (initialize)
          {
            childSeeds_.resize(0);
            childDofs_.resize(0);
          }

          childSeeds_.push_back(lfSon.entity().seed());
          childDofs_.push_back(std::vector<double>(lfSon.size()));
          std::vector<double>& chDofs = childDofs_.back();
          const unsigned int numDofs = lfSon.size();
          for (unsigned int i=0;i<numDofs; ++i )
            chDofs[i] = lfSon[i];
        }

        template <class LFFather>
        void restrictFinalize( LFFather &lfFather ) const
        {
          std::vector< EntityType > childEntities(childSeeds_.size());
          std::vector< BasisFunctionSetType > childBasisSets(childSeeds_.size());
          const unsigned int cSize = childSeeds_.size();
          for (unsigned int i=0; i<cSize; ++i)
          {
            childEntities[i] = space_.gridPart().entity( childSeeds_[i] );
            childBasisSets[i] = space_.basisFunctionSet( childEntities[i] );
          }


          // bind interpolation object to father
          auto guard = bindGuard( interpolation_, lfFather.entity() );

          typedef Impl::SonsWrapper<BasisFunctionSetType, LFFather> SonsWrapperType;
          SonsWrapperType sonsWrapper( lfFather, childEntities, childBasisSets, childDofs_ );

          // call interpolation
          interpolation_( sonsWrapper, lfFather );
        }

        //! prolong data to children
        template< class LFFather, class LFSon, class LocalGeometry >
        void prolongLocal ( const LFFather &lfFather, LFSon &lfSon,
                            const LocalGeometry &geometryInFather, bool initialize ) const
        {
          assert( lfFather.numDofs() == lfSon.numDofs() );

          typedef Impl::FatherWrapper<LocalGeometry,LFFather> FatherWrapperType;
          FatherWrapperType fatherWrapper(geometryInFather,lfFather);

          // bind interpolation object to father
          auto guard = bindGuard( interpolation_, lfSon.entity() );

          // call interpolation
          interpolation_( fatherWrapper, lfSon );
        }

        //! do discrete functions need a communication after restriction / prolongation?
        bool needCommunication () const { return true; }

      protected:
        template< class Entity >
        static DomainFieldType calcWeight ( const Entity &father, const Entity &son )
        {
          return son.geometry().volume() / father.geometry().volume();
        }
        const DiscreteFunctionSpaceType &space_;
        mutable LocalInterpolation< DiscreteFunctionSpaceType > interpolation_;

        mutable std::vector< EntitySeedType > childSeeds_;
        mutable std::vector< std::vector<double> > childDofs_;

        mutable std::vector< double > tmpDofs_;
      };
    } // namespce Impl

    template< class LFEMap, class FunctionSpace, class Storage >
    class DefaultLocalRestrictProlong< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
    : public Impl::DefaultLocalRestrictProlongLFE
               < LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
    {
    public:
      typedef Impl::DefaultLocalRestrictProlongLFE
               < LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > > BaseType;
      using BaseType::BaseType;
    };
    template< class LFEMap, class FunctionSpace, class Storage >
    class DefaultLocalRestrictProlong< DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
    : public Impl::DefaultLocalRestrictProlongLFE
               < DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >

    {
    public:
      typedef Impl::DefaultLocalRestrictProlongLFE
               < DiscontinuousLocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > > BaseType;
      using BaseType::BaseType;
    };
  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LOCALRESTRICTPROLONG_HH
