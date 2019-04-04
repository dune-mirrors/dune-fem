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
          lfFather_.evaluate( localGeo_.global(x), y);
        }

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

        SonsWrapper (
          const std::vector< EntityType >& childEntities,
          const std::vector< BasisFunctionSetType >& childBasisSets,
          const std::vector< std::vector<double> >& childDofs )
          : childEntities_(childEntities), childBasisSets_(childBasisSets), childDofs_(childDofs)
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
            auto y = childEntities_[i].geometryInFather().local(x);
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

      private:
        const std::vector< EntityType >& childEntities_;
        const std::vector< BasisFunctionSetType >& childBasisSets_;
        const std::vector< std::vector<double> >& childDofs_;
      };
    } // namespce Impl

    // a detailed description is given in MR308
    // https://gitlab.dune-project.org/dune-fem/dune-fem/merge_requests/308
    template< class LFEMap, class FunctionSpace, template< class > class Storage >
    struct DefaultLocalRestrictProlong< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
    {
      typedef LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > DiscreteFunctionSpaceType;
      typedef DefaultLocalRestrictProlong< DiscreteFunctionSpaceType > ThisType;

      typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
      typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename DiscreteFunctionSpaceType::EntityType EntityType;
      typedef typename EntityType::LocalGeometry LocalGeometryType;
      typedef typename EntityType::EntitySeed EntitySeedType;

      DefaultLocalRestrictProlong (const DiscreteFunctionSpaceType &space)
      : space_( space ), childSeeds_(0), childDofs_(0)
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
        const int numDofs = lfFather.numDofs();
        assert( lfFather.numDofs() == lfSon.numDofs() );

        if (initialize)
        {
          childSeeds_.resize(0);
          childDofs_.resize(0);
        }

        childSeeds_.push_back(lfSon.entity().seed());
        childDofs_.push_back(std::vector<double>(lfSon.size()));
        for (unsigned int i=0;i<lfSon.size();++i)
          childDofs_.back()[i] = lfSon[i];
      }

      template <class LFFather>
      void restrictFinalize( LFFather &lfFather ) const
      {
        const int numDofs = lfFather.numDofs();
        std::vector< EntityType > childEntities(childSeeds_.size());
        std::vector< BasisFunctionSetType > childBasisSets(childSeeds_.size());
        for (unsigned int i=0; i<childSeeds_.size();++i)
        {
          childEntities[i] = space_.gridPart().entity( childSeeds_[i] );
          childBasisSets[i] = space_.basisFunctionSet( childEntities[i] );
        }
        space_.interpolation(lfFather.entity())
          ( Impl::SonsWrapper<BasisFunctionSetType, LFFather>( childEntities, childBasisSets, childDofs_ ),
            lfFather.localDofVector() );
      }

      //! prolong data to children
      template< class LFFather, class LFSon, class LocalGeometry >
      void prolongLocal ( const LFFather &lfFather, LFSon &lfSon,
                          const LocalGeometry &geometryInFather, bool initialize ) const
      {
        const int numDofs = lfFather.numDofs();
        assert( lfFather.numDofs() == lfSon.numDofs() );
        DynamicVector<double> sonDofs( numDofs );
        space_.interpolation(lfSon.entity())
          ( Impl::FatherWrapper<LocalGeometry,LFFather>(geometryInFather,lfFather),
            sonDofs );
        for (int i=0;i<numDofs;++i)
          lfSon[i] = sonDofs[i];
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
      mutable std::vector< EntitySeedType > childSeeds_;
      mutable std::vector< std::vector<double> > childDofs_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LOCALRESTRICTPROLONG_HH
