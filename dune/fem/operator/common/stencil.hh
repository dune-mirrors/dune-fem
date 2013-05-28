// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_STENCIL_HH
#define DUNE_FEM_STENCIL_HH

#include <iostream>
#include <set>
#include <map>

namespace Dune 
{
  namespace Fem 
  {
    template <class DomainSpace, class RangeSpace>
    struct Stencil
    {
      // Domain = Row
      typedef typename DomainSpace::IteratorType        DomainIteratorType;
      typedef typename DomainIteratorType::Entity       DomainEntityType;
      typedef typename DomainSpace::BlockMapperType     DomainBlockMapper;
      typedef typename DomainBlockMapper::GlobalKeyType DomainGlobalKeyType;

      // Range = Column
      typedef typename RangeSpace::IteratorType         RangeIteratorType;
      typedef typename RangeIteratorType::Entity        RangeEntityType;
      typedef typename RangeSpace::BlockMapperType      RangeBlockMapper;
      typedef typename RangeBlockMapper::GlobalKeyType  RangeGlobalKeyType;

      typedef std::set<RangeGlobalKeyType> LocalStencilType;
      typedef std::map<DomainGlobalKeyType,LocalStencilType> GlobalStencilType;

    private:

      struct FillFunctor
      {
        typedef DomainGlobalKeyType GlobalKey;
        FillFunctor(GlobalStencilType &stencil) 
        : stencil_(stencil),
          localStencil_(0)
        {}
        void set(const std::size_t, const DomainGlobalKeyType &domainGlobal)
        {
          localStencil_ = &(stencil_[ domainGlobal ]);
        }
        void operator() ( const std::size_t, const RangeGlobalKeyType &rangeGlobal)
        {
          localStencil_->insert( rangeGlobal );
        }
        private:
        GlobalStencilType &stencil_;
        LocalStencilType *localStencil_;
      };

      typedef typename Dune::Fem::MatrixFunctor<RangeBlockMapper,RangeEntityType,FillFunctor > MFunctor;

    public:
      Stencil(const DomainSpace &dSpace, const RangeSpace &rSpace)
        : domainBlockMapper_( dSpace.blockMapper() )
        , rangeBlockMapper_( rSpace.blockMapper() )
      {
      }

      //! create entries for element and neighbors
      void fill ( const DomainEntityType &dEntity, const RangeEntityType &rEntity )
      {
        domainBlockMapper_.mapEach(dEntity, 
                  MFunctor( rangeBlockMapper_, rEntity, FillFunctor(globalStencil_) ) );
      }

      const LocalStencilType &localStencil(const DomainGlobalKeyType &key) const
      { 
        return globalStencil_[key]; 
      }
      const GlobalStencilType &globalStencil() const
      { 
        return globalStencil_; 
      }
      int maxNZ() const
      {
        int ret = 0;
        typedef typename GlobalStencilType::const_iterator StencilIteratorType;
        const GlobalStencilType &glStencil = globalStencil();
        StencilIteratorType end = glStencil.end();
        for ( StencilIteratorType it = glStencil.begin(); it != end; ++it)
          ret = std::max(ret,it->second.size());
        return ret;
      }

    private:
      const DomainBlockMapper &domainBlockMapper_;
      const RangeBlockMapper &rangeBlockMapper_;
      GlobalStencilType globalStencil_;
    };

    template <class DomainSpace, class RangeSpace>
    struct SimpleStencil : public Stencil<DomainSpace,RangeSpace>
    {
      typedef Stencil<DomainSpace,RangeSpace> BaseType;
      typedef typename DomainSpace::IteratorType        DomainIteratorType;
      typedef typename DomainIteratorType::Entity       DomainEntityType;
      typedef typename DomainSpace::BlockMapperType     DomainBlockMapper;
      typedef typename DomainBlockMapper::GlobalKeyType DomainGlobalKeyType;

      typedef typename RangeSpace::IteratorType         RangeIteratorType;
      typedef typename RangeIteratorType::Entity        RangeEntityType;
      typedef typename RangeSpace::BlockMapperType      RangeBlockMapper;
      typedef typename RangeBlockMapper::GlobalKeyType  RangeGlobalKeyType;

      typedef std::set<RangeGlobalKeyType> LocalStencilType;
      typedef std::map<DomainGlobalKeyType,LocalStencilType> GlobalStencilType;

      SimpleStencil(int maxNZ) 
      : maxNZ_(maxNZ)
      {}
      int maxNZ() const
      {
        return maxNZ_;
      }
      void fill ( const DomainEntityType &dEntity, const RangeEntityType &rEntity )
      {
      }
      const LocalStencilType &localStencil(const DomainGlobalKeyType &key) const
      { 
      }
      const GlobalStencilType &globalStencil() const
      { 
      }
    private:
      int maxNZ_;
    };
  
  } // namespace Fem

} // namespace Dune

#endif // #if defined DUNE_FEM_STENCIL_HH

