#ifndef DUNE_SUBSPACE_HH
#define DUNE_SUBSPACE_HH

//- System includes

//- Dune includes

//- Local includes
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/common/dofmapperinterface.hh>
#include <dune/fem/space/common/dofstorage.hh>

namespace Dune {

  //- Forward declarations
  template <class CombinedSpaceImp>
  class SubSpace;
  template <class CombinedSpaceImp>
  class SubBaseFunctionSet;
  template <class CombinedSpaceImp>
  class SubMapper;

  template <class CombinedSpaceImp>
  struct SubSpaceTraits {
  private:
    typedef CombinedSpaceImp CombinedSpaceType;
    typedef typename CombinedSpaceType::Traits CombinedTraits;

    typedef typename CombinedTraits::RangeFieldType DofType;
    typedef typename CombinedTraits::RangeType CombinedRangeType;
    typedef typename CombinedTraits::JacobianRangeType CombinedJacobianRangeType;
    
    typedef typename CombinedTraits::BaseFunctionSetType CombinedBaseFunctionSetType;
    typedef typename CombinedTraits::MapperType CombinedMapperType;
    typedef typename CombinedTraits::ContainedMapperType ContainedMapperType;
    
    enum { CombinedDimRange = CombinedTraits::DimRange };

  public:
    typedef typename CombinedTraits::GridPartType GridPartType;
    
    // Assumption: only scalar contained function spaces
    enum { DimDomain = CombinedTraits::DimDomain,
           DimRange = 1 };

    typedef FunctionSpace<
      DofType, DofType, DimDomain, DimRange> FunctionSpaceType;
    typedef SubSpace<CombinedSpaceImp> DiscreteFunctionSpaceType;
    typedef SubBaseFunctionSet<CombinedSpaceImp> BaseFunctionSetImp; 
    typedef SimpleBaseFunctionProxy<BaseFunctionSetImp> BaseFunctionSetType;
    typedef SubMapper<CombinedSpaceImp> MapperType;

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename CombinedTraits::GridType GridType;
    typedef typename CombinedTraits::IndexSetType IndexSetType;
    typedef typename CombinedTraits::IteratorType IteratorType;
    typedef typename CombinedTraits::DofConversionType DofConversionType;

  public:
    //- Friends
    friend class SubSpace<CombinedSpaceImp>;
    friend class SubBaseFunctionSet<CombinedSpaceImp>;
    friend class SubMapper<CombinedSpaceImp>;
  };

  template <class CombinedSpaceImp>
  class SubSpace : 
    public DiscreteFunctionSpaceDefault<SubSpaceTraits<CombinedSpaceImp> >
  {
  public:
    //- Typedefs and enums
    typedef CombinedSpaceImp CombinedSpaceType;
    
    typedef SubSpace<CombinedSpaceType> ThisType;
    typedef SubSpaceTraits<CombinedSpaceType> Traits;
    typedef DiscreteFunctionSpaceDefault<Traits> BaseType;

    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::IteratorType IteratorType;
    typedef typename Traits::MapperType MapperType;
    typedef typename Traits::BaseFunctionSetImp  BaseFunctionSetImp;
    typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;

    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::DomainType DomainType;

  public:
    //- Public methods
    //! constructor
    SubSpace(CombinedSpaceType& spc, int component);

    //! destructor
    ~SubSpace();

    //! is data continuous?
    bool continuous() const { return spc_.continuous(); }

    enum { polynomialOrder = CombinedSpaceType::polynomialOrder};

    /** @copydoc DiscreteFunctionSpaceInterface::order */
    int order () const
    {
      return spc_.order();
    }

    //! access to gridPart
    GridPartType& gridPart() { return spc_.gridPart(); }
    //! access to gridPart
    const GridPartType& gridPart() const { return spc_.gridPart(); }

    //! begin iterator
    IteratorType begin() const { return spc_.begin(); }

    //! end iterator
    IteratorType end() const { return spc_.end(); }

    //! total number of dofs
    int size() const { return mapper_.size(); }

    //! map a local dof number to a global one
    template <class EntityType>
    int mapToGlobal(EntityType& en, int local) const {
      return mapper_.mapToGlobal(en, local);
    }

    //! access to baseFunctionSet for given Entity
    template <class EntityType>
    const BaseFunctionSetType baseFunctionSet(EntityType& en) const 
    {
      return baseFunctionSet(en.geometry().type()); 
    }

    //! access to baseFunctionSet for given Entity
    const BaseFunctionSetType baseFunctionSet(const GeometryType geomType) const 
    {
      assert(baseSetMap_.find( geomType ) != baseSetMap_.end());
      return BaseFunctionSetType(baseSetMap_[geomType]);
    }

    //! access to grid (const version)
    const GridType& grid() const { return spc_.grid(); }

    //! access to mapper
    MapperType& mapper() const { return mapper_; }

  private:
    //- Forbidden methods
    SubSpace(const ThisType& other);
    ThisType& operator=(const ThisType& other);

  private:
    //- Data members
    const CombinedSpaceType& spc_;
    mutable MapperType mapper_;
    int component_;

    typedef std::map< const GeometryType, BaseFunctionSetImp* > BaseFunctionMapType; 
    mutable BaseFunctionMapType baseSetMap_; 
  };



  // Idea: wrap contained base function set, since this is exactly what you 
  // need here (except for when you go for ranges of subfunctions...)
  template< class CombinedSpaceImp >
  struct SubBaseFunctionSetTraits {
    typedef SubSpaceTraits<CombinedSpaceImp> SpaceTraits;
    typedef typename SpaceTraits::FunctionSpaceType FunctionSpaceType;
    typedef SubBaseFunctionSet<CombinedSpaceImp> BaseFunctionSetType;
  };
  template< class CombinedSpaceImp >
  class SubBaseFunctionSet
  : public BaseFunctionSetDefault< SubBaseFunctionSetTraits< CombinedSpaceImp > >
  {
  public:
    //! type of the associated CombinedSpace
    typedef CombinedSpaceImp CombinedSpaceType;
    
    //! type of the traits
    typedef SubBaseFunctionSetTraits< CombinedSpaceType > Traits;
    typedef typename Traits::SpaceTraits SpaceTraits;
  private:
    typedef BaseFunctionSetDefault< Traits > BaseType;

  public:
    typedef typename SpaceTraits :: DiscreteFunctionSpaceType 
      DiscreteFunctionSpaceType;

    typedef typename SpaceTraits :: DomainType DomainType;
    typedef typename SpaceTraits :: RangeType RangeType;
    typedef typename SpaceTraits :: JacobianRangeType JacobianRangeType;

  private:
    typedef typename SpaceTraits :: CombinedRangeType CombinedRangeType;
    typedef typename SpaceTraits :: CombinedJacobianRangeType
      CombinedJacobianRangeType;
    enum { CombinedDimRange = SpaceTraits::CombinedDimRange };

    typedef typename SpaceTraits :: CombinedBaseFunctionSetType
      CombinedBaseFunctionSetType;

  protected:
    const CombinedBaseFunctionSetType bSet_;
    const int component_;
    
  public:
    // using BaseType :: evaluate;
    // using BaseType :: jacobian;

  public:
    SubBaseFunctionSet( const CombinedBaseFunctionSetType bSet,
                        int component)
    : bSet_( bSet ),
      component_( component )
    {
    }

    inline GeometryType geometryType () const
    {
      return bSet_.geometryType();
    }
    
    /** \copydoc Dune::BaseFunctionSetInterface::numBaseFunctions */
    inline int numBaseFunctions () const
    {
      return bSet_.numDifferentBaseFunctions();
      assert( bSet_.numBaseFunctions() % CombinedDimRange == 0 );
      return bSet_.numBaseFunctions() / CombinedDimRange;
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const
    {
      bSet_.evaluateScalar(baseFunction,diffVariable,x,phi);
      return;
      // Assumption: dimRange == 1
      CombinedRangeType tmp;
      bSet_.evaluate( baseFunction, diffVariable, x, tmp );
      phi[ 0 ] = tmp[ component_ ];
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const PointType &x,RangeType &phi) const */
    template< class PointType >
    void evaluate ( const int baseFunction,
                    const PointType &x,
                    RangeType &phi ) const
    {
      bSet_.evaluateScalar(baseFunction,x,phi);
      return;
      // Assumption: dimRange == 1
      CombinedRangeType tmp;
      bSet_.evaluate( baseFunction, x, tmp );
      phi[ 0 ] = tmp[ component_ ];
    }

    /** \copydoc Dune::BaseFunctionSetInterface::jacobian(const int baseFunction,const PointType &x,JacobianRangeType &phi) const */
    template< class PointType >
    void jacobian ( const int baseFunction,
                    const PointType &x,
                    JacobianRangeType &phi ) const
    {
      // Assumption: dimRange == 1
      CombinedJacobianRangeType tmp;
      bSet_.jacobian( baseFunction, x, tmp );
      phi[ 0 ] = tmp[ component_ ];
    }
  };



  template <class CombinedSpaceImp>
  class SubMapper : 
    public DofMapperDefault<SubMapper<CombinedSpaceImp> > 
  {
  public:
    //- Typedefs and enums
    typedef CombinedSpaceImp CombinedSpaceType;
    
    typedef SubMapper<CombinedSpaceType> ThisType;
    typedef SubSpaceTraits<CombinedSpaceType> Traits;

    typedef typename Traits::CombinedMapperType CombinedMapperType;
    typedef typename Traits::ContainedMapperType ContainedMapperType;
    typedef typename Traits::DofConversionType DofConversionType;

  public:
    //- Public methods
    SubMapper(const CombinedSpaceType& spc,
              int component) :
      spc_(spc),
      mapper_(spc.containedSpace().mapper()),
      component_(component),
      utilGlobal_(spc.containedSpace(),
                  spc.myPolicy() == PointBased ? 
                  spc.numComponents() :
                  spc.size()/spc.numComponents())
    {}
    SubMapper(const ThisType& other) :
      spc_(other.spc_),
      mapper_(other.mapper_),
      component_(other.component_),
      utilGlobal_(other.utilGlobal_) {}
    /*
    ThisType& operator=(const ThisType& other) {
      spc_ = other.spc_;
      mapper_ = other.mapper_;
      component_ = other.component_;
      utilGlobal_ = other.utilGlobal_;
    }
    */
    //! Total number of degrees of freedom
    int size() const {
      return mapper_.size();
    }
    int operator()(int containedGlobal) const {
      utilGlobal_.newSize(mapper_.size());
      return utilGlobal_.combinedDof(
          containedGlobal, component_);
    }

  private:
    //- Data members
    const CombinedSpaceType& spc_;
    const ContainedMapperType& mapper_;
    int component_;

    mutable DofConversionType utilGlobal_;
  };
} // end namespace Dune

#include "subspace.cc"

#endif
