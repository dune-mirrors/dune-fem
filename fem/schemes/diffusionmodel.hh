#ifndef ELLIPTC_MODEL_HH
#define ELLIPTC_MODEL_HH

#include <cassert>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fempy/quadrature/fempyquadratures.hh>

#define VirtualDiffusionModelMethods(POINT) \
  virtual void source ( const POINT &x,\
                const DRangeType &value,\
                const DJacobianRangeType &gradient,\
                RRangeType &flux ) const = 0;\
  virtual void linSource ( const DRangeType& uBar,\
                   const DJacobianRangeType &gradientBar,\
                   const POINT &x,\
                   const DRangeType &value,\
                   const DJacobianRangeType &gradient,\
                   RRangeType &flux ) const = 0;\
  virtual void diffusiveFlux ( const POINT &x,\
                       const DRangeType &value,\
                       const DJacobianRangeType &gradient,\
                       RJacobianRangeType &flux ) const = 0;\
  virtual void linDiffusiveFlux ( const DRangeType& uBar,\
                          const DJacobianRangeType& gradientBar,\
                          const POINT &x,\
                          const DRangeType &value,\
                          const DJacobianRangeType &gradient,\
                          RJacobianRangeType &flux ) const = 0;\
  virtual void fluxDivergence( const POINT &x,\
                         const DRangeType &value,\
                         const DJacobianRangeType &jacobian,\
                         const DHessianRangeType &hessian,\
                         RRangeType &flux) const = 0;\
  virtual void alpha(const POINT &x,\
             const DRangeType &value,\
             RRangeType &val) const = 0;\
  virtual void linAlpha(const DRangeType &uBar,\
                const POINT &x,\
                const DRangeType &value,\
                RRangeType &val) const = 0;\
  virtual void dirichlet( int bndId, const POINT &x,\
                RRangeType &value) const = 0;

#define WrapperDiffusionModelMethods(POINT) \
  virtual void source ( const POINT &x,\
                const DRangeType &value,\
                const DJacobianRangeType &gradient,\
                RRangeType &flux ) const \
  { impl().source(x, value, gradient, flux); } \
  virtual void linSource ( const DRangeType& uBar,\
                   const DJacobianRangeType &gradientBar,\
                   const POINT &x,\
                   const DRangeType &value,\
                   const DJacobianRangeType &gradient,\
                   RRangeType &flux ) const \
  { impl().linSource(uBar, gradientBar, x, value, gradient, flux); } \
  virtual void diffusiveFlux ( const POINT &x,\
                       const DRangeType &value,\
                       const DJacobianRangeType &gradient,\
                       RJacobianRangeType &flux ) const \
  { impl().diffusiveFlux(x, value, gradient, flux); } \
  virtual void linDiffusiveFlux ( const DRangeType& uBar,\
                          const DJacobianRangeType& gradientBar,\
                          const POINT &x,\
                          const DRangeType &value,\
                          const DJacobianRangeType &gradient,\
                          RJacobianRangeType &flux ) const \
  { impl().linDiffusiveFlux(uBar, gradientBar, x, value, gradient, flux); } \
  virtual void fluxDivergence( const POINT &x,\
                         const DRangeType &value,\
                         const DJacobianRangeType &jacobian,\
                         const DHessianRangeType &hessian,\
                         RRangeType &flux) const \
  { impl().fluxDivergence(x, value, jacobian, hessian, flux); } \
  virtual void alpha(const POINT &x,\
             const DRangeType &value,\
             RRangeType &val) const \
  { impl().alpha(x, value, val); } \
  virtual void linAlpha(const DRangeType &uBar,\
                const POINT &x,\
                const DRangeType &value,\
                RRangeType &val) const \
  { impl().linAlpha(uBar, x, value, val); } \
  virtual void dirichlet( int bndId, const POINT &x,\
                RRangeType &value) const \
  { impl().dirichlet(bndId, x, value); }

template< class GridPart, int dimDomain, int dimRange=dimDomain, class RangeField = double >
struct DiffusionModel
{
  typedef GridPart GridPartType;
  static const int dimD = dimDomain;
  static const int dimR = dimRange;
  typedef DiffusionModel<GridPartType, dimD, dimR, RangeField> ModelType;
  typedef RangeField RangeFieldType;

  typedef Dune::Fem::FunctionSpace< double, RangeFieldType,
              GridPart::dimensionworld, dimD > DFunctionSpaceType;
  typedef Dune::Fem::FunctionSpace< double, RangeFieldType,
              GridPart::dimensionworld, dimR > RFunctionSpaceType;
  typedef typename DFunctionSpaceType::DomainType DomainType;
  typedef typename DFunctionSpaceType::RangeType DRangeType;
  typedef typename DFunctionSpaceType::JacobianRangeType DJacobianRangeType;
  typedef typename DFunctionSpaceType::HessianRangeType DHessianRangeType;
  typedef typename DFunctionSpaceType::DomainFieldType DDomainFieldType;
  typedef typename RFunctionSpaceType::RangeType RRangeType;
  typedef typename RFunctionSpaceType::JacobianRangeType RJacobianRangeType;
  typedef typename RFunctionSpaceType::HessianRangeType RHessianRangeType;
  typedef typename RFunctionSpaceType::DomainFieldType rDomainFieldType;

  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::IntersectionType IntersectionType;
  typedef typename EntityType::Geometry::LocalCoordinate LocalDomainType;

  template <class F,int d>
  using Traits = Dune::FemPy::FempyQuadratureTraits<F,d>;
  typedef typename Dune::Fem::CachingQuadrature< GridPartType, 0, Traits >::
                   QuadraturePointWrapperType Point;
  typedef typename Dune::Fem::CachingQuadrature< GridPartType, 1, Traits >::
                   QuadraturePointWrapperType IntersectionPoint;
  typedef typename Dune::Fem::ElementQuadrature< GridPartType, 0, Traits >::
                   QuadraturePointWrapperType ElementPoint;
  typedef typename Dune::Fem::ElementQuadrature< GridPartType, 1, Traits >::
                   QuadraturePointWrapperType ElementIntersectionPoint;

  /*
  static const bool isLinear;
  static const bool isSymmetric;
  */

public:
  DiffusionModel( )
  { }
  virtual ~DiffusionModel() {}

  virtual std::string name() const = 0;

  virtual bool init( const EntityType &entity) const = 0;

  VirtualDiffusionModelMethods(Point)
  VirtualDiffusionModelMethods(ElementPoint)
  VirtualDiffusionModelMethods(IntersectionPoint)
  VirtualDiffusionModelMethods(ElementIntersectionPoint)
  VirtualDiffusionModelMethods(LocalDomainType)

  virtual bool hasDirichletBoundary () const = 0;
  virtual bool hasNeumanBoundary () const = 0;
  virtual bool isDirichletIntersection( const IntersectionType& inter, Dune::FieldVector<int, dimRange> &dirichletComponent ) const = 0;
};

template < class ModelImpl >
struct DiffusionModelWrapper : public DiffusionModel<typename ModelImpl::GridPartType, ModelImpl::dimD, ModelImpl::dimR,
                                      typename ModelImpl::RRangeFieldType>
{
  typedef ModelImpl Impl;
  typedef typename ModelImpl::GridPartType GridPartType;
  static const int dimD  = ModelImpl::dimD;
  static const int dimR  = ModelImpl::dimR;
  typedef DiffusionModel<GridPartType, dimD, dimR, typename ModelImpl::RRangeFieldType> Base;

  typedef typename Base::Point Point;
  typedef typename Base::IntersectionPoint IntersectionPoint;
  typedef typename Base::ElementPoint ElementPoint;
  typedef typename Base::ElementIntersectionPoint ElementIntersectionPoint;
  typedef typename Base::LocalDomainType LocalDomainType;
  typedef typename Base::DomainType DomainType;
  typedef typename Base::DRangeType DRangeType;
  typedef typename Base::DJacobianRangeType DJacobianRangeType;
  typedef typename Base::DHessianRangeType DHessianRangeType;
  typedef typename Base::RRangeType RRangeType;
  typedef typename Base::RJacobianRangeType RJacobianRangeType;
  typedef typename Base::RHessianRangeType RHessianRangeType;
  typedef typename Base::EntityType EntityType;
  typedef typename Base::IntersectionType IntersectionType;

  template< class... Args, std::enable_if_t< std::is_constructible< ModelImpl, Args &&... >::value, int > = 0 >
  explicit DiffusionModelWrapper ( Args &&... args )
    : impl_( std::forward< Args >( args )... )
  {}

  ~DiffusionModelWrapper()
  {
  }

  WrapperDiffusionModelMethods(Point);
  WrapperDiffusionModelMethods(ElementPoint);
  WrapperDiffusionModelMethods(IntersectionPoint);
  WrapperDiffusionModelMethods(ElementIntersectionPoint);
  WrapperDiffusionModelMethods(LocalDomainType);

  // other virtual functions
  virtual std::string name() const
  {
    return impl().name();
  }
  virtual bool hasDirichletBoundary () const
  {
    return impl().hasDirichletBoundary();
  }
  virtual bool hasNeumanBoundary () const
  {
    return impl().hasNeumanBoundary();
  }
  virtual bool isDirichletIntersection( const IntersectionType& inter, Dune::FieldVector<int, dimR> &dirichletComponent ) const
  {
    return impl().isDirichletIntersection(inter, dirichletComponent);
  }
  virtual bool init( const EntityType &entity) const
  {
    return impl().init(entity);
  }
  const ModelImpl &impl() const
  {
    return impl_;
  }
  ModelImpl &impl()
  {
    return impl_;
  }
  private:
  ModelImpl impl_;
};

template < class ModelTraits >
struct DiffusionModelEngine : public DiffusionModel<typename ModelTraits::GridPartType, ModelTraits::dimD, ModelTraits::dimR, typename ModelTraits::RRangeFieldType>
{
  typedef typename ModelTraits::GridPartType GridPartType;
  static const int dimD  = ModelTraits::dimD;
  static const int dimR  = ModelTraits::dimR;
  typedef DiffusionModel<GridPartType, dimD, dimR> Base;

  typedef typename Base::Point Point;
  typedef typename Base::IntersectionPoint IntersectionPoint;
  typedef typename Base::ElementPoint ElementPoint;
  typedef typename Base::ElementIntersectionPoint ElementIntersectionPoint;
  typedef typename Base::LocalDomainType LocalDomainType;
  typedef typename Base::DomainType DomainType;
  typedef typename Base::DRangeType DRangeType;
  typedef typename Base::DJacobianRangeType DJacobianRangeType;
  typedef typename Base::DHessianRangeType DHessianRangeType;
  typedef typename Base::RRangeType RRangeType;
  typedef typename Base::RJacobianRangeType RJacobianRangeType;
  typedef typename Base::RHessianRangeType RHessianRangeType;
  typedef typename Base::EntityType EntityType;
  typedef typename Base::IntersectionType IntersectionType;

  typedef typename ModelTraits::ModelType ModelImpl;

  DiffusionModelEngine(ModelImpl &model) : impl_(model) {}
  ~DiffusionModelEngine()
  {
  }

  WrapperDiffusionModelMethods(Point);
  WrapperDiffusionModelMethods(ElementPoint);
  WrapperDiffusionModelMethods(IntersectionPoint);
  WrapperDiffusionModelMethods(ElementIntersectionPoint);
  WrapperDiffusionModelMethods(LocalDomainType);

  // other virtual functions
  virtual std::string name() const
  {
    return impl().name();
  }
  virtual bool hasDirichletBoundary () const
  {
    return impl().hasDirichletBoundary();
  }
  virtual bool hasNeumanBoundary () const
  {
    return impl().hasNeumanBoundary();
  }
  virtual bool isDirichletIntersection( const IntersectionType& inter, Dune::FieldVector<int, dimR> &dirichletComponent ) const
  {
    return impl().isDirichletIntersection(inter, dirichletComponent);
  }
  virtual bool init( const EntityType &entity) const
  {
    return impl().init(entity);
  }
  const ModelImpl &impl() const
  {
    return impl_;
  }
  ModelImpl &impl()
  {
    return impl_;
  }
  private:
  ModelImpl &impl_;
};


#endif // #ifndef ELLIPTC_MODEL_HH
