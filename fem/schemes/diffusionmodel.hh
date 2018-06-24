#ifndef ELLIPTC_MODEL_HH
#define ELLIPTC_MODEL_HH

#include <cassert>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#define VirtualDiffusionModelMethods(POINT) \
  virtual void source ( const POINT &x,\
                const RangeType &value,\
                const JacobianRangeType &gradient,\
                RangeType &flux ) const = 0;\
  virtual void linSource ( const RangeType& uBar,\
                   const JacobianRangeType &gradientBar,\
                   const POINT &x,\
                   const RangeType &value,\
                   const JacobianRangeType &gradient,\
                   RangeType &flux ) const = 0;\
  virtual void diffusiveFlux ( const POINT &x,\
                       const RangeType &value,\
                       const JacobianRangeType &gradient,\
                       JacobianRangeType &flux ) const = 0;\
  virtual void linDiffusiveFlux ( const RangeType& uBar,\
                          const JacobianRangeType& gradientBar,\
                          const POINT &x,\
                          const RangeType &value,\
                          const JacobianRangeType &gradient,\
                          JacobianRangeType &flux ) const = 0;\
  virtual void fluxDivergence( const POINT &x,\
                         const RangeType &value,\
                         const JacobianRangeType &jacobian,\
                         const HessianRangeType &hessian,\
                         RangeType &flux) const = 0;\
  virtual void alpha(const POINT &x,\
             const RangeType &value,\
             RangeType &val) const = 0;\
  virtual void linAlpha(const RangeType &uBar,\
                const POINT &x,\
                const RangeType &value,\
                RangeType &val) const = 0;\
  virtual void dirichlet( int bndId, const POINT &x,\
                RangeType &value) const = 0;

#define WrapperDiffusionModelMethods(POINT) \
  virtual void source ( const POINT &x,\
                const RangeType &value,\
                const JacobianRangeType &gradient,\
                RangeType &flux ) const \
  { impl().source(x, value, gradient, flux); } \
  virtual void linSource ( const RangeType& uBar,\
                   const JacobianRangeType &gradientBar,\
                   const POINT &x,\
                   const RangeType &value,\
                   const JacobianRangeType &gradient,\
                   RangeType &flux ) const \
  { impl().linSource(uBar, gradientBar, x, value, gradient, flux); } \
  virtual void diffusiveFlux ( const POINT &x,\
                       const RangeType &value,\
                       const JacobianRangeType &gradient,\
                       JacobianRangeType &flux ) const \
  { impl().diffusiveFlux(x, value, gradient, flux); } \
  virtual void linDiffusiveFlux ( const RangeType& uBar,\
                          const JacobianRangeType& gradientBar,\
                          const POINT &x,\
                          const RangeType &value,\
                          const JacobianRangeType &gradient,\
                          JacobianRangeType &flux ) const \
  { impl().linDiffusiveFlux(uBar, gradientBar, x, value, gradient, flux); } \
  virtual void fluxDivergence( const POINT &x,\
                         const RangeType &value,\
                         const JacobianRangeType &jacobian,\
                         const HessianRangeType &hessian,\
                         RangeType &flux) const \
  { impl().fluxDivergence(x, value, jacobian, hessian, flux); } \
  virtual void alpha(const POINT &x,\
             const RangeType &value,\
             RangeType &val) const \
  { impl().alpha(x, value, val); } \
  virtual void linAlpha(const RangeType &uBar,\
                const POINT &x,\
                const RangeType &value,\
                RangeType &val) const \
  { impl().linAlpha(uBar, x, value, val); } \
  virtual void dirichlet( int bndId, const POINT &x,\
                RangeType &value) const \
  { impl().dirichlet(bndId, x, value); }

template< class GridPart, int dimR, class RangeField = double >
struct DiffusionModel
{
  typedef GridPart GridPartType;
  static const int dimRange = dimR;
  typedef DiffusionModel<GridPartType, dimRange, RangeField> ModelType;
  typedef RangeField RangeFieldType;

  typedef Dune::Fem::FunctionSpace< double, RangeFieldType,
              GridPart::dimensionworld, dimRange > FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::IntersectionType IntersectionType;
  typedef typename EntityType::Geometry::LocalCoordinate LocalDomainType;

  typedef typename Dune::Fem::CachingQuadrature< GridPartType, 0 >::
                   QuadraturePointWrapperType Point;
  typedef typename Dune::Fem::CachingQuadrature< GridPartType, 1 >::
                   QuadraturePointWrapperType IntersectionPoint;
  typedef typename Dune::Fem::ElementQuadrature< GridPartType, 0 >::
                   QuadraturePointWrapperType ElementPoint;
  typedef typename Dune::Fem::ElementQuadrature< GridPartType, 1 >::
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
struct DiffusionModelWrapper : public DiffusionModel<typename ModelImpl::GridPartType, ModelImpl::dimRange, typename ModelImpl::RangeFieldType>
{
  typedef ModelImpl Impl;
  typedef typename ModelImpl::GridPartType GridPartType;
  static const int dimRange  = ModelImpl::dimRange;
  typedef DiffusionModel<GridPartType, dimRange, typename ModelImpl::RangeFieldType> Base;

  typedef typename Base::Point Point;
  typedef typename Base::IntersectionPoint IntersectionPoint;
  typedef typename Base::ElementPoint ElementPoint;
  typedef typename Base::ElementIntersectionPoint ElementIntersectionPoint;
  typedef typename Base::LocalDomainType LocalDomainType;
  typedef typename Base::DomainType DomainType;
  typedef typename Base::RangeType RangeType;
  typedef typename Base::JacobianRangeType JacobianRangeType;
  typedef typename Base::HessianRangeType HessianRangeType;
  typedef typename Base::EntityType EntityType;
  typedef typename Base::IntersectionType IntersectionType;

  template< class... Args, std::enable_if_t< std::is_constructible< ModelImpl, Args &&... >::value, int > = 0 >
  explicit DiffusionModelWrapper ( Args &&... args )
    : impl_( std::forward< Args >( args )... )
  {}

  ~DiffusionModelWrapper()
  {
    std::cout << "In DiffusionModelWrapper destructor" << std::endl;
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
  virtual bool isDirichletIntersection( const IntersectionType& inter, Dune::FieldVector<int, dimRange> &dirichletComponent ) const
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
struct DiffusionModelEngine : public DiffusionModel<typename ModelTraits::GridPartType, ModelTraits::dimRange, typename ModelTraits::RangeFieldType>
{
  typedef typename ModelTraits::GridPartType GridPartType;
  static const int dimRange  = ModelTraits::dimRange;
  typedef DiffusionModel<GridPartType, dimRange> Base;

  typedef typename Base::Point Point;
  typedef typename Base::IntersectionPoint IntersectionPoint;
  typedef typename Base::ElementPoint ElementPoint;
  typedef typename Base::ElementIntersectionPoint ElementIntersectionPoint;
  typedef typename Base::LocalDomainType LocalDomainType;
  typedef typename Base::DomainType DomainType;
  typedef typename Base::RangeType RangeType;
  typedef typename Base::JacobianRangeType JacobianRangeType;
  typedef typename Base::HessianRangeType HessianRangeType;
  typedef typename Base::EntityType EntityType;
  typedef typename Base::IntersectionType IntersectionType;

  typedef typename ModelTraits::ModelType ModelImpl;

  DiffusionModelEngine(ModelImpl &model) : impl_(model) {}
  ~DiffusionModelEngine()
  {
    std::cout << "In DiffusionModelWrapper destructor" << std::endl;
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
  virtual bool isDirichletIntersection( const IntersectionType& inter, Dune::FieldVector<int, dimRange> &dirichletComponent ) const
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
