#ifndef ELLIPTC_MODEL_HH
#define ELLIPTC_MODEL_HH

#include <cassert>
#include <cmath>

#include <algorithm>

template< class FunctionSpace, class GridPart >
struct Model
{
  typedef FunctionSpace FunctionSpaceType;
  typedef GridPart GridPartType;

  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPartType::IntersectionType IntersectionType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  static const int dimRange = FunctionSpaceType::dimRange;

  Model ( ) {}

  bool isDirichletIntersection ( const IntersectionType &intersection,
                                 Dune::FieldVector< int, dimRange > &dirichletComponent ) const
  {
    std::fill( dirichletComponent.begin(), dirichletComponent.end(), 0 );
    dirichletComponent[0] = 1;
    return (intersection.geometry().center()[1]>0.5);
  }
  template< class Point >
  RangeType dirichlet ( int bndId, const Point &x ) const
  {
    auto xgl = entity().geometry().global( x );
    return RangeType( 3.*xgl[0]*xgl[1] );
  }
  void bind ( const EntityType &entity ) { entity_ = entity; }
  void unbind () {}
  const EntityType &entity () const { return entity_; }

protected:
  EntityType entity_;
};

#endif // #ifndef ELLIPTC_MODEL_HH
