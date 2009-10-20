#ifndef DUNE_FEM_MODELTRAITS_HH
#define DUNE_FEM_MODELTRAITS_HH

#include <dune/common/fvector.hh>

namespace Dune {

  // *** ModelTraits class
  template <class GridPart,int dimRange2,
      int dimRange1=dimRange2* GridPart::GridType::dimensionworld>
  class ModelTraits {
  public:
    typedef GridPart GridPartType;
    typedef typename GridPartType :: GridType GridType;
    enum { dimDomain = GridType::dimensionworld };
    enum { dimRange = dimRange2, dimGradRange = dimRange1 };
    typedef FieldVector<double, dimDomain> DomainType;
    typedef FieldVector<double, dimDomain-1> FaceDomainType;
    typedef FieldVector<double,dimRange> RangeType;
    typedef FieldVector<double,dimGradRange> GradientType;
    typedef FieldMatrix<double,dimRange,dimDomain> FluxRangeType;
    typedef FieldMatrix<double,dimGradRange,dimDomain> DiffusionRangeType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;
   };

} // end namespace Dune

#endif
