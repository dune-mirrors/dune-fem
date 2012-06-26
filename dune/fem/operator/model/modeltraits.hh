#ifndef DUNE_FEM_MODELTRAITS_HH
#define DUNE_FEM_MODELTRAITS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune 
{

  namespace Fem
  {

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
    typedef typename GridPartType::template Codim<0>::EntityType EntityType;
   };

  } // end namespace Fem

  // #if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: ModelTraits ;

  // #endif // DUNE_FEM_COMPATIBILITY
  

} // end namespace Dune

#endif
