#ifndef DUNE_SIMPLEXPOINTS_HH
#define DUNE_SIMPLEXPOINTS_HH

#include <cassert>

//- Dune includes
#include <dune/common/fvector.hh>

// include pardg quadratures 
#include <dune/fem/solver/pardg.hh>

namespace Dune {

  //! Adapter to the quadratures defined by parDG.
  template <int dim>
  class ParDGSimplexPointsAdapter {
  public:
    enum { numCorners = dim+1 };
    typedef typename pardg::Quadrature<dim> ParDGQuadratureType;
    typedef FieldVector<double, dim> CoordinateType;

  public:
    //! Constructor.
    ParDGSimplexPointsAdapter(int order) :
      quad_(ParDGQuadratureType::quadrature(order)),
      order_(order)
    {
    }

    //! Number of quadrature points.
    int numPoints() const 
    {
      return quad_.number_of_points();
    }

    //! The actual order of the quadrature.
    int order() const 
    {
      return order_;
    }

    //! Access to the ith quadrature point.
    CoordinateType point(int i) const 
    {
      assert(i >= 0 && i < numPoints());
      CoordinateType result;
      for (size_t j = 0; j < dim; ++j) 
      {
        result[j] = quad_.x(i)[j];
      }
      return result;
    }

    //! Access to the ith quadrature weight.
    double weight(int i) const 
    {
      assert(i >= 0 && i < numPoints());
      // scale with volume of reference element!
      return quad_.w(i);
    }
    
  private:
    const ParDGQuadratureType& quad_;
    const int order_;
  };

} // end namespace Dune
#endif
