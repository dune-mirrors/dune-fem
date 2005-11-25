#ifndef DUNE_SIMPLEXPOINTS_HH
#define DUNE_SIMPLEXPOINTS_HH

#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/quadrature/fixedorder/ugquadratures.hh>

namespace Dune {

  template <int dim>
  class UGSimplexPointsAdapter {
  public:
    enum { numCorners = dim+1 };
    typedef UG_Quadratures::QUADRATURE UGQuadratureType;
    typedef FieldVector<double, dim> CoordinateType;

  public:
    UGSimplexPointsAdapter(int order) :
      quad_(0),
      refVol_(0.0)
    {
      quad_ = UG_Quadratures::GetQuadratureRule(dim, numCorners, order);

      switch (dim) {
      case 1:
        refVol_ = 1.0;
      case 2:
        refVol_ = 0.5;
      case 3:
        refVol_ = 1.0/3.0;
      default:
        DUNE_THROW(NotImplemented,
                   "Quadrature points for simplices only up to dim == 3");
      }
    }

    int numPoints() const 
    {
      return quad_->nip;
    }

    int order() const 
    {
      return quad_->order;
    }

    CoordinateType point(int i) const 
    {
      assert(i >= 0 && i < numPoints());
      CoordinateType result;
      for (size_t j = 0; j < dim; ++j) {
        result[j] = quad_->local[i][j];
      }
      return result;
    }

    double weight(int i) const 
    {
      assert(i >= 0 && i < numPoints());
      // scale with volume of reference element!
      return refVol_*quad_->weight[i];
    }
    
  private:
    UGQuadratureType* quad_;
    double refVol_;
  };

  /* Not used: UGSimplexPointsAdapter does the same
  class UGTrianglePointsAdapter {
  public:
    enum { dim = 2 };
    enum { numCorners = 3 };
    typedef UG_Quadratures::QUADRATURE UGQuadratureType;
    typedef FieldVector<double, dim> CoordinateType;

  public:
    UGTrianglePointsAdapter(int order) :
      quad_(0) 
    {
      quad_ = UG_Quadratures::GetQuadratureRule(dim, numCorners, order);
    }

    int numPoints() const 
    {
      return quad_->nip;
    }

    int order() const 
    {
      return quad_->order;
    }

    CoordinateType point(int i) const 
    {
      assert(i >= 0; i < numPoints());
      CoordinateType result;
      result[0] = quad_->local[i][0];
      result[1] = quad_->local[i][1];
      return result;
    }

    double weight(int i) const 
    {
      assert(i >= 0; i < numPoints());
      // scale with volume of reference element!
      return 0.5*quad_->weight[i];
    }
    
  private:
    UGQuadratureType* ugQuad_;
  };

  class UGTetrahedronPointsAdapter {
  public:
    enum { dim = 2 };
    enum { numCorners = 3 };
    typedef UG_Quadratures::QUADRATURE UGQuadratureType;
    typedef FieldVector<double, dim> CoordinateType;

  public:
    UGTetrahedronPointsAdapter(int order) :
      quad_(0) 
    {
      quad_ = UG_Quadratures::GetQuadratureRule(dim, numCorners, order);
    }

    int numPoints() const 
    {
      return quad_->nip;
    }

    int order() const 
    {
      return quad_->order;
    }

    CoordinateType point(int i) const 
    {
      assert(i >= 0; i < numPoints());
      CoordinateType result;
      result[0] = quad_->local[i][0];
      result[1] = quad_->local[i][1];
      result[2] = quad_->local[i][2];
      return result;
    }

    double weight(int i) const 
    {
      assert(i >= 0; i < numPoints());
      // scale with volume of reference element!
      return 0.166666666666666666666667*quad_->weight[i];
    }
  private:    
    UGQuadratureType* ugQuad_;
  };
  */

} // end namespace Dune

#endif
