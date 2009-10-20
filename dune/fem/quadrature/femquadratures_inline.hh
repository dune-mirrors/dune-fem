#ifndef DUNE_FEM_FEMQUADRATURES_INLINE_HH
#define DUNE_FEM_FEMQUADRATURES_INLINE_HH

#include "femquadratures.hh"

namespace Dune
{

#define SimplexPointsAdapter ParDGSimplexPointsAdapter

// only if we use dune-fem quadratures  
  template <class ct, int dim>
  SimplexQuadrature<ct, dim>::SimplexQuadrature(const GeometryType&, int order, size_t id) :
    QuadratureImp<ct, dim>(id),
    order_((order <= 0) ? 1 : order)
  {
    SimplexPointsAdapter<dim> points(order);

    order_ = points.order();

    for (int i = 0; i < points.numPoints(); ++i) {
      this->addQuadraturePoint(points.point(i), points.weight(i));
    }
  }

  template <class ct, int dim>
  CubeQuadrature<ct, dim>::CubeQuadrature(const GeometryType&, int order, size_t id) :
    QuadratureImp<ct, dim>(id),
    order_((order <= 0) ? 1 : order)
  {
    typedef FieldVector<ct, dim> CoordinateType;

    const GaussPts& gp = GaussPts::instance();

    // find the right Gauss Rule from given order
    int m = 0;
    for (int i = 0; i <= GaussPts::MAXP; i++) {
      if (gp.order(i)>=order_) {
        m = i;
        break;
      }
    }

    if (m==0) DUNE_THROW(NotImplemented, "order not implemented");
    order_ = gp.order(m);

    // fill in all the gauss points
    int n = gp.power(m,dim);
    for (int i = 0; i < n; i++) {
      // compute multidimensional coordinates of Gauss point
      int x[dim];
      int z = i;
      for (int k=0; k<dim; ++k) {
        x[k] = z%m;
        z = z/m;
      }
      
      // compute coordinates and weight
      double weight = 1.0;
      FieldVector<ct, dim> local;
      for (int k = 0; k < dim; k++) {
        local[k] = gp.point(m,x[k]);
        weight *= gp.weight(m,x[k]);
      }

      // put in container
      this->addQuadraturePoint(local, weight);
    }

  }

  template <>
  inline CubeQuadrature<double, 0>::CubeQuadrature(const GeometryType&, int order, size_t id) :
    QuadratureImp<double, 0>(id),
    order_((order <= 0) ? 1 : order)
  {
	    typedef double ct;
    typedef FieldVector<ct, 0> CoordinateType;

    order_ = 20;

    // fill in all the gauss points
   // compute coordinates and weight
   double weight = 1.0;
   FieldVector<ct, 0> local(0);

   // put in container
   this->addQuadraturePoint(local, weight);

  }

  template <class ct>
  LineQuadrature<ct>::LineQuadrature(const GeometryType&, int order, size_t id) :
    QuadratureImp<ct, 1>(id),
    order_((order <= 0) ? 1 : order)
  {
    typedef FieldVector<ct, 1> CoordinateType;

    const GaussPts& gp = GaussPts::instance();
    
    int m=0;
    for (int i = 0; i <= GaussPts::MAXP; i++) {
      if (gp.order(i)>=order_) {
        m = i;
        break;
      }
    }
    if (m==0) DUNE_THROW(NotImplemented, "order not implemented");
    order_ = gp.order(m);
	  
    // fill in all the gauss points
		int n = gp.power(m,1);		
    for (int i = 0; i < n; ++i) {
      CoordinateType local(0.0);

      local[0] = gp.point(m, i);
      double weight = gp.weight(m, i);
      
      this->addQuadraturePoint(local, weight);
    }
  }

  template <class ct>
  TriangleQuadrature<ct>::TriangleQuadrature(const GeometryType&, int order, size_t id) :
    QuadratureImp<ct, 2>(id),
    order_((order <= 0) ? 1 : order)
  {
    SimplexPointsAdapter<2> points(order);

    order_ = points.order();

    for (int i = 0; i < points.numPoints(); ++i) {
      this->addQuadraturePoint(points.point(i), points.weight(i));
    }
  }

  template <class ct>
  QuadrilateralQuadrature<ct>::QuadrilateralQuadrature(const GeometryType&, int order, size_t id) :
    QuadratureImp<ct, 2>(id),
    order_((order <= 0) ? 1 : order)
  {
    typedef FieldVector<ct, 2> CoordinateType;

    const GaussPts& gp = GaussPts::instance();
    const int dim = 2;

    // find the right Gauss Rule from given order
    int m = 0;
    for (int i = 0; i <= GaussPts::MAXP; i++) {
      if (gp.order(i)>=order_) {
        m = i;
        break;
      }
    }

    if (m==0) DUNE_THROW(NotImplemented, "order not implemented");
    order_ = gp.order(m);

    // fill in all the gauss points
    int n = gp.power(m,dim);
    for (int i = 0; i < n; i++) {
      // compute multidimensional coordinates of Gauss point
      int x[dim];
      int z = i;
      for (int k=0; k<dim; ++k) {
        x[k] = z%m;
        z = z/m;
      }
      
      // compute coordinates and weight
      double weight = 1.0;
      FieldVector<ct, dim> local;
      for (int k = 0; k < dim; k++) {
        local[k] = gp.point(m,x[k]);
        weight *= gp.weight(m,x[k]);
      }

      // put in container
      this->addQuadraturePoint(local, weight);
    }
  }

  template <class ct>
  TetraQuadrature<ct>::TetraQuadrature(const GeometryType&, int order, size_t id) :
    QuadratureImp<ct, 3>(id),
    order_((order <= 0) ? 1 : order)
  {
    SimplexPointsAdapter<3> points(order);

    order_ = points.order();

    for (int i = 0; i < points.numPoints(); ++i) {
      this->addQuadraturePoint(points.point(i), points.weight(i));
    }
  }

  template <class ct>
  HexaQuadrature<ct>::HexaQuadrature(const GeometryType&, int order, size_t id) :
    QuadratureImp<ct, 3>(id),
    order_((order <= 0) ? 1 : order)
  {
    typedef FieldVector<ct, 3> CoordinateType;

    const GaussPts& gp = GaussPts::instance();
    const int dim = 3;

    // find the right Gauss Rule from given order
    int m = 0;
    for (int i = 0; i <= GaussPts::MAXP; i++) {
      if (gp.order(i)>=order_) {
        m = i;
        break;
      }
    }

    if (m==0) DUNE_THROW(NotImplemented, "order not implemented");
    order_ = gp.order(m);

    // fill in all the gauss points
    int n = gp.power(m,dim);
    for (int i = 0; i < n; i++) {
      // compute multidimensional coordinates of Gauss point
      int x[dim];
      int z = i;
      for (int k=0; k<dim; ++k) {
        x[k] = z%m;
        z = z/m;
      }
      
      // compute coordinates and weight
      double weight = 1.0;
      FieldVector<ct, dim> local;
      for (int k = 0; k < dim; k++) {
        local[k] = gp.point(m,x[k]);
        weight *= gp.weight(m,x[k]);
      }

      // put in container
      this->addQuadraturePoint(local, weight);
    }
  }

  template <class ct>
  PrismQuadrature<ct>::PrismQuadrature(const GeometryType&, int order, size_t id) :
    QuadratureImp<ct, 3>(id),
    order_((order <= 0) ? 1 : order)
  {
    SimplexPointsAdapter<2> simplexPoints(order);
    int simplexOrder = simplexPoints.order();

		const GaussPts& gp = GaussPts::instance();
    // find the right Gauss Rule from given order
    int m = 0;
    for (int i = 0; i <= GaussPts::MAXP; i++) {
      if (gp.order(i)>=order_) {
        m = i;
        break;
      }
    }
    if (m==0) DUNE_THROW(NotImplemented, "order not implemented");
    
		int gaussOrder = gp.order(m);
    int minOrder = ((simplexOrder < gaussOrder) ? simplexOrder : gaussOrder);
    order_ = minOrder;
  
    int numSimplexPoints = simplexPoints.numPoints();
    int numGaussPoints = gp.power(m,1);

		FieldVector<ct, 3> local;
		double weight, simplexWeight;
    
    for (int i = 0; i < numSimplexPoints; ++i) {
			local[0] = simplexPoints.point(i)[0];
			local[1] = simplexPoints.point(i)[1];
			simplexWeight = simplexPoints.weight(i);
      for (int j = 0; j < numGaussPoints; ++j) {
				local[2] = gp.point(m,j);
				weight = simplexWeight;
				weight *= gp.weight(m,j);
        this->addQuadraturePoint(local, weight);
      }    
		}
  }

  template <class ct>
  PyramidQuadrature<ct>::PyramidQuadrature(const GeometryType&, int order, size_t id) :
    QuadratureImp<ct, 3>(id),
    order_((order <= 0) ? 1 : order)
  {
    const PyramidPoints& points = PyramidPoints::instance();

    int m = 0;
    for (int i = 0; i < PyramidPoints::numQuads; ++i) {
      if (points.order(i) >= order_) {
        m = i;
        break;
      }
    }

    if (m==0) DUNE_THROW(NotImplemented, "order not implemented");
    order_ = points.order(m);

    // fill in the points
    for (int i = 0; i < points.numPoints(m); ++i) {
      this->addQuadraturePoint(points.point(m, i), points.weight(m, i));
    }
  }
  
} // end namespace Dune

#endif
