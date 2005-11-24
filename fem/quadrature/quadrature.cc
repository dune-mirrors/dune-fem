#include "gausspoints.hh"

namespace Dune {

  //- class QuadratureImp
  template <typename ct, int dim>
  QuadratureImp<ct, dim>::QuadratureImp(size_t id) :
    points_(),
    weights_(),
    id_(id)
  {}
    
  template <typename ct, int dim>
  void QuadratureImp<ct, dim>::
  addQuadraturePoint(const CoordinateType& point, ct weight)
  {
    points_.push_back(point);
    weights_.push_back(weight);
  }

  template <class ct, int dim>
  CubeQuadrature<ct, dim>::CubeQuadrature(int order, size_t id) :
    QuadratureImp<ct, dim>(id),
    order_(order)
  {
    typedef FieldVector<ct, dim> CoordinateType;

    const GaussPoints& gp = GaussPoints::instance();
 
    // find the right Gauss Rule from given order
    int m = 0;
    for (int i = 0; i <= GaussPoints::MAXP; i++) {
      if (gp.order(i)>=order) {
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
  LineQuadrature<ct>::LineQuadrature(int order, size_t id) :
    QuadratureImp<ct, 1>(id),
    order_(order)
  {
    typedef FieldVector<ct, 1> CoordinateType;

    const GaussPoints& gp = GaussPoints::instance();
    
    int m=0;
    for (int i = 0; i <= GaussPoints::MAXP; i++) {
      if (gp.order(i)>=order) {
        m = i;
        break;
      }
    }
    if (m==0) DUNE_THROW(NotImplemented, "order not implemented");
    order_ = gp.order(m);
	  
    // fill in all the gauss points
    for (int i = 0; i < m; ++i) {
      CoordinateType local(0.0);

      local[0] = gp.point(m, 0);
      double weight = gp.weight(m, 0);
      
      this->addQuadraturePoint(local, weight);
    }
  }

  template <class ct>
  TriangleQuadrature<ct>::TriangleQuadrature(int order, size_t id) :
    QuadratureImp<ct, 2>(id),
    order_(order)
  {

  }

  template <class ct>
  QuadrilateralQuadrature<ct>::QuadrilateralQuadrature(int order, size_t id) :
    QuadratureImp<ct, 2>(id),
    order_(order)
  {
    typedef FieldVector<ct, 2> CoordinateType;

    const GaussPoints& gp = GaussPoints::instance();
    const int dim = 2;

    // find the right Gauss Rule from given order
    int m = 0;
    for (int i = 0; i <= GaussPoints::MAXP; i++) {
      if (gp.order(i)>=order) {
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
  TetrahedronQuadrature<ct>::TetrahedronQuadrature(int order, size_t id) :
    QuadratureImp<ct, 3>(id),
    order_(order)
  {

  }

  template <class ct>
  HexahedronQuadrature<ct>::HexahedronQuadrature(int order, size_t id) :
    QuadratureImp<ct, 3>(id),
    order_(order)
  {
    typedef FieldVector<ct, 3> CoordinateType;

    const GaussPoints& gp = GaussPoints::instance();
    const int dim = 3;

    // find the right Gauss Rule from given order
    int m = 0;
    for (int i = 0; i <= GaussPoints::MAXP; i++) {
      if (gp.order(i)>=order) {
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
  PrismQuadrature<ct>::PrismQuadrature(int order, size_t id) :
    QuadratureImp<ct, 3>(id),
    order_(order)
  {

  }

  template <class ct>
  PyramidQuadrature<ct>::PyramidQuadrature(int order, size_t id) :
    QuadratureImp<ct, 3>(id),
    order_(order)
  {

  }
} // end namespace Dune
