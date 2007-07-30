namespace Dune {

  //- class QuadratureImp
  template <typename ct, int dim>
  IntegrationPointListImp<ct, dim>::IntegrationPointListImp(size_t id) :
    points_(),
    id_(id)
  {}
    
  template <typename ct, int dim>
  void IntegrationPointListImp<ct, dim>::
  addIntegrationPoint(const CoordinateType& point)
  {
    points_.push_back(point);
  }

  //- class QuadratureImp
  template <typename ct, int dim>
  QuadratureImp<ct, dim>::QuadratureImp(size_t id) :
    BaseType(id),
    weights_()
  {}
    
  template <typename ct, int dim>
  void QuadratureImp<ct, dim>::
  addQuadraturePoint(const CoordinateType& point, ct weight)
  {
    // add integration point to list 
    this->addIntegrationPoint(point);
    // store weight 
    weights_.push_back(weight);
  }

// only if we use dune-fem quadratures  
#ifndef USE_DUNE_QUADRATURES
  template <class ct, int dim>
  SimplexQuadrature<ct, dim>::SimplexQuadrature(const GeometryType&, int order, size_t id) :
    QuadratureImp<ct, dim>(id),
    order_((order <= 0) ? 1 : order)
  {
#ifdef HAVE_ALBERTA_FOUND
    AlbertaSimplexPointsAdapter<dim> points(order);
#else
    UGSimplexPointsAdapter<dim> points(order);
#endif

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
    for (int i = 0; i < m; ++i) {
      CoordinateType local(0.0);

      local[0] = gp.point(m, 0);
      double weight = gp.weight(m, 0);
      
      this->addQuadraturePoint(local, weight);
    }
  }

  template <class ct>
  TriangleQuadrature<ct>::TriangleQuadrature(const GeometryType&, int order, size_t id) :
    QuadratureImp<ct, 2>(id),
    order_((order <= 0) ? 1 : order)
  {
#ifdef HAVE_ALBERTA_FOUND
    AlbertaSimplexPointsAdapter<2> points(order);
#else
    UGSimplexPointsAdapter<2> points(order);
#endif

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
#ifdef HAVE_ALBERTA_FOUND
    AlbertaSimplexPointsAdapter<3> points(order);
#else 
    UGSimplexPointsAdapter<3> points(order);
#endif

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
    const PrismPoints& points = PrismPoints::instance();

    int m = 0;
    for (int i = 0; i < PrismPoints::numQuads; ++i) {
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
#endif // end #infdef USE_DUNE_QUADRATURES
  
  template <class ct, int dim>
  TestQuadrature<ct, dim>::TestQuadrature(const GeometryType& geo, int order) :
    QuadratureImp<ct, dim>(IdProvider::instance().newId()),
    geo_(geo),
    order_(order)
  {}

  template <class ct, int dim>
  void TestQuadrature<ct, dim>::
  newQuadraturePoint(const CoordinateType& c, ct w) 
  {
    this->addQuadraturePoint(c, w);
  }

} // end namespace Dune
