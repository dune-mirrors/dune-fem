namespace Dune
{
  
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
