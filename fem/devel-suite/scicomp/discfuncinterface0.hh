#include <dune/common/array.hh>
#include <dune/common/fvector.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/discretefunction/dfadapt.hh>
#include <dune/grid/io/visual/grapedatadisplay.hh>
#include <memory> 
struct LocalDiscFuncInterface;
// pure virtual interface 
struct DiscFuncInterface {
  typedef std::auto_ptr<LocalDiscFuncInterface> LocalDiscFuncInterfaceType;
  typedef GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  virtual double& operator[](int) = 0;
  virtual double operator[](int) const = 0;
  virtual int size() const = 0;
  virtual LocalDiscFuncInterfaceType
          operator[](const EntityType&) = 0; // why pointer?
  virtual double evaluate(const DomainType&) const {}
  virtual int map(const EntityType&,int) const = 0;
 protected:
  friend class LocalDiscFuncInterface;
};
// teils pure virtual interface 
// teils schon Implementierung/Datenspeicherung und Defaultimplementierung
// Alternative waere eine virtual getEntity(),getDF() methode einzubauen
// und die Defaultimplementierung damit durchzufuehren; fuehrt
// allerdings zu mehr Funktionsaufrufen. Nachteil dieser Technik ist, dass
// alle implementierungen diese Struktur erben und damit ein Festlegung 
// getroffen wird...
struct LocalDiscFuncInterface {
  typedef  GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  LocalDiscFuncInterface(DiscFuncInterface& df,const EntityType& en) :
    df_(df), en_(en) {}
  virtual double evaluate(const DomainType&) const = 0;
  virtual int size() = 0;
  virtual double& operator[](int i) { // possible default implementation
    return df_[df_.map(en_,i)];
  }
  virtual double operator[](int i) const { // possible default implementation
    return df_[df_.map(en_,i)];
  }
 protected:
  DiscFuncInterface& df_;
  const EntityType& en_;
};
// **************************************************************
// **************************************************************
class LocalContDiscFunc;
class ContDiscFunc : public DiscFuncInterface {
 public:
  typedef DiscFuncInterface::LocalDiscFuncInterfaceType 
  LocalDiscFuncInterfaceType;
  typedef  GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  typedef  GridType :: Traits :: LeafIndexSet IndexSetType;
  ContDiscFunc(const GridType& grid) : 
    grid_(grid),
    set_(grid.leafIndexSet()),
    dofs_(set_.size(dimDomain)),
    size_(set_.size(dimDomain))
  {
  }
  virtual double& operator[](int i) {
    return dofs_[i];
  }
  virtual double operator[](int i) const {
    return dofs_[i];
  }
  virtual int size() const {
    return size_;
    // return set_.size(dimDomain); // ganz schlecht!!!
  }
  virtual LocalDiscFuncInterfaceType
  operator[](const EntityType& en);      // warum nicht hier definieren?
  // only for testing : converts array to DiscreteFunction instance and calls 
  // grape
  void disp() {
    Dune::GrapeDataDisplay < GridType > disp (this->grid_);
    typedef Dune::LeafGridPart< GridType > GridPartType;
    typedef Dune::FunctionSpace < double, double , dimDomain , 1 > FuncSpace;
    typedef Dune::LagrangeDiscreteFunctionSpace < FuncSpace , GridPartType, 1 > FuncSpaceType ;
    typedef Dune::DFAdapt < FuncSpaceType > DiscFuncType;
    GridPartType gridpart( this->grid_ );
    FuncSpaceType space ( gridpart );
    DiscFuncType dfgrape ( "test", space , (double *) &dofs_[0] );
    disp.dataDisplay(dfgrape);
  }
 protected:
  virtual int map(const EntityType& en,int i) const {
    return set_.subIndex<3>(en,i);
  }
  const GridType& grid_;
  const IndexSetType& set_;
  Dune::Array<double> dofs_;
  int size_;
  friend class LocalContDiscFunc;
};
class LocalContDiscFunc : public LocalDiscFuncInterface {
 public:
  typedef  GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  LocalContDiscFunc(ContDiscFunc& df,const EntityType& en) :
    LocalDiscFuncInterface(df,en) {}
  virtual double evaluate(const DomainType& x) const {
    return 
      (*this)[0]*x[0] +
      (*this)[1]*x[1] +
      (*this)[2]*x[2];
  }
  virtual int size() {
    return dimDomain+1;
  }
};
ContDiscFunc :: LocalDiscFuncInterfaceType
ContDiscFunc :: operator[](const ContDiscFunc :: EntityType& en) {
  return LocalDiscFuncInterfaceType(new LocalContDiscFunc(*this,en)); 
}
typedef DiscFuncInterface IDiscFuncType;
