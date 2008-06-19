#include <dune/common/array.hh>
#include <dune/common/fvector.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/discretefunction/dfadapt.hh>
#include <dune/grid/io/visual/grapedatadisplay.hh>
struct LocalDiscFuncInterface;
// pure virtual interface 
struct DiscFuncInterface {
  typedef LocalDiscFuncInterface LocalDiscFuncInterfaceType;
  typedef GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  virtual double& operator[](int) = 0;
  virtual double operator[](int) const = 0;
  virtual int size() const = 0;
  virtual LocalDiscFuncInterfaceType*
          operator[](const EntityType&) = 0; // why pointer?
  virtual double evaluate(const DomainType&) const {}
  virtual int map(const EntityType&,int) const = 0;
 protected:
  friend class LocalDiscFuncInterface;
};
// pure virtual interface 
struct LocalDiscFuncInterface {
  typedef  GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  virtual double evaluate(const DomainType&) const = 0;
  virtual int size() = 0;
  virtual double& operator[](int i) = 0;
  virtual double operator[](int i) const = 0;
};
/**************************************************************
 * Implementierung einer stueckweisen linearen stetigen 
 * diskreten Funktion
 **************************************************************/
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
  virtual LocalDiscFuncInterfaceType*
  operator[](const EntityType& en);      // warum nicht hier definieren?
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
    df_(df), en_(en) {}
  virtual double evaluate(const DomainType& x) const {
    return 
      (*this)[0]*x[0] +
      (*this)[1]*x[1] +
      (*this)[2]*x[2];
  }
  virtual int size() {
    return dimDomain+1;
  }
  virtual double& operator[](int i) { 
    return df_[df_.map(en_,i)];
  }
  virtual double operator[](int i) const { 
    return df_[df_.map(en_,i)];
  }
 private:
  ContDiscFunc& df_;
  const EntityType& en_;
};
ContDiscFunc :: LocalDiscFuncInterfaceType*
ContDiscFunc :: operator[](const ContDiscFunc :: EntityType& en) {
  return new LocalContDiscFunc(*this,en); 
}
 
// Typ des Interfaces
typedef DiscFuncInterface IDiscFuncType;
