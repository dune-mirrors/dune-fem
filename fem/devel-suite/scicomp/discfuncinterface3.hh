#include <dune/common/array.hh>
#include <dune/common/fvector.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/discretefunction/dfadapt.hh>
#include <dune/grid/io/visual/grapedatadisplay.hh>
template <class DiscFuncImpl> struct DiscFuncTraits;
template <class LocalDiscFuncImpl>
struct LocalDiscFuncInterface;
// Using BN
template <class DiscFuncImpl>
struct DiscFuncInterface {
  typedef GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  typedef DiscFuncImpl ImplType;
  // hier typename
  typedef typename DiscFuncTraits<ImplType>::LocalDiscFuncType LocalFuncImpl;
  typedef typename DiscFuncTraits<DiscFuncImpl>::LocalDiscFuncInterfaceType 
          LocalDiscFuncInterfaceType;
  double& operator[](int i) {
    return asImp()[i];
  }
  double operator[](int i) const {
    return asImp()[i];
  }
  int size() const {
    return asImp().size();
  }
  LocalDiscFuncInterfaceType operator[](const EntityType& en) {
    return asImp()[en];
  }
  double evaluate(const DomainType&) const {return 0.;}
  int map(const EntityType&,int) const {
    return asImp().map();
  }
 protected:
  const ImplType& asImp() const {
    return static_cast<const ImplType&>(*this);
  }
  ImplType& asImp() {
    return static_cast<ImplType&>(*this);
  }
};
// Using Engine
template <class DiscFuncImpl>
struct LocalDiscFuncInterface {
  typedef  GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  typedef DiscFuncImpl DiscFuncType;
  typedef typename DiscFuncTraits<DiscFuncImpl>::LocalDiscFuncType ImplType;
  LocalDiscFuncInterface(ImplType impl) :
    impl_(impl), df_(impl.df_), en_(impl.en_) {}
  LocalDiscFuncInterface(const LocalDiscFuncInterface& ldf) :
    impl_(ldf.impl_), df_(ldf.df_), en_(ldf.en_) {}
  double evaluate(const DomainType& x) {
    return impl_.evaluate(x);
  }
  int size() {
    return impl_.size();
  }
  double& operator[](int i) { // possible default implementation
    return df_[df_.map(en_,i)];
  }
  double operator[](int i) const { // possible default implementation
    return df_[df_.map(en_,i)];
  }
 private:
  ImplType impl_;
  DiscFuncType& df_;
  const EntityType& en_;
};
// **************************************************************
class ContDiscFunc;
class LocalContDiscFunc;
// durch template Speziallisierung die Implmementierungstypen
// bekannt geben....
template <> 
struct DiscFuncTraits<ContDiscFunc> {
  typedef ContDiscFunc DiscFuncType;
  typedef LocalContDiscFunc LocalDiscFuncType;
  typedef LocalDiscFuncInterface<DiscFuncType> LocalDiscFuncInterfaceType;
};
class LocalContDiscFunc;
class ContDiscFunc : public DiscFuncInterface<ContDiscFunc> {
 public:
  // hier kein typename...
  typedef DiscFuncTraits<ContDiscFunc>::LocalDiscFuncType LocalFuncImpl;
  typedef DiscFuncTraits<ContDiscFunc>::LocalDiscFuncInterfaceType 
  LocalDiscFuncInterfaceType;
  typedef GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  typedef  GridType :: Traits :: LeafIndexSet IndexSetType;
  ContDiscFunc(const GridType& grid) : 
    grid_(grid),
    set_(grid.leafIndexSet()),
    dofs_(set_.size(dimDomain)),
    size_(set_.size(dimDomain))
  {}
  double& operator[](int i) {
    return dofs_[i];
  }
  double operator[](int i) const {
    return dofs_[i];
  }
  int size() const {
    return size_;
    return set_.size(dimDomain);
  }
  LocalDiscFuncInterfaceType operator[](const EntityType& en);      
  // warum nicht hier definieren?
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
 public:
  int map(const EntityType& en,int i) const {
    return set_.subIndex<dimDomain>(en,i);
  }
  const GridType& grid_;
  const IndexSetType& set_;
  Dune::Array<double> dofs_;
  int size_;
  friend class LocalContDiscFunc;
};
class LocalContDiscFunc {
 public:
  typedef GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  LocalContDiscFunc(ContDiscFunc& df,const EntityType& en) :
    df_(df), en_(en) {
      // caching
      dofs_[0] = &df_[df_.map(en_,0)];
      dofs_[1] = &df_[df_.map(en_,1)];
      dofs_[2] = &df_[df_.map(en_,2)];
    }
  LocalContDiscFunc(const LocalContDiscFunc& ldf) :
    df_(ldf.df_), en_(ldf.en_) {
      dofs_[0] = ldf.dofs_[0];
      dofs_[1] = ldf.dofs_[1];
      dofs_[2] = ldf.dofs_[2];
    }
  double evaluate(const DomainType& x) {
    return 
      (*(dofs_[0]))*x[0] +
      (*(dofs_[1]))*x[1] +
      (*(dofs_[2]))*x[2];
    // (*this)[0]*x[0] +
    // (*this)[1]*x[1] +
    // (*this)[2]*x[2];
  }
  int size() {
    return dimDomain+1;
  }
 protected:
  friend class LocalDiscFuncInterface<ContDiscFunc>;
  ContDiscFunc& df_;
  const EntityType& en_;
 private:
  double *dofs_[3];
};
ContDiscFunc :: LocalDiscFuncInterfaceType
ContDiscFunc :: operator[](const ContDiscFunc :: EntityType& en) {
  // demonstrate auto_ptr here
  return LocalDiscFuncInterfaceType(LocalContDiscFunc(*this,en));
}
typedef DiscFuncInterface<ContDiscFunc> IDiscFuncType;
