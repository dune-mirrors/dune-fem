#include <dune/common/array.hh>
#include <dune/common/fvector.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/discretefunction/dfadapt.hh>
#include <dune/grid/io/visual/grapedatadisplay.hh>
#include <stack> 
template <class DiscFuncImpl> struct LDFWrapper;
template <class DiscFuncImpl> struct DiscFuncTraits;
template <class DiscFuncImpl> struct LocalDiscFuncInterface;
// Using BN
template <class DiscFuncImpl>
struct DiscFuncInterface {
  typedef DiscFuncImpl ImplType;
  // hier typename
  typedef typename DiscFuncTraits<ImplType>::LocalDiscFuncType LocalFuncImpl;
  typedef typename DiscFuncTraits<DiscFuncImpl>::LocalDiscFuncInterfaceType 
          LocalDiscFuncInterfaceType;
  typedef GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
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
  friend class LDFWrapper<ImplType>;
  inline LocalDiscFuncInterface<ImplType>* 
         pop(const EntityType& en); 
  inline void push(LocalDiscFuncInterface<ImplType>* ldf);
  std::stack<LocalDiscFuncInterface<ImplType>*> ldfstack_;
  const ImplType& asImp() const {
    return static_cast<const ImplType&>(*this);
  }
  ImplType& asImp() {
    return static_cast<ImplType&>(*this);
  }
};
template <class DiscFuncImpl>
struct LocalDiscFuncInterface {
  typedef DiscFuncImpl DiscFuncType;
  typedef typename DiscFuncTraits<DiscFuncType>::LocalDiscFuncType ImplType;
  typedef  GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  LocalDiscFuncInterface(DiscFuncType& df,const EntityType& en) :
    df_(df), en_(&en) {}
  double evaluate(const DomainType& x) const {
    return asImp().evaluate(x);
  }
  int size() {
    return asImp().size();
  }
  double& operator[](int i) { // possible default implementation
    return df_[df_.map(*en_,i)];
  }
  double operator[](int i) const { // possible default implementation
    return df_[df_.map(*en_,i)];
  }
 protected:
  void set(const EntityType& en) {
    asImp().setImpl(en);
  }
  void setImpl(const EntityType& en) {
    en_ = &en;                      // Anweisung geht nicht fuer const Entity&
  }
  friend class LDFWrapper<DiscFuncType>;
  friend class DiscFuncInterface<DiscFuncType>;
  DiscFuncType& df_;
  const EntityType* en_;
 protected:
  const ImplType& asImp() const {
    return static_cast<const ImplType&>(*this);
  }
  ImplType& asImp() {
    return static_cast<ImplType&>(*this);
  }
};
template <class DiscFuncImpl>
struct LDFWrapper {
  LDFWrapper(LocalDiscFuncInterface<DiscFuncImpl>* ldf) : ldf_(ldf) {}
  ~LDFWrapper() {
    ldf_->df_.push(ldf_);
  }
  typedef  GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  double evaluate(const DomainType& x) const {
    return ldf_->evaluate(x);
  }
  int size() {
    return ldf_->size();
  }
  double& operator[](int i) { // possible default implementation
    return ldf_->operator[](i);
  }
  double operator[](int i) const { // possible default implementation
    return ldf_->operator[](i);
  }
 private:
  LocalDiscFuncInterface<DiscFuncImpl>* ldf_;
};
// **************************************************************
// **************************************************************
class ContDiscFunc;
class LocalContDiscFunc;
// durch template Speziallisierung die Implmementierungstypen
// bekannt geben....
template <> 
struct DiscFuncTraits<ContDiscFunc> {
  typedef ContDiscFunc DiscFuncType;
  typedef LocalContDiscFunc LocalDiscFuncType;
  typedef LDFWrapper<ContDiscFunc> LocalDiscFuncInterfaceType;
};
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
class LocalContDiscFunc : public LocalDiscFuncInterface<ContDiscFunc> {
  typedef LocalDiscFuncInterface<ContDiscFunc> BaseType;
 public:
  typedef  GridType::Codim<0>::Entity EntityType;
  enum { dimDomain = GridType::dimensionworld };
  typedef Dune::FieldVector<double, dimDomain> DomainType;
  LocalContDiscFunc(ContDiscFunc& df,const EntityType& en) :
    BaseType(df,en) {
      dofs_[0] = &(*this)[0];
      dofs_[1] = &(*this)[1];
      dofs_[2] = &(*this)[2];
    }
  void setImpl(const EntityType& en) {
    BaseType::setImpl(en);
    dofs_[0] = &(*this)[0];
    dofs_[1] = &(*this)[1];
    dofs_[2] = &(*this)[2];
  }
  double evaluate(const DomainType& x) const {
    return 
      (*(dofs_[0]))*x[0] +
      (*(dofs_[1]))*x[1] +
      (*(dofs_[2]))*x[2];
    //(*this)[0]*x[0] +
    //(*this)[1]*x[1] +
    //(*this)[2]*x[2];
  }
  int size() {
    return dimDomain+1;
  }
 private:
  double *dofs_[3];
};
template <class ImplType>
inline LocalDiscFuncInterface<ImplType>* 
DiscFuncInterface<ImplType>::pop(const EntityType& en) {
  LocalDiscFuncInterface<ImplType>* ret = 0;
  if (!ldfstack_.empty()) {
    ret = ldfstack_.top();
    ldfstack_.pop();
    ret->set(en);
  }
  return ret;
}
template <class ImplType>
inline void 
DiscFuncInterface<ImplType>::push(LocalDiscFuncInterface<ImplType>* ldf) {
  ldf->en_ = 0;   // Sicherheit
  ldfstack_.push(ldf);
}

ContDiscFunc :: LocalDiscFuncInterfaceType
ContDiscFunc :: operator[](const ContDiscFunc :: EntityType& en) {
  LocalDiscFuncInterface<ContDiscFunc>* ret=pop(en);       
  if (!ret)
    ret = new LocalContDiscFunc(*this,en); 
  return ret;   // achtung impl. Type-Umwandlung!
}

// So heist das Interface (statisch, will ich ein anderes muss ich hier was tun)
typedef DiscFuncInterface<ContDiscFunc> IDiscFuncType;
