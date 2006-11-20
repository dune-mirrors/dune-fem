#include <iostream>
#include "combineinterface.hh"

template <class Impl>
struct InterfaceA {
  int a(int i) {return asImp().a(i);}
  Impl& asImp() {return static_cast<Impl&>(*this);}
};
struct InterfaceAImpl1 : public InterfaceA<InterfaceAImpl1> {
  int a(int i) {return i*i;}
};
struct InterfaceAImpl2 : public InterfaceA<InterfaceAImpl2> {
  int a(int i) {return -i*i;}
};
// now how to combine two interface implementations 
template <class I1,class I2>
struct InterfaceAPair : 
  public InterfaceA<InterfaceAPair<I1,I2> >,
  public Dune::PairOfInterfaces<I1,I2> {
  // always use the following interface
  InterfaceAPair(I1 i1, I2 i2) : Dune::PairOfInterfaces<I1,I2>(i1,i2) {}
  // the result of the interface method should be summed up...
  int a(int i) {
    return this->first().a(i) + this->second().a(i);
  }
};
template <class A>
int test1(InterfaceA<A>& impl) {
  return impl.a(10);
}
template <class A1,class A2>
int test2(InterfaceA<A1>& impl1,
	  InterfaceA<A2>& impl2) {
  typedef InterfaceA<A1> I1Type;
  typedef InterfaceA<A2> I2Type;
  // Can now further combine the interfaces also using reference, pointer, or
  // simple types
  typedef Dune::CombineInterface<InterfaceAPair,I1Type*,I1Type&,I1Type&,
    InterfaceAImpl1&,InterfaceAImpl1&,InterfaceAImpl2&,InterfaceAImpl2&,
    I2Type&,I2Type> CIType;
  InterfaceAImpl1 a1;
  InterfaceAImpl2 a2;
  CIType ci(&impl1,impl1,impl1,a1,a1,a2,a2,impl2,impl2);
  return ci.a(10);
}
int main() {
  typedef Dune::CombineInterface<InterfaceAPair,
    InterfaceAImpl1&,InterfaceAImpl2&,InterfaceAImpl2&> InterfaceAImpl122;
  InterfaceAImpl1 i1;
  InterfaceAImpl2 i2,i3;
  InterfaceAImpl122 i122(i1,i2,i3);
  // InterfaceAImpl122 is again an interface implementation...
  typedef Dune::CombineInterface<InterfaceAPair,
    InterfaceAImpl1,InterfaceAImpl122,InterfaceAImpl1,InterfaceAImpl122> InterfaceAImpl11221122;
  InterfaceAImpl11221122 i11221122(i1,i122,i1,i122);
  // now use the interface method:
  std::cout << test1(i1) << " " << test1(i2) << " " << test1(i3) 
            << " " << test1(i122) << " " << test1(i11221122) << std::endl;
  std::cout << "Example of recombination: " 
	    << "result is A11111221122112211221122 (13xA1 + 10xA2)" << std::endl;
  std::cout << test2(i1,i11221122) << std::endl;
}
