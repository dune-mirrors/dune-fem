#include <iostream>
#include "combineinterface.hh"
using namespace Dune;
// Ein Interface mit einer void Methode
struct GInter {
  void hallo() {
    std::cout << "Grid Interface :-)" << std::endl;
  }
};
// ... nur zum testen eine zweite Implementierung ...
struct GIntertest {
  void hallo() {
    std::cout << "Gtest: Bye" << std::endl;
  }
};
// Ein "lokales" Interface, welches durch ein "globale" zugaenglich ist
// Implementierung 1
template <class Global>
struct LocalInter1 {
  LocalInter1(const LocalInter1&) {std::cout << "LocalInterface1-copy" << std::endl;}
  LocalInter1() {std::cout << "Local1()" << std::endl;}
  void test() {
    std::cout << "Local 1" << std::endl;
  }
};
// Implementierung 2
template <class Global>
struct LocalInter2 {
  LocalInter2(const LocalInter2&) {std::cout << "LocalInterface2-copy" << std::endl;}
  LocalInter2() {std::cout << "Local2()" << std::endl;}
  void test() {
    std::cout << "Local 2" << std::endl;
  }
};
// Und jetzt das globale Interface (2 Implementierungen)
template <class G>
struct GlobalInter1 {
  // das Template-Argument mit Referenz beim Anlegen und Zugriffsmethod
  typedef G GType;
  G &g_;
  GlobalInter1(G& g) : g_(g) {}
  G& g() {return g_;}
  // eine einfache interface methode
  int a(int i) const {
    // std::cerr << "In a1 " << i << std::endl;
    int k=4242; // nur test
    return 10*i;
  }
  enum {q=10};
  // Zugriff auf lokales Interface
  typedef LocalInter1<GlobalInter1<G> > LocalType;
  typedef LocalInter1<GlobalInter1<G> > RefLocalType; // dieser Typ ist nur
                                                      // zur Komb. noetig
  LocalType local() {
    return LocalType();
  }
  LocalType& reflocal() {
    return local_;
  }
  LocalType local_;
};
template <class G>
struct GlobalInter2 {
  typedef G GType;
  G &g_;
  GlobalInter2(G& g) : g_(g) {}
  G& g() {return g_;}
  int a(int i) const {
    // std::cerr << "In a2 " << i << std::endl;
    return 20*i;
  }
  enum {q=20};
  typedef LocalInter2<GlobalInter2<G> > LocalType;
  LocalType& reflocal() {
    return local_;
  }
  LocalType local_;
};
/*****************************************************************/
// Jetzt das neue
/*****************************************************************/
// Klasse zur Kombinieren zweier lokaler Interface-Klassen
template <class T1,class T2> 
struct CombLocalInter : public PairOfInterfaces<T1,T2> {
  // diese Zeilen muessen immer dabei sein
  CombLocalInter(T1 t1,T2 t2) : PairOfInterfaces<T1,T2>(t1,t2) {}
  // Interface methode
  void test() {
    this->first().test();
    this->second().test();
  }
};
// Kombination von "globalen" Interfaceklassen
template <class T1,class T2> 
struct CombGlobalInter : public PairOfInterfaces<T1,T2> {
 private:
  CombGlobalInter(const CombGlobalInter&) ;
  CombGlobalInter();
 public:
  // erstmal die Typen:
  typedef PairOfInterfaces<T1,T2> BaseType;
  typedef typename BaseType::T1Type T1Type;
  typedef typename BaseType::T2Type T2Type;
  typedef typename BaseType::T1Type::GType GType;// wobei GType fuer alle
                                                      // Interfaces gleich sein muss
  typedef CombLocalInter<typename T1Type::LocalType&,
                         typename T2Type::LocalType&> LocalType;
  // ob der GType und die Referenz in T1 und T2 gleich sind kann leicht
  // getestet werden ...
  typedef typename T2Type::GType G2Type;  
  // ... aber erst im Konstruktor 
  CombGlobalInter(T1 t1,T2 t2) : PairOfInterfaces<T1,T2>(t1,t2),
    local_(t1.reflocal(),t2.reflocal()) {
    // same type?
    IsTrue<is_same<GType,G2Type>::value>::yes();
    // same reference?
    // if (&(t1.g())!=&(t2.g())) abort();
  }
  // Nun die Interface-Methoden
  G2Type& g() {return this->first().g();}
  int a(int i) {
    // std::cerr << "In a<T1,T2> " << i << std::endl;
    return this->first().a(i) + this->second().a(i);
  }
  enum {q=T1Type::q+T2Type::q};
  LocalType& reflocal() {
    // Hier werden die lokalen Objekte sehr haeufig Kopiert (>130), da die
    // Methode local() keine Referenzen zurueckgeben, bzw.
    // die LocalType kein Referenztype ist
    return local_;
  }
  LocalType local_;
};
/*****************************************************/
int main() {
  GInter g;
  GInter g1;
  typedef GlobalInter1<GInter> II1;
  typedef GlobalInter2<GInter> II2;
  II1 o1(g);
  II2 o2(g),o3(g);
  // Combiniere II1 und II2 ueber 
  typedef CombineInterface<CombGlobalInter,II1&,II2&> co12Type;
  co12Type co12(o1,o2);
  // Combiniere II12 mit II2 (ginge auch mit
  //                          CombineInterface<CombGlobalInter,co12Type,II2&>
  typedef CombineInterface<CombGlobalInter,co12Type&,II2&> co122Type;
  co122Type co122(co12,o3);
  // Jetzt legen wir richtig los
  typedef CombineInterface<CombGlobalInter,co12Type&,II2&,co12Type&> co12212Type;
  co12212Type co12212(co12,o3,co12);
  // und noch mehr (im Moment Maximum = 4)
  typedef CombineInterface<CombGlobalInter,II2&,co12212Type&,co12212Type&,II2&>
    co212212122122Type;
  // co212212122122 ist mir zu lang...
  typedef co212212122122Type coType;
  coType co(o3,co12212,co12212,o2);
  // und jetzt benutzen
  std::cout << "Combined Type is: 2,1,2,2,1,2,1,2,2,1,2,2" 
            << std::endl
            << "  a(-29) = (20+10+20+20+10+20+10+20+20+10+20+20)*(-29) = " << -5080
	    << std::endl
	    << "  a(-29) = " << co.a(-29)
	    << std::endl
	    << "  q = 20+10+20+20+10+20+10+20+20+10+20+20 = " << 200
	    << std::endl
            << "  q = " << coType::q 
	    << std::endl; 
  // Gemeinsamme Typen und Referenzen koennen durchgereicht werden
  coType::GType& grid=co.g();
  grid.hallo();
  // oder durch Kombination von Referenzobjekten
  coType::LocalType& rlocal=co.reflocal();
  rlocal.test();
}
