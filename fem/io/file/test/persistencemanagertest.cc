#include <iostream>
#include <cmath>
#include <config.h>
#include <dune/fem/io/file/persistencemanager.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

using namespace Dune;
using namespace std;
struct TestClass1 : public PersistentObject {
  TestClass1(double pa) : a(pa) {}
  void set(double pa) {a=pa;}
  virtual void backup() const {
    PersistenceManager::backupValue("class1",a);
  }
  virtual void restore() {
    PersistenceManager::restoreValue("class1",a);
  }
  double a;
};
struct TestClass2 : public PersistentObject {
  TestClass2(double pb) : b(pb) {}
  virtual void backup() const {
    ofstream out(PersistenceManager::uniqueFileName().c_str());
    out << b; 
  }
  virtual void restore() {
    ifstream in(PersistenceManager::uniqueFileName().c_str());
    in >> b; 
  }
  double b;
};
struct TestClassA : public AutoPersistentObject {
  TestClassA() : a(42.), b(42.) {}
  TestClassA(double pa,double pb) : a(pa), b(pb) {}
  void calc() {
    b=b*b;
    a=a+b;
  }
  virtual void backup() const {
    PersistenceManager::backupValue("classA",a);
    ofstream out(PersistenceManager::uniqueFileName().c_str());
    out << b; 
  }
  virtual void restore() {
    PersistenceManager::restoreValue("classA",a);
    ifstream in(PersistenceManager::uniqueFileName().c_str());
    in >> b; 
  }
  double a,b;
};
int main (int argc, char **argv) {
   Parameter::append(argc,argv);
   if(argc != 2) {
    fprintf(stderr,"usage: %s restart? (0/1) \n",argv[0]);
    exit(1);
   }
  int restart = atoi( argv[1] );
  
  TestClass1 test1a(1.0),test1b(-sqrt(2.));
  double a=-100,b=-200;
  std::string s="TEST";
  {
    TestClassA testA(42,42);
    int c=-42;
    TestClass1 test1c(0.42);
    persistenceManager << test1a << a << s << test1b << c << test1c;
    persistenceManager << c << test1a << a << s << test1b;
    persistenceManager >> c >> c >> test1c;
  }
  TestClass2 test2(1.0);
  TestClassA testA(5.0,2.0);
  TestClassA* testAptr = new TestClassA[10]; 
  persistenceManager << b << test2;

  typedef YaspGrid<2> GridType;
  GridPtr<GridType> gridptr("2dgrid.dgf");
  GridType& grid=*gridptr;
  GridTimeProvider<GridType> timeProv(grid);

  double param;
  persistenceManager << param;

  double* aPtr = new double;
  persistenceManager << *aPtr;
  
  if (restart) {
    PersistenceManager::restore("backup");
    Parameter::get("test",1e5,param);
  } else {
    a=10;
    b=20;
    s="HALLO";
    test1a.set(-10.0);
    testA.calc();
    testA.calc();
    Parameter::get("test",1e5,param);
    *aPtr = -17.;
    PersistenceManager::backup("backup");
  }
  {
    TestClassA testa(testA);
    delete [] testAptr;
  }
  PersistenceManager::backup("backup.end");
  cout << "WERTE: " << a << " " << b << " " << s << " ,  "
       << test1a.a << " " << test1b.a << " ,  " 
       << test2.b << " ,  " 
       << testA.a << " " << testA.b  << "  ,  "
       << param << "  ,  "
       << timeProv.time() << " "
       << *aPtr << " "
       << std::endl;
}
