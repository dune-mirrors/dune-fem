#include <iostream>
#include <cmath>
#include <config.h>
#include <dune/fem/io/file/persistencemanager.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

using namespace Dune;
using namespace Fem;
using namespace std;

/* A test class derived from PersistentObject
 * storing one double value. This value is written to
 * the global checkpoint file.
 */
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
/* A test class derived from PersistentObject
 * storing one double value. This value is written to
 * a seperate file - the filename is obtained from the PersistenceManager.
 */
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
/* This is a class derived from AutopersistentObject.
 * It is automatically added to the PersistenceManager on
 * construction.
 */
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

/* The main test program */
int main (int argc, char **argv) {

  MPIManager :: initialize( argc, argv );
  Parameter::append(argc,argv);
  int restart = (argc>2)? atoi( argv[1] ) : 0;

  // test adding and removing objects
  TestClass1 test1a(1.0),test1b(-sqrt(2.));  // derived from PeristentObject
  double a=-100,b=-200;                      // build in types can also be handeled
  std::string s="TEST";
  {
    TestClassA testA(42,42);   // this object is automatically persistent
    int c=-42;
    TestClass1 test1c(0.42);
    // first add objects
    persistenceManager << test1a << a << s << test1b << c << test1c;
    persistenceManager << c << test1a << a << s << test1b;
    // now remove a few
    persistenceManager >> c >> c >> test1c;
    // note that objects can be added/removed more than once - they are nevertheless
    // stored/restored only once
  }   // note that testA is removed here - but this is allowed before the first call to backup/restore

  // some more objects
  TestClass2 test2(1.0);
  TestClassA testA(5.0,2.0);
  TestClassA* testAptr = new TestClassA[10];
  persistenceManager << b << test2;

  // the TimeProvider is an example of an AutoPersistent object
  typedef YaspGrid<2> GridType;
  GridPtr<GridType> gridptr("2dgrid.dgf");
  GridType& grid=*gridptr;
  GridTimeProvider<GridType> timeProv(grid);
  GridTimeProvider<GridType> timeProv1(1,grid);

  double param;
  persistenceManager << param;

  double* aPtr = new double;
  persistenceManager << *aPtr;

  if (restart) {
    // restore
    PersistenceManager::restore("backup");
    Parameter::get("test",1e5,param);
  } else {
    // backup
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
    std::cout << "Try to construct an auto-persistent object after the first call to backup:"
              << std::endl;
    TestClassA testa(testA);
    std::cout << "Delete an auto-persistent object between calls to backup:"
              << std::endl;
    delete [] testAptr;
  }
  PersistenceManager::backup("backup.end");

  // write all objects
  cout << "Data:     "
       << a << " " << b << " " << s << " ,  "
       << test1a.a << " " << test1b.a << " ,  "
       << test2.b << " ,  "
       << testA.a << " " << testA.b  << "  ,  "
       << param << "  ,  "
       << timeProv.time() << " "
       << timeProv1.time() << "  ,  "
       << *aPtr << " "
       << std::endl;
  cout << "Expected: "
       << 10 << " " << 20 << " " << "HALLO" << " ,  "
       << -10 << " " << -sqrt(2) << " ,  "
       << 1 << " ,  "
       << 25 << " " << 16  << "  ,  "
       << 100000 << "  ,  "
       << 0 << " "
       << 1 << "  ,  "
       << -17 << " "
       << std::endl;

  return 0;
}
