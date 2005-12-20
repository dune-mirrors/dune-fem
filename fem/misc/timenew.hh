#ifndef NEW_TIMESTEPPING_HH
#define NEW_TIMESTEPPING_HH

#include <memory>

#include "inverseoperatorfactory.hh"
#include "timeutility.hh"

#include <strstream>
//#include "../visual/dx/dxdata.hh"

namespace Dune {

// Solve D_1 u' + D_0 u = f
// Assume D_1 and D_0 to be linear operators, f independent from u
template <class DiscreteFunctionImp>
class Loop : public TimeProvider {
public:
  // Typedef
  typedef DiscreteFunctionImp DiscreteFunctionType;

public:  
  // Constructor
  Loop(DiscreteFunctionImp& sol,
       double dt = 0.01,
       double startTime = 0.0,
       int startLevel = 0) :
    TimeProvider(startTime),
    sol_(sol), 
    dt_(dt),
    level_(startLevel) { 
    std::cout << "Loop::Loop: dt = " << dt_ << std::endl;
  }
  
  // Routines
  void solve();
  void solveUntil(int endLevel);
  void solveUntil(double endTime);
  
  // Accessor functions
  int level() const { return level_; }
  DiscreteFunctionType& solution() { return sol_; }
  
protected:
  virtual void next() = 0;
  double dt() { return dt_; }
  
  DiscreteFunctionType& sol_;

private: 
  double dt_;
  int level_;
};

template <class DiscreteFunctionImp>
void Loop<DiscreteFunctionImp>::solve() {
  // Preparations 
  this->resetTimeStepEstimate();
  
  // Call next()
  next();
  
  // Update time and level
  this->augmentTime(dt_);
  ++level_;

  // Output info
  std::cout << this->time() << " time step finished! (dt = " << dt_ << ")" << std::endl;
}

template <class DiscreteFunctionImp>
void Loop<DiscreteFunctionImp>::solveUntil(int endLevel) {
  assert(endLevel > level_);
  while (level_ <= endLevel) {
    solve();
  }
}

template <class DiscreteFunctionImp>
void Loop<DiscreteFunctionImp>::solveUntil(double endTime) {
  assert(endTime > this->time());
  while (time() <= endTime) {
    solve();
  }
}


template <class DiscreteFunctionImp, class MappingImp>
class ExplicitEuler : public Loop<DiscreteFunctionImp>
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
  typedef MappingImp MappingType;
public:
  ExplicitEuler(DiscreteFunctionImp& sol,
                MappingType& timeOp,
                MappingType& spaceOp,
                InverseOperatorFactory<DiscreteFunctionType>& factory,
                double dt,
                double starttime = 0.0, 
                int startlevel = 0) :
    Loop<DiscreteFunctionImp>(sol, dt, starttime, startlevel),
    timeOp_(timeOp),
    spaceOp_(spaceOp),
    factory_(factory)
  {
    std::cout << "ExplicitEuler::ExplicitEuler: dt = " << dt << std::endl;
    // Set up inverse operator
    std::auto_ptr<MappingType> tmp(factory.createOperator(timeOp));
    inv_ = tmp;
  }

  virtual ~ExplicitEuler() {}

  // do one timestep 
  void next()
  {
    std::ostringstream filestream;
    filestream << "upd" << this->level();

    DiscreteFunctionType* temp = 
      new DiscreteFunctionType("euler", this->sol_.getFunctionSpace());
    DiscreteFunctionType* upd = 
      new DiscreteFunctionType("upd", this->sol_.getFunctionSpace());

    //DXWriter<typename DiscreteFunctionImp::DiscreteFunctionSpaceType, false> dx(upd->getFunctionSpace(), filestream.str());
  


    //DXWriter<FuncSpaceType,false> dx(pd_.space(), stream.str());

    //dx.write(sol_, "bef");

    spaceOp_(this->sol_, *upd);

    //dx.write(*upd, "data");

    (*(this->inv_))(*upd, *temp);

    this->sol_.add(*temp, this->dt());

    //dx.write(sol_, "aft");
 
    delete temp;
    delete upd;
  }

private:
  ExplicitEuler();
  ExplicitEuler(const ExplicitEuler&);
  ExplicitEuler& operator= (const ExplicitEuler&);

private:
  MappingType& timeOp_;
  MappingType& spaceOp_;

  InverseOperatorFactory<DiscreteFunctionType>& factory_;
  std::auto_ptr<MappingType> inv_;
};

} // end namespace Dune
#endif
