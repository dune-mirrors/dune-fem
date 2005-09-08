#ifndef __TIMESTEPPING_HH__
#define __TIMESTEPPING_HH__

#include <memory>
#include "description.hh"
#include "inverseoperatorfactory.hh"

#include <strstream>
//#include "../../dune-dxwriter/dxdata.hh"

using namespace Dune;

// Solve D_1 u' + D_0 u = f
// Assume D_1 and D_0 to be linear operators, f independent from u
template <class ProblemDescriptionType>
class Loop {
public:
  // Typedef
  typedef typename ProblemDescriptionType::DiscreteFunctionType DiscreteFunctionType;
  
  // Constructor
  Loop(ProblemDescriptionType& pd,
       double dt = 0.01,
       double startTime = 0.0,
       int startLevel = 0) :
    pd_(pd),
    sol_(pd.initialData()),
    dt_(dt),
    time_(startTime),
    level_(startLevel) { 
    std::cout << "Loop::Loop: dt = " << dt_ << std::endl;
  }
  
  // Routines
  void solve();
  void solveUntil(int endLevel);
  void solveUntil(double endTime);
  
  // Accessor functions
  double time() const { return time_; }
  int level() const { return level_; }
  DiscreteFunctionType&  solution() { return sol_; }
  
protected:
  virtual void next() = 0;
  double dt() { return dt_; }
  
  DiscreteFunctionType& sol_;
  ProblemDescriptionType& pd_;
private:
 
  double dt_;
  double time_;
  int level_;
};

template <class ProblemDescriptionType>
void Loop<ProblemDescriptionType>::solve() {
  // Preparations (nothing to do right now)
  
  // Call next()
  next();
  
  // Update time and level
  time_ += dt_;
  ++level_;

  // Output info
  std::cout << time_ << " time step finished! (dt = " << dt_ << ")" << std::endl;
}

template <class ProblemDescriptionType>
void Loop<ProblemDescriptionType>::solveUntil(int endLevel) {
  assert(endLevel > level_);
  while (level_ <= endLevel) {
    solve();
  }
}

template <class ProblemDescriptionType>
void Loop<ProblemDescriptionType>::solveUntil(double endTime) {
  assert(endTime > time_);
  while (time_ <= endTime) {
    solve();
  }
}


template <class ProblemDescriptionType>
class ExplicitEuler : public Loop<ProblemDescriptionType>
{
  typedef typename ProblemDescriptionType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
  typedef typename ProblemDescriptionType::MappingType MappingType;
public:
  ExplicitEuler(ProblemDescriptionType& pd,
                InverseOperatorFactory<DiscreteFunctionType>& factory,
                double dt,
                double starttime = 0.0, 
                int startlevel = 0) :
    Loop<ProblemDescriptionType>(pd, dt, starttime, startlevel),
    factory_(factory)
  {
    std::cout << "ExplicitEuler::ExplicitEuler: dt = " << dt << std::endl;
    // Set up inverse operator
    MappingType& op = pd.D1();
    std::auto_ptr<MappingType> tmp(factory.createOperator(op));
    inv_ = tmp;
  }

  virtual ~ExplicitEuler() {}

  // do one timestep 
  void next()
  {
    int level = this->pd_.level();

    DiscreteFunctionType * rhs = &this->pd_.rhs(this->time());
    DiscreteFunctionType * temp = this->pd_.temporary();

    //std::stringstream stream;
    //stream << "temp" << this->level();
    //DXWriter<FuncSpaceType,false> dx(pd_.space(), stream.str());

    //dx.write(sol_, "bef");

    this->pd_.D0()(this->sol_, *temp);

    //dx.write(*temp, "upd");

    rhs->add(*temp, -1.0);

    (*(this->inv_))(*rhs, *temp);

    //dx.write(*temp, "inv");

    this->sol_.add(*temp, this->dt());

    //dx.write(sol_, "aft");
 
    delete temp;
  }

private:
  InverseOperatorFactory<DiscreteFunctionType>& factory_;
  std::auto_ptr<MappingType> inv_;
};

template <class ProblemDescriptionType>
class ImplicitEuler : public Loop<ProblemDescriptionType>
{
  typedef typename ProblemDescriptionType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
  typedef typename ProblemDescriptionType::MappingType MappingType;
public:
  ImplicitEuler(ProblemDescriptionType& pd,
                InverseOperatorFactory<DiscreteFunctionType>& factory,
                double dt,
                double starttime = 0.0, 
                int startlevel = 0) :
    Loop<ProblemDescriptionType>(pd, dt, starttime, startlevel),
    factory_(factory)
  {
    // Set up inverse operator
    MappingType& op = pd.D0();
    op *= dt;
    op += pd.D1();

    std::auto_ptr<MappingType> inv(factory.createOperator(op));
    inv_ = inv;
  }

  // do one timestep 
  void next()
  {
    // * Probably not used
    //int level = pd_.level();

    DiscreteFunctionType * rhs = &this->pd_.rhs(this->time()+this->dt());
    DiscreteFunctionType * temp = this->pd_.temporary();

    // Calculate D_1 u
    this->pd_.D1()(this->sol_, *temp);
   
    // Calculate inv(D_1 + dt D_0)(D_1 u)
    (*inv_)(*temp, this->sol_);

    // Add rhs contribution
    (*inv_)(*rhs, *temp);
    this->sol_.add(*temp, this->dt());
    
    // Temporary output
    /*    
    Descr::LaplaceOperatorType* lap = 
      dynamic_cast<Descr::LaplaceOperatorType*>(&pd_.D0(pd_.level()));
    Descr::MassMatrixOpType* mop = 
      dynamic_cast<Descr::MassMatrixOpType*>(&pd_.D1(pd_.level()));

    const SparseRowMatrix<double>* mat = lap->getMatrix();
    std::cout << "Laplace Operator" << std::endl;
    mat->print(std::cout);

    const SparseRowMatrix<double>* mat2 = mop->getMatrix();
    std::cout << "Mass Matrix Operator" << std::endl;
    mat2->print(std::cout);
    */
    //pd_.rightHandSideOperator(level)(*sol,*rhs);
    //pd_.leftHandSideOperator(level)(*rhs,*sol);

    delete temp;
  }

private:
  InverseOperatorFactory<DiscreteFunctionType>& factory_;
  std::auto_ptr<MappingType> inv_;
};


#endif
