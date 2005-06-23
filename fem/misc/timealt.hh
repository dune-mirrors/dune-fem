#ifndef __TIMESTEPPING_HH__
#define __TIMESTEPPING_HH__

//- System includes
#include <memory>
#include <limits>
//#include "description.hh"

//- Local includes
#include "inverseoperatorfactory.hh"

// * temporary include
//#include "../../dune-dxwriter/dxdata.hh"

// * includes for time measurement 
#include "time.hh"


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
    time_(startTime),
    endTime_(std::numeric_limits<double>::max()),
    level_(startLevel) { 
    std::cout << "Loop::Loop: dt = " << dt_ << std::endl;
    setDt(dt);

  }

  virtual ~Loop() {}
  
  // Routines
  void solve();
  void solveUntil(int endLevel);
  void solveUntil(double endTime);
  
  // Accessor functions
  double time() const { return time_; }
  int level() const { return level_; }
  DiscreteFunctionType& solution() { return sol_; }
  
protected:
  virtual void next() = 0;

  void setDt(double dt) {
    // Prevent overshoot at the end of integration
    dt_ = dt;
    pd_.setDt(dt_);
  }

  double getDt() const {
    return pd_.getDt();
  }

  ProblemDescriptionType& pd_;
  DiscreteFunctionType& sol_;
private:
  //std::auto_ptr<TimeStepSizeStrategy> dt_;

  double dt_;
  double time_;
  double endTime_;
  int level_;
};

template <class ProblemDescriptionType>
void Loop<ProblemDescriptionType>::solve() {
  // Preparations (nothing to do right now)

  // Call next()
  next();
  
  // Update time and level
  time_ += getDt();
  ++level_;

  // Output info
  //std::cout << time_ << " time step finished! (dt = " << getDt() << ")" << std::endl;
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
  endTime_ = endTime;
  while (time_ < endTime_) {
    solve();
  }
  endTime_ = std::numeric_limits<double>::max();
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
    factory_(factory),
    temp_(pd.temporary()),
    cfl_(pd.cfl()),
    //    inv_(factory.createOperator(pd.D1()),
    outfile_(pd.runFile().c_str(), std::ios::out)
  {
    std::cout << "ExplicitEuler::ExplicitEuler: dt = " << dt << std::endl;
    // Set up inverse operator
    MappingType& op = pd.D1();
    std::auto_ptr<MappingType> tmp(factory.createOperator(op));
    inv_ = tmp;
  }

  ~ExplicitEuler() {
    // delete temp_;
  }

  // do one timestep 
  void next()
  {
    rt_struct runtimes;
    runtimes.shortinfo();
    runtimes.upd = 0.0;
    struct tms com_start, com_end,
      step_start, step_end, 
      upd_start, upd_end,
      adapt_start, adapt_end,
      all_start, all_end;
    
  // Overall measurement
    times(&all_start);

    times(&com_start);
    // Comm
    this->pd_.communicate();
    times(&com_end);
    runtimes.com = runtimes.timediff(com_start,com_end);
    
    // Numerics
    times(&step_start);
    calculateUpdate();
    times(&step_end);
    runtimes.flx = runtimes.timediff(step_start, step_end);

    // Get (global) timestep size
    this->setDt(cfl_*this->pd_.getTimeStep());

    // Update
    times(&upd_start);
    this->sol_.add(*temp_, this->getDt());
    times(&upd_end);
    runtimes.upd = runtimes.timediff(upd_start, upd_end);

    // Adaptation
    times(&adapt_start);
    this->pd_.markEntities(*temp_);
    this->pd_.adapt();
    times(&adapt_end);
    runtimes.adp = runtimes.timediff(adapt_start, adapt_end);

    // end overall measurement
    times(&all_end);
    runtimes.all = runtimes.timediff(all_start, all_end);

    // Output
    int count = 0;
    typedef typename ProblemDescriptionType::GridType::LeafIterator LeafIterator;
    LeafIterator endit = this->pd_.grid().leafend(this->pd_.grid().maxlevel());
    for (LeafIterator it = 
           this->pd_.grid().leafbegin(this->pd_.grid().maxlevel());
         it != endit; ++it, ++count);

    outfile_ << this->time() << " " << this->getDt() << " "
             << count << " " << runtimes << std::endl;
  }

private:
  // Calculates the update and stores it in DiscreteFunction temp_
  void calculateUpdate() {
    DiscreteFunctionType * rhs = &(this->pd_.rhs(this->time()));
    // * not needed
    //temp_->clear();
    
    // Calculate D_0 u
    this->pd_.D0()(this->sol_, *temp_);
    
    // Get new timestep size
    // * newly outside
    //this->setDt(cfl_*this->pd_.getTimeStep());
    
    // Add rhs and D_0 u
    rhs->add(*temp_, -1.0);
    
    // Invert everything
    (*(this->inv_))(*rhs, *temp_);
  }

  InverseOperatorFactory<DiscreteFunctionType>& factory_;
  std::auto_ptr<MappingType> inv_;
  std::auto_ptr<DiscreteFunctionType> temp_;
  
  double cfl_;
  std::ofstream outfile_;
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

    DiscreteFunctionType * rhs = &(this->pd_.rhs(this->time()+this->getDt()));
    DiscreteFunctionType * temp = this->pd_.temporary();

    // Calculate D_1 u
    this->pd_.D1()(this->sol_, *temp);
   
    // Calculate inv(D_1 + dt D_0)(D_1 u)
    (*(this->inv_))(*temp, this->sol_);

    // Add rhs contribution
    (*(this->inv_))(*rhs, *temp);
    this->sol_.add(*temp, this->getDt());
    
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
