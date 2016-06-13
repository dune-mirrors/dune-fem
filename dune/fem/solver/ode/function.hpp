#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <cassert>


namespace pardg
{

class Function
{
public:
  Function();
  virtual ~Function();

  virtual void operator()(const double *u, double *f, int i = 0) = 0;
  virtual int dim_of_argument(int i = 0) const = 0;
  virtual int dim_of_value(int i = 0) const = 0;
  virtual void activateLinear () {}
  virtual void deactivateLinear () {}

  // for time dependent functions
  void operator()(double t, const double *u, double *f, int i = 0);

  double& time();
  double time() const;

  int& flag();
  int flag() const;

private:
  int _flag;
  double _time;
};


} // namespace pardg



// class Function inline implementation

inline
pardg::Function::Function() :
   _flag(0), _time(0.0)
{
}


inline
pardg::Function::~Function()
{
}


inline
void pardg::Function::operator()(double t, const double *u, double *f, int i)
{
  _time = t;
  operator()(u, f, i);
}


inline
double pardg::Function::time() const
{
  return _time;
}


inline
double& pardg::Function::time()
{
  return _time;
}


inline
int pardg::Function::flag() const
{
  return _flag;
}


inline
int& pardg::Function::flag()
{
  return _flag;
}

#endif
