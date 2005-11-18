#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <cassert>

namespace DuneODE {

class Function
{
public:
  Function(int num_of_parameters = 0);
  virtual ~Function();

  virtual void operator()(const double *u, double *f, int i = 0) = 0;
  virtual int dim_of_argument(int i = 0) const = 0;
  virtual int dim_of_value(int i = 0) const = 0;

  // for time dependent functions
  void operator()(double t, const double *u, double *f, int i = 0);

  double& time();
  double time() const;

  int& flag();
  int flag() const;

  double& parameter(int i);
  double parameter(int i) const;

private:
  const int num_of_parameters;
  int _flag;
  double _time, *_parameter;
};



// class Function inline implementation

inline
Function::Function(int num_of_parameters) : 
  num_of_parameters(num_of_parameters), _flag(0), _time(0.0), 
  _parameter(new double[num_of_parameters])
{
  assert(_parameter);
}


inline
Function::~Function()
{
  delete[] _parameter;
}


inline
void Function::operator()(double t, const double *u, double *f, int i)
{
  _time = t;
  operator()(u, f, i);
}


inline
double Function::time() const
{ 
  return _time;
}


inline
double& Function::time()
{ 
  return _time;
}


inline
int Function::flag() const
{ 
  return _flag;
}


inline
int& Function::flag()
{ 
  return _flag;
}


inline
double Function::parameter(int i) const
{
  assert(i >= 0 && i < num_of_parameters);
  return _parameter[i];
}


inline
double& Function::parameter(int i)
{ 
  assert(i >= 0 && i < num_of_parameters);
  return _parameter[i];
}

}




#endif
