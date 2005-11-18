#ifndef LIMITER_HPP
#define LIMITER_HPP

namespace DuneODE {

class Limiter
{
public:
  virtual void operator()(double *u) = 0;
};


}
#endif
