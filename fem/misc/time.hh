#ifndef ADI_TIME_HH
#define ADI_TIME_HH

#include <iostream>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifndef CLK_TCK
#define CLK_TCK 1.
#endif

class rt_struct {
public:
  bool fulloutput;
  double flx,upd,adp,com,all;
  void clear() {flx=0,upd=0,adp=0,com=0,all=0;fulloutput=true;}
  double timediff(struct tms start,struct tms end) {
    clock_t tend   = end.tms_utime; // + end.tms_stime;
    clock_t tstart = start.tms_utime; // + start.tms_stime;
    return ((float)(tend - tstart))/((float)CLK_TCK);
  }
  rt_struct &longinfo() {fulloutput=true;return *this;}
  rt_struct &shortinfo() {fulloutput=false;return *this;}
};


inline std::ostream & operator << (std::ostream &out, const rt_struct& rt) {
  if (rt.fulloutput) {
    out.precision (3) ;
    out << " [com|flx|upd|adp|all] "
     << rt.com << " " << rt.flx << " " << rt.upd << " " << rt.adp << " "
     << rt.all ;
  }
  else {
    out.precision(10);
    out << rt.com << " " << rt.flx << " " << rt.upd << " " << rt.adp  << " "
        << rt.all;
  }
  return out;
} 

#endif
