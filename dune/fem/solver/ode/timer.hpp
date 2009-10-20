#ifndef TIMER_HPP
#define TIMER_HPP

#include <unistd.h>
#include <time.h>
#include <sys/times.h>



class Timer
{
public:
  Timer();
  void start();
  void stop();
  void reset();
  double process_time() const;
  double real_time() const;

private:
  struct tms process_start_time, process_elapsed_time;
  time_t real_start_time;
  double real_elapsed_time;
  bool running;
};


inline 
Timer::Timer()
{
  reset();
}


inline 
void Timer::reset()
{
  process_elapsed_time.tms_utime = 0;
  process_elapsed_time.tms_stime = 0;
  process_elapsed_time.tms_cutime = 0;
  process_elapsed_time.tms_cstime = 0;
  real_elapsed_time = 0.0;
  running = false;
}


inline
void Timer::start()
{
  if (!running){
    times(&process_start_time);
    time(&real_start_time);
    running = true;
  }
}


inline
void Timer::stop()
{
  if (running){
    struct tms process_stop_time;
    times(&process_stop_time);
    process_elapsed_time.tms_utime 
      += process_stop_time.tms_utime - process_start_time.tms_utime;
    process_elapsed_time.tms_stime 
      += process_stop_time.tms_stime - process_start_time.tms_stime;
    process_elapsed_time.tms_cutime 
      += process_stop_time.tms_cutime - process_start_time.tms_cutime;
    process_elapsed_time.tms_cstime 
      += process_stop_time.tms_cstime - process_start_time.tms_cstime;
    
    time_t real_stop_time;
    time(&real_stop_time);
    real_elapsed_time += difftime(real_stop_time, real_start_time);

    running = false;
  }
}


// process time in seconds
inline
double Timer::process_time() const
{
  const double clocks_per_second = sysconf(_SC_CLK_TCK);
  
  double total_time = process_elapsed_time.tms_utime 
    +process_elapsed_time.tms_stime
    +process_elapsed_time.tms_cutime 
    +process_elapsed_time.tms_cstime;
 
  if (running){
    struct tms process_stop_time;
    times(&process_stop_time);
    total_time += (process_stop_time.tms_utime - process_start_time.tms_utime)
      +(process_stop_time.tms_stime - process_start_time.tms_stime)
      +(process_stop_time.tms_cutime - process_start_time.tms_cutime)
      +(process_stop_time.tms_cstime - process_start_time.tms_cstime);
  }

  return total_time / clocks_per_second;
}


// real time in seconds
inline
double Timer::real_time() const
{
  double total_time = real_elapsed_time;

  if (running){
    time_t real_stop_time;
    time(&real_stop_time);
    total_time += difftime(real_stop_time, real_start_time);
  }

  return total_time;
}

#endif
