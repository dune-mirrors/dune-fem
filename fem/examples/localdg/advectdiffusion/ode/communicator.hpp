#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

#include "thread.hpp"
#include <iostream>
#include <vector>
#include <cassert>
// #include <pthread.h>
// #include <mpe.h>

namespace DuneODE {
#include "mpi.h"
  class Communicator {
  public:
    Communicator() {}
    Communicator(int argc, char *argv[]) {}
    ~Communicator() {}
    void set_output(std::ostream &os) {}
    int id() const {return 0;}
    int size() const {return 1;}
    int master() const {return 0;}
    bool rb_empty(int i) {return true;}
    template<class T> void put(int dest, const T& content) {}
    template<class T> void put_all(const T& content) {}
    template<class T> bool get(int source, T& content) {return true;}
    template<class T> void put(int dest, const T* content, int num) {}
    template<class T> void put_all(const T* content, int num) {}
    template<class T> bool get(int source, T* content, int num) {return true;}
    
    void send_request(int dest) {}
    void receive_request(int source) {}
    void send_receive_request(int process) {}
    void send_request_all() {}
    void receive_request_all() {}
    void send_receive_request_all() {}
    void start_communication() {}
    bool finish_communication() {return true;}
    bool start_and_finish_communication() {return true;}
    void barrier() {}
    
    // MPI-wrapper
    void allreduce(int n, double *in, double *out, MPI_Op op) const {
      for (int i=0;i<n;i++)
        out[i]=in[i];
    }
    
    // conversion to MPI Communicator
    operator MPI_Comm() const {}
    
    // Timings
    double time() const {return 0.;}
    double communication_time() const {return 0.;}
    double run_time() const {return 0.;}
    const double* idle_times() const {return 0;}
    void reset_timers() {}
    
    // I/O
    void write(const char filename[]) {}
    void read(const char filename[]) {}
  };

/*
//extern "C"{ void* communicate(void *ptr);}


class Communicator : public Thread
{
public:
  Communicator(int argc, char *argv[]);
  ~Communicator();
  void set_output(std::ostream &os);
  int id() const;
  int size() const;
  int master() const; // master is the process for special operations
  bool rb_empty(int i) const;
  template<class T> void put(int dest, const T& content);
  template<class T> void put_all(const T& content);
  template<class T> bool get(int source, T& content);
  template<class T> void put(int dest, const T* content, int num);
  template<class T> void put_all(const T* content, int num);
  template<class T> bool get(int source, T* content, int num);

  void send_request(int dest);
  void receive_request(int source);
  void send_receive_request(int process);
  void send_request_all();
  void receive_request_all();  
  void send_receive_request_all();
  void start_communication();
  bool finish_communication();
  bool start_and_finish_communication();
  void barrier();

  // MPI-wrapper
  void allreduce(int n, double *in, double *out, MPI_Op op) const;

  // conversion to MPI Communicator
  operator MPI_Comm() const; 

  // Timings
  double time() const;
  double communication_time() const;
  double run_time() const;
  const double* idle_times() const;
  void reset_timers();

  // I/O
  void write(const char filename[]);
  void read(const char filename[]);

protected:
  // from Thread
  virtual void run();
private:
  void write1(const char filename[]);
  void read1(const char filename[]);
  void write2(const char filename[]);
  void read2(const char filename[]);

  MPI_Comm comm;
  int _id, num_of_processes, _io_version;
  std::vector<char> *send_buffer, *receive_buffer;
  int *receive_buffer_size;
  bool *sreq, *rreq;
  pthread_t comm_thread;
  int *read_pos;
  double *idle_timer, comm_time, comm_start, start_time;
  std::ostream *os;

  // logging
  int id_comm_begin, id_comm_end, id_idle_begin, id_idle_end;
};




inline
int Communicator::id() const
{
  return _id;
}


inline
int Communicator::size() const
{
  return num_of_processes;
}  


inline
int Communicator::master() const
{
  return 0;
}  


inline
Communicator::operator MPI_Comm() const
{
  return comm;
}


// wrapper for MPI_Allreduce
inline
void Communicator::allreduce(int dim, double *in, 
			     double *out, MPI_Op op) const
{
  double idle_start = time();
  MPI_Allreduce(in, out, dim, MPI_DOUBLE, op, comm);
  idle_timer[_id] += time() - idle_start;
}


inline
double Communicator::time() const
{
  return MPI_Wtime();
}


inline
double Communicator::communication_time() const
{
  return comm_time;
}


inline 
double Communicator::run_time() const
{
  return (MPI_Wtime() - start_time);
}


inline 
const double* Communicator::idle_times() const
{
  MPI_Allgather(idle_timer+_id,1, MPI_DOUBLE, idle_timer,1, MPI_DOUBLE, comm);
  return idle_timer;
}


inline
void Communicator::reset_timers()
{
  idle_timer[_id] = 0.0;
  comm_time = 0.0;
  start_time = MPI_Wtime();
}


template<class T>
inline
void Communicator::put(int dest, const T& content)
{
  if ( !(dest>=0 && dest<num_of_processes) ){
    std::cout << "dest = " << dest << std::endl;
  }
  assert(dest>=0 && dest<num_of_processes);
  char *c = (char *)(&content);
  for(int i=0; i<sizeof(T); i++) send_buffer[dest].push_back(c[i]);
}


template<class T>
inline
void Communicator::put(int dest, const T* content, int num)
{
  assert(dest>=0 && dest<num_of_processes);
  char *c = (char *)content;
  for(int i=0; i<sizeof(T)*num; i++) send_buffer[dest].push_back(c[i]);
}


template<class T>
inline
void Communicator::put_all(const T& content)
{
  char *c = (char *)(&content);
  for(int k=0; k<num_of_processes; k++){
    for(int i=0; i<sizeof(T); i++) send_buffer[k].push_back(c[i]);
  }
}


template<class T>
inline
void Communicator::put_all(const T* content, int num)
{
  char *c = (char *)content;
  for(int k=0; k<num_of_processes; k++){
    for(int i=0; i<sizeof(T)*num; i++) send_buffer[k].push_back(c[i]);
  }
}


// todo: return
template<class T>
inline
bool Communicator::get(int source, T& content)
{
  assert(source>=0 && source<num_of_processes);
  int size = sizeof(T);

  if ( read_pos[source]+size <= receive_buffer_size[source] ){ 
    memcpy(&content, &(receive_buffer[source][ read_pos[source] ]), size);
    read_pos[source] += size;
    return true;
  }
  else return false;
}


template<class T>
inline
bool Communicator::get(int source, T* content, int num)
{
  assert(source>=0 && source<num_of_processes);
  int size = sizeof(T)*num;

  if ( read_pos[source]+size <= receive_buffer_size[source] ){ 
    memcpy(content, &(receive_buffer[source][ read_pos[source] ]), size);
    read_pos[source] += size;
    return true;
  }
  else return false;
}


inline
bool Communicator::rb_empty(int i) const
{
  return (read_pos[i] >= receive_buffer_size[i]);
}


inline
void Communicator::send_request(int i)
{
  sreq[i] = true;
}


inline
void Communicator::receive_request(int i)
{
  rreq[i] = true;
}


inline
void Communicator::send_receive_request(int i)
{
  send_request(i);
  receive_request(i);
}


inline
void Communicator::send_request_all()
{
  for(int i=0; i<num_of_processes; i++) sreq[i] = true;
}


inline
void Communicator::receive_request_all()
{
  for(int i=0; i<num_of_processes; i++) rreq[i] = true;
}


inline
void Communicator::send_receive_request_all()
{
  send_request_all();
  receive_request_all();
}


inline
void Communicator::start_communication()
{
  MPE_Log_event(id_idle_begin, 0, NULL);
  double idle_start = time();
  MPI_Barrier(comm);
  MPE_Log_event(id_idle_end, 0, NULL);
  idle_timer[_id] += time() - idle_start;
  comm_start = MPI_Wtime();
  start();
}


// todo: return true if successful
inline
bool Communicator::finish_communication()
{
  MPE_Log_event(id_idle_begin, 0, NULL);
  double idle_start = time();
  join();
  MPE_Log_event(id_idle_end, 0, NULL);
  idle_timer[_id] += time() - idle_start;
  return true;
}


inline
bool Communicator::start_and_finish_communication()
{
  start_communication();
  return finish_communication();
}


inline
void Communicator::barrier()
{
  if (os){
    *os << "======================================" << std::endl;
  }
  MPI_Barrier(comm);
}
*/

}


#endif
