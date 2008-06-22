// (c) Dennis Diehl 2007 
#ifndef MPICOMMUNICATOR_HPP
#define MPICOMMUNICATOR_HPP

#include "thread.hpp"
#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <mpi.h>
#include "random.hpp"
#include "timer.hpp"


namespace pardg
{


// forward declaration
class Communicator;


class Buffer
{
public:
  Buffer();
  ~Buffer();

  int size() const; // returns the size without pre_data
  int capacity() const; // returns the capacity without pre_data
  void clear();
  bool empty() const;
  bool all_read() const;

protected:
  void resize(int new_size);
  void put(const char *a, int n);
  void get(char *a, int n);
  
  static const int pre_data_size = sizeof(int);

  int source, dest, _size, _capacity, read_pos, size_left;
  char *data;
};

 

class SendBuffer : public Buffer
{
public:
  SendBuffer();

  // send: size+data
  void send(MPI_Comm &comm, MPI_Request *send_request, 
      int tag, std::ostream *os = NULL);

  // send: data, size is already known
  void send_rest(MPI_Comm &comm, MPI_Request *send_request, 
     int tag, std::ostream *os = NULL);

  bool rest_to_send();

  template<class T> void put(const T& content);
  template<class T> void put(const T* content, int num);


private:
  friend class Communicator;

  int rcapacity;
};





class ReceiveBuffer : public Buffer
{
public:
  ReceiveBuffer();

  void receive(MPI_Comm &comm, MPI_Request *receive_request, int tag);

  void receive_rest(MPI_Comm &comm, MPI_Request *receive_request, 
        int tag, std::ostream *os = NULL);

  bool rest_to_receive(std::ostream *os = NULL);

  template<class T> void get(T& content);
  template<class T> void get(T* content, int num);

private:
  friend class Communicator;
};


class Communicator
{
  Communicator() : 
    comm(MPI_COMM_WORLD), num_send(0), num_recv(0), counter(0), os(NULL)
  {
    MPI_Comm_size(comm, &num_of_processes);
    MPI_Comm_rank(comm, &_id);

    // buffers
    send_buffer = new SendBuffer[num_of_processes];
    receive_buffer = new ReceiveBuffer[num_of_processes];
    sreq = new bool[num_of_processes];
    rreq = new bool[num_of_processes];
    mpi_send_request = new MPI_Request[num_of_processes];
    mpi_receive_request = new MPI_Request[num_of_processes];
    assert(send_buffer && receive_buffer && sreq && rreq
     && mpi_send_request && mpi_receive_request);

    for(int i=0; i<num_of_processes; i++){
      send_buffer[i].source = _id;
      send_buffer[i].dest = i;
      receive_buffer[i].source = i;
      receive_buffer[i].dest = _id;
      sreq[i] = false;
      rreq[i] = false;
    }


    // timers
    timer.start();
    comm_time = 0.0;
    start_time = MPI_Wtime();

    // I/O version
    _io_version = 1;

    // random
    rng.seed(_id+17);
    perm = new int[num_of_processes];
    assert(perm);
    for(int i=0; i<num_of_processes; i++) perm[i] = i;
    rng.permute<int>(perm, num_of_processes);

    // safe communication
    safe_communication = true;

    // tag
    tag = min_tag;
  }
public:
  // static method to give one instance needed
  static Communicator & instance () 
  { 
    int flag;
    flag = 0;
    MPI_Initialized(&flag);
    if (!flag) {
      std::cerr << "the MPI-Communicator in fem/solver/ode "
                << "is used before MPI_Init was called!" << std::endl
                << "Solution 1: configure dune-fem with the " 
                << "--disable-mpi option, or " << std::endl
                << "Solution 2: start your main programm with " << std::endl
                << "   Dune::MPIHelper & mpihelper = "
                << "Dune::MPIHelper::instance(argc,argv);" << std::endl;
      assert(0);
      abort(); 
    }
    static Communicator comm;
    return comm;
  }

  virtual ~Communicator();
  void set_output(std::ostream &os);
  void set_io_version(int version);
  void set_unsafe();
  int id() const;
  int size() const;
  int master() const; // master is the process for special operations
  bool all_read(int i) const;
  bool rb_empty(int i) const; // same as all_read, todo: remove this
  template<class T> void put(int dest, const T& content);
  template<class T> void put_all(const T& content);
  template<class T> void get(int source, T& content);
  template<class T> void put(int dest, const T* content, int num);
  template<class T> void put_all(const T* content, int num);
  template<class T> void get(int source, T* content, int num);

  void send_request(int dest);
  void receive_request(int source);
  void send_receive_request(int process);
  void send_request_all();
  void receive_request_all();  
  void send_receive_request_all();
  virtual void start_communication(const char comment[] = "");
  virtual bool finish_communication();
  bool start_and_finish_communication(const char comment[] = "");
  void barrier();
  void reset_buffers();

  // MPI-wrappers
  void allreduce(int n, double *in, double *out, MPI_Op op);
  void allreduce(int n, int *in, int *out, MPI_Op op);

  // conversion to MPI Communicator
  operator MPI_Comm&(); 
  operator MPI_Comm*(); 

  // Timings
  double time() const;
  void timers(double &run, double &idle, double &comm) const;
  void reset_timers();

  // I/O
  void write(const char filename[]);
  int read(const char filename[], int cycle);

  static const int perm_threshold = 1000;

 protected:
  void write1(const char filename[]);
  int read1(const char filename[], int cycle);
  void write2(const char filename[]);
  int read2(const char filename[], int cycle);
  void inc_tag();

  MPI_Comm comm;
  int _id, num_of_processes, _io_version;

  SendBuffer *send_buffer;
  ReceiveBuffer *receive_buffer;
  bool *sreq, *rreq;
  MPI_Request *mpi_send_request, *mpi_receive_request;
  int num_send, num_recv;

  Timer timer;
  double idle_time, comm_time, comm_start, start_time;
  static const int min_tag = 0;
  static const int max_tag = 10000;
  int counter, tag;
  std::ostream *os;

  Random rng; // random number generator
  int *perm;  // permutation array
  bool safe_communication;
};



// todo: check if the 
// const char *comment 
// can cause a segfault!!!!!!!
class CommunicatorThreaded : public Communicator, public Thread
{
 public:
  CommunicatorThreaded(int argc, char *argv[]);

  // from Communicator
  virtual void start_communication(const char comment[] = "");
  virtual bool finish_communication();

 protected:
  // from Thread
  virtual void run();

 private:
  const char *comment;
};




} // namespace pardg




// class Communicator inline implementation
template<class T> 
inline
void pardg::Communicator::put(int dest, const T& content)
{
  send_buffer[dest].put<T>(content);
}


template<class T> 
inline
void pardg::Communicator::put_all(const T& content)
{
  for(int i=0; i<num_of_processes; i++){
    send_buffer[i].put<T>(content);
  }
}


template<class T> 
inline
void pardg::Communicator::get(int source, T& content)
{
  receive_buffer[source].get<T>(content);
}


template<class T> 
inline 
void pardg::Communicator::put(int dest, const T* content, int num)
{
  send_buffer[dest].put<T>(content, num);
}


template<class T> 
inline
void pardg::Communicator::put_all(const T* content, int num)
{
  for(int i=0; i<num_of_processes; i++){
    send_buffer[i].put<T>(content, num);
  }
}


template<class T> 
inline
void pardg::Communicator::get(int source, T* content, int num)
{
  receive_buffer[source].get<T>(content, num);
}


inline
int pardg::Communicator::id() const
{
  return _id;
}


inline
int pardg::Communicator::size() const
{
  return num_of_processes;
}  


inline
int pardg::Communicator::master() const
{
  return 0;
}  


inline
pardg::Communicator::operator MPI_Comm&()
{
  return comm;
}


inline
pardg::Communicator::operator MPI_Comm*()
{
  return &comm;
}


// wrappers for MPI_Allreduce
inline
void pardg::Communicator::allreduce(int dim, double *in, double *out, 
            MPI_Op op)
{
  double idle_start = time();
  MPI_Allreduce(in, out, dim, MPI_DOUBLE, op, comm);
  idle_time += time() - idle_start;
}


inline
void pardg::Communicator::allreduce(int dim, int *in, int *out, MPI_Op op)
{
  double idle_start = time();
  MPI_Allreduce(in, out, dim, MPI_INT, op, comm);
  idle_time += time() - idle_start;
}


inline
void pardg::Communicator::barrier()
{
  MPI_Barrier(comm);
}



inline
double pardg::Communicator::time() const
{
  return MPI_Wtime();
}


inline
bool pardg::Communicator::all_read(int i) const
{
  return receive_buffer[i].all_read();
}


inline
bool pardg::Communicator::rb_empty(int i) const
{
  return all_read(i);
}


inline
void pardg::Communicator::send_request(int i)
{
  sreq[i] = true;
}


inline
void pardg::Communicator::receive_request(int i)
{
  rreq[i] = true;
}


inline
void pardg::Communicator::send_receive_request(int i)
{
  send_request(i);
  receive_request(i);
}


inline
void pardg::Communicator::send_request_all()
{
  for(int i=0; i<num_of_processes; i++) sreq[i] = true;
}


inline
void pardg::Communicator::receive_request_all()
{
  for(int i=0; i<num_of_processes; i++) rreq[i] = true;
}


inline
void pardg::Communicator::send_receive_request_all()
{
  send_request_all();
  receive_request_all();
}


inline
void pardg::Communicator::inc_tag()
{
  tag++;
  if (tag > max_tag) tag = min_tag;
}




// class Buffer inline implementation
inline
int pardg::Buffer::size() const
{
  return _size;
}


inline
int pardg::Buffer::capacity() const
{
  return _capacity;
}


inline 
void pardg::Buffer::clear()
{
  _size = 0;
  read_pos = 0;
}


inline
bool pardg::Buffer::empty() const
{
  return _size == 0;
}


inline
bool pardg::Buffer::all_read() const
{
  return read_pos == _size;
}


// resize data to pre_data_size+new_capacity
inline
void pardg::Buffer::resize(int new_capacity)
{
  if (false){
    std::cout << "resize " << source << "," << dest << "    "
        << _capacity << " -> " << new_capacity << std::endl;
  }

  char *tmp = new char[pre_data_size + new_capacity];
  assert(tmp);

  if (new_capacity > _capacity){ // enlarge
    memcpy(tmp, data, _capacity+pre_data_size);
  }
  else{ // shrink
    memcpy(tmp, data, new_capacity+pre_data_size);
  }

  delete[] data;
  data = tmp;
  _capacity = new_capacity;
  assert(_size <= _capacity);
}



inline 
void pardg::Buffer::put(const char *a, int n)
{
  // if necessary resize, allocate 20 percent more than necessary
  const double eta = 0.2; 
  if (_size+n > _capacity) resize(static_cast<int>( (1.0+eta) * (_size+n) ) );
  memcpy(data+pre_data_size+_size, a, n);
  _size += n;
}


inline 
void pardg::Buffer::get(char *a, int n)
{
  //std::cout << "get: " << n << std::endl; 
  assert(read_pos+n <= _size);
  memcpy(a, data+pre_data_size+read_pos, n);
  read_pos += n;
}





// class SendBuffer inline implementation
template<class T>
inline
void pardg::SendBuffer::put(const T& content)
{
  pardg::Buffer::put( (const char*)(&content), sizeof(T) );
}


template<class T>
inline
void pardg::SendBuffer::put(const T* content, int num)
{
  pardg::Buffer::put( (const char*)(content), num*sizeof(T) );
}


inline
bool pardg::SendBuffer::rest_to_send()
{
  return (size_left > 0);
}



// class ReceiveBuffer inline implementation
template<class T>
inline
void pardg::ReceiveBuffer::get(T& content)
{
  pardg::Buffer::get( (char*)(&content), sizeof(T) );
}


template<class T>
inline
void pardg::ReceiveBuffer::get(T* content, int num)
{
  pardg::Buffer::get( (char*)(content), num*sizeof(T) );
}





#endif
