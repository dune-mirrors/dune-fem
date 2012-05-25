#include "communicator.hpp"

using namespace pardg;



// class Buffer
Buffer::Buffer() : source(0), dest(0), _size(0), _capacity(0), read_pos(0),
		   data(new char[pre_data_size])
{ 
  assert(data);
}


Buffer::~Buffer()
{
  delete[] data;
}




// class SendBuffer
SendBuffer::SendBuffer() : rcapacity(0)
{}




void SendBuffer::send(MPI_Comm &comm, MPI_Request *send_request,
		      int tag, std::ostream *os)
{
  memcpy(data, &_size, pre_data_size);

  if (_size > rcapacity){
    size_left = _size - rcapacity;
    MPI_Isend(data, pre_data_size+rcapacity, MPI_CHAR, dest, tag, comm, 
	      send_request);

    if (os){
      *os << "send: " << source << "->" << dest << "   "
	  << "1. size: " << rcapacity << "    "
	  << "capacity: " << _capacity << "    "
	  << "rcapacity: " << rcapacity << "    "
	  << std::endl;
    }   
  }
  else{
    size_left = 0;
    MPI_Isend(data, pre_data_size+_size, MPI_CHAR, dest, tag, comm, 
	      send_request);

    if (os){
      *os << "send: " << source << "->" << dest << "   "
	  << "size: " << _size << "    "
	  << "capacity: " << _capacity << "    "
	  << "rcapacity: " << rcapacity << "    "
	  << std::endl;
    }
  }
}


void SendBuffer::send_rest(MPI_Comm &comm, MPI_Request *send_request,
			   int tag, std::ostream *os)
{
  MPI_Isend(data+pre_data_size+rcapacity, size_left, MPI_CHAR, dest, tag, 
	    comm, send_request);

  if (os){
    *os << "send: " << source << "->" << dest << "   "
	<< "2. size: " << size_left << "    "
	<< "capacity: " << _capacity << "    "
	<< "rcapacity: " << _size << "    "
	<< std::endl;
  }
 
  rcapacity = _size;
  size_left = 0;
}



// class ReceiveBuffer
ReceiveBuffer::ReceiveBuffer() {} 


void ReceiveBuffer::receive(MPI_Comm &comm, MPI_Request *receive_request,
			    int tag)
{
  // it doesn't make sense to receive something before
  // everything is read from the buffer
  assert( all_read() );

  read_pos = 0;

  MPI_Irecv(data, pre_data_size+_capacity, MPI_CHAR, source, tag, comm, 
	    receive_request);    
}


bool ReceiveBuffer::rest_to_receive(std::ostream *os)
{
  memcpy(&_size, data, pre_data_size);
  if (_size > _capacity){
    size_left = _size - _capacity;
    resize(_size);
    
    if (os){
      *os << "receive: " << dest << "<-" << source << "   "
	  << "1. size: " << _size - size_left // old capacity
	  << "    "
	  << "capacity: " << _size - size_left
	  << std::endl;
    } 

    return true;
  }
  else{
    size_left = 0;

    if (os){
      *os << "receive: " << dest << "<-" << source << "   "
	  << "size: " << _size << "    "
	  << "capacity: " << _capacity
	  << std::endl;
    } 

    return false;
  }
}


void ReceiveBuffer::receive_rest(MPI_Comm &comm, MPI_Request *receive_request, 
				 int tag, std::ostream *os)
{
  MPI_Irecv(data+pre_data_size+(_capacity-size_left), size_left, MPI_CHAR, 
	    source, tag, comm, receive_request);    

  if (os){
    *os << "receive: " << dest << "<-" << source << "   "
	<< "2. size: " << size_left << "   "
	<< "capacity: " << _capacity
	<< std::endl;
  } 
}


