#include "communicator.hpp"
#include <fstream>

namespace DuneODE {

Communicator::Communicator(int argc, char *argv[]) : 
  comm(MPI_COMM_WORLD), os(NULL)
{
  // MPI-1
  //MPI_Init(&argc, &argv);

  // MPI-2
  int thread_level;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &thread_level);
  //assert(thread_level == MPI_THREAD_SERIALIZED);

  MPE_Init_log();
  //MPE_Log_event(0, 0, "inittt");

  MPI_Comm_size(comm, &num_of_processes);
  MPI_Comm_rank(comm, &_id);
  send_buffer = new std::vector<char>[num_of_processes];
  receive_buffer = new std::vector<char>[num_of_processes];
  read_pos = new int[num_of_processes];
  receive_buffer_size = new int[num_of_processes];
  sreq = new bool[num_of_processes];
  rreq = new bool[num_of_processes];
 
  assert(send_buffer && receive_buffer && read_pos 
	 && receive_buffer_size && sreq && rreq);

  for(int i=0; i<num_of_processes; i++){
    read_pos[i] = 0;
    receive_buffer_size[i] = 0;
    sreq[i] = false;
    rreq[i] = false;
  }

  // timers
  idle_timer = new double[num_of_processes];
  idle_timer[_id] = 0.0;
  comm_time = 0.0;
  start_time = MPI_Wtime();

  // I/O version
  _io_version = 1;

  // MPE logging
  id_comm_begin = MPE_Log_get_event_number();
  id_comm_end = MPE_Log_get_event_number();
  MPE_Describe_state(id_comm_begin, id_comm_end, "communicate", "blue");
  id_idle_begin = MPE_Log_get_event_number();
  id_idle_end = MPE_Log_get_event_number();
  MPE_Describe_state(id_idle_begin, id_idle_end, "idle", "red");
}


Communicator::~Communicator()
{
  delete[] send_buffer;
  delete[] receive_buffer;
  delete[] read_pos;
  delete[] receive_buffer_size;
  delete[] sreq;
  delete[] rreq;
  delete[] idle_timer;

  MPE_Finish_log("nsk2d_log");
  MPI_Finalize();
}


void Communicator::set_output(std::ostream &os)
{
  this->os = &os;
}


void Communicator::run()
{
  //MPE_Log_event(id_comm_begin, 0, NULL);

  int num_send = 0;
  int num_recv = 0;
  for(int i=0; i<num_of_processes; i++){
    if (sreq[i]) num_send++;
    if (rreq[i]) num_recv++;
  }
		   
  // send
  MPI_Request *send_request = new MPI_Request[num_send];
  assert(send_request);
  
  for(int dest=0, i=0; dest<num_of_processes; dest++) {
    if (sreq[dest]){
      std::vector<char> &buffer = send_buffer[dest];

      MPI_Issend( &(buffer[0]), buffer.size(), 
		 MPI_CHAR, dest, 0, comm, send_request+i );

      if (os){
	*os << "send: " << id() << "->" << dest << "   "
	    << "size: " << buffer.size() 
	    << std::endl;
      }
      i++;
    }
  }

  // probe
  MPI_Status *status = new MPI_Status[num_recv];
  assert(status);
  for(int source=0, i=0; source<num_of_processes; source++) {
    if (rreq[source]){
      MPI_Probe( source, 0, comm, status+i );
      i++;
    }
  }

  // resize buffers & receive
  MPI_Request *receive_request = new MPI_Request[num_recv];
  assert(receive_request);
  for(int source=0, i=0; source<num_of_processes; source++) {
    if (rreq[source]){
      std::vector<char> &buffer = receive_buffer[source];
      int buffer_size;
      MPI_Get_count( status+i, MPI_CHAR, &buffer_size );
      receive_buffer_size[ source ] = buffer_size;
      read_pos[ source ] = 0;
      
      if (buffer_size > buffer.capacity()) buffer.reserve(buffer_size);
 
      MPI_Irecv( &(buffer[0]), buffer_size, MPI_CHAR, 
		 source, 0, comm, receive_request+i );    

      if (os){
	*os << "receive: " << id() << "<-" << source << "   "
	    << "size: " << buffer_size
	    << std::endl;
      } 
      rreq[source] = false;
      i++;
    }
  }

  // wait 
  MPI_Status *send_status = new MPI_Status[num_send];
  MPI_Status *receive_status = new MPI_Status[num_recv];
  assert(send_status && receive_status);
  MPI_Waitall( num_send, send_request, send_status );
  MPI_Waitall( num_recv, receive_request, receive_status );

  // destroy arrays
  delete[] send_request;
  delete[] receive_request;
  delete[] status;
  delete[] send_status;
  delete[] receive_status;

  // clear buffers
  for(int dest=0; dest<num_of_processes; dest++){
    if (sreq[dest]){
      send_buffer[dest].clear();
      sreq[dest] = false;
    }
  }  

  // update timers
  comm_time += MPI_Wtime() - comm_start;
  //MPE_Log_event(id_comm_end, 0, NULL);
}




void Communicator::write2(const char filename[])
{
  MPI_Status status;
  MPI_File file;

  int *length = new int[num_of_processes]; // length in bytes
  assert(length);
  length[_id] = send_buffer[_id].size();
  MPI_Allgather(length+_id, 1, MPI_INT, length, 1, MPI_INT, comm);

  MPI_File_open(comm, const_cast<char*>(filename), 
		MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

  // master process writes header: num_proc, length[num_proc]
  if (_id == master()){
    MPI_File_write(file, &num_of_processes, 1, MPI_INT, &status);
    MPI_File_write(file, length, num_of_processes, MPI_INT, &status);
  }
  
  // write the data
  int offset = (num_of_processes+1) * sizeof(int); // header
  for(int i=0; i<_id; i++) offset += length[i] * sizeof(char);
  MPI_File_seek(file, offset, MPI_SEEK_SET);
  MPI_File_write(file, &send_buffer[_id][0], length[_id], MPI_CHAR, &status);

  if (os){
    *os << "send: " << id() << "->" << filename << "   "
	<< "size: " << length[_id]
	<< std::endl;
  } 
  
  // clear buffer and close file
  send_buffer[_id].clear();
  MPI_File_close(&file);
}


// assumption: file.num_of_processes < num_of_preocesses
// todo: extend this..
void Communicator::read2(const char filename[])
{
  MPI_Status status;
  MPI_File file;

  MPI_File_open(comm, const_cast<char*>(filename), 
		MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
   
  // read the header: num_proc, length[num_proc]
  int old_num_of_processes;
  MPI_File_read(file, &old_num_of_processes, 1, MPI_INT, &status);
  int *length = new int[old_num_of_processes];
  assert(length);
  MPI_File_read(file, length, old_num_of_processes, MPI_INT, &status);
  
  if (_id < old_num_of_processes){
    int offset = (old_num_of_processes+1) * sizeof(int); // header
    for(int i=0; i<_id; i++) offset += length[i] * sizeof(char);
    MPI_File_seek(file, offset, MPI_SEEK_SET);

    // resize the buffer
    std::vector<char> &buffer = receive_buffer[_id];
    if (buffer.size() < length[_id]) buffer.resize(length[_id]);

    MPI_File_read(file, &buffer[0], length[_id], MPI_CHAR, &status);

    if (os){
      *os << "receive: " << id() << "<-" << filename << "   "
	  << "size: " << length[_id]
	  << std::endl;
    } 
  }

  // close file
  MPI_File_close(&file);
  delete[] length;
}





void Communicator::write1(const char filename[])
{
  // todo: only master process needs this
  int *length = new int[num_of_processes]; // length in bytes
  assert(length);
  length[_id] = send_buffer[_id].size();
  MPI_Allgather(length+_id, 1, MPI_INT, length, 1, MPI_INT, comm);


  if (_id != master()){
    // send buffer
    MPI_Send(&send_buffer[_id][0], length[_id], MPI_CHAR, master(), 0,comm);
  }
  else{
    MPI_Status status;
    std::ofstream file(filename);
    
    // buffer
    int max_length = 0;
    for(int i=0; i<num_of_processes; i++){
      if (length[i] > max_length) max_length = length[i];
    }
    char *buffer = new char[max_length];
    assert(buffer);
    

    // write header
    file.write( (char *)(&num_of_processes), sizeof(int) );
    file.write( (char *)(length), sizeof(int)*num_of_processes);

    // receiver and write data
    for(int i=0; i<num_of_processes; i++){
      if (i != _id){
	// receive and write
	MPI_Recv(buffer, length[i], MPI_CHAR, i, 0, comm, &status);
	file.write(buffer, length[i]);
      }
      else{
	// write send_buffer
	file.write(&send_buffer[_id][0], length[_id]);
      }
    }
      
    delete[] buffer;
  }


  if (os){
    *os << "send: " << id() << "->" << filename << "   "
	<< "size: " << length[_id]
	<< std::endl;
  } 

  // clear buffer
  send_buffer[_id].clear();
  delete[] length;
}



// todo:....
void Communicator::read1(const char filename[])
{
  assert(0);
}



void Communicator::write(const char filename[])
{
  if (_io_version == 1) write1(filename);
  else if (_io_version == 2) write2(filename);
  else assert(0);
}


void Communicator::read(const char filename[])
{
  if (_io_version == 1) read1(filename);
  else if (_io_version == 2) read2(filename);
  else assert(0);
}

}
