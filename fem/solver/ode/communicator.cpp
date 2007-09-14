// (c) Dennis Diehl 2007
#include <fstream>
#include "mpicommunicator.hpp"

namespace pardg {

Communicator::~Communicator()
{
}


void Communicator::set_output(std::ostream &os)
{
  this->os = &os;
}


// version = 1 or 2
void Communicator::set_io_version(int version)
{
  assert(version==1 || version==2);
  _io_version = version;
}


// for unsafe? communication, i.e. without calling this communication
// might be safe
void Communicator::set_unsafe()
{
  safe_communication = false;
}



// must be called on all processes!
void Communicator::reset_buffers()
{
  for(int i=0; i<num_of_processes; i++){
    send_buffer[i].clear();
    send_buffer[i].resize(0);
    send_buffer[i].rcapacity = 0;
    receive_buffer[i].clear();
    receive_buffer[i].resize(0);
  }
}


void Communicator::timers(double &run, double &idle, double &comm) const
{
  const int version = 3;

  if (version == 1){
    run = timer.real_time();
    idle = run - timer.process_time();
    comm = comm_time;
  }
  else if (version == 2){ // works if all cpus are free
    run = time() - start_time;
    idle = idle_time;
    comm = comm_time;
  }
  else if (version == 3){ // experimental
    run = time() - start_time;
    idle = idle_time * timer.process_time() / timer.real_time();
    comm = comm_time;
  }

}


void Communicator::reset_timers()
{
  timer.reset();
  timer.start();
  idle_time = comm_time = 0.0;
  start_time = MPI_Wtime();
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

  if ( !file ){
    std::cerr << __FILE__ << ", " << __LINE__ << ": "
	      << "error writing file: " 
	      << filename
	      << std::endl;
  }

  // master process writes header: num_proc, length[num_proc]
  if (_id == master()){
    MPI_File_write(file, &num_of_processes, 1, MPI_INT, &status);
    MPI_File_write(file, length, num_of_processes, MPI_INT, &status);
  }
  
  // write the data
  int offset = (num_of_processes+1) * sizeof(int); // header
  for(int i=0; i<_id; i++) offset += length[i] * sizeof(char);
  MPI_File_seek(file, offset, MPI_SEEK_SET);
  MPI_File_write(file, send_buffer[_id].data + Buffer::pre_data_size, 
		 length[_id], MPI_CHAR, &status);

  if (os){
    *os << "send: " << id() << "->" << filename << "   "
	<< "size: " << length[_id]
	<< std::endl;
  } 
  
  // clear buffer and close file
  send_buffer[_id].clear();
  MPI_File_close(&file);
}




int Communicator::read2(const char filename[], int cycle)
{
  int exold_num_of_processes, *length;
 
  receive_buffer[_id].clear();

  MPI_Status status;
  MPI_File file;

  MPI_File_open(comm, const_cast<char*>(filename), 
		MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

  if ( !file ){
    std::cerr << __FILE__ << ", " << __LINE__ << ": "
	      << "error reading file: " 
	      << filename
	      << std::endl;
    exit(-1);
  }
   
  // read the header: num_proc, length[num_proc]
  int old_num_of_processes;
  MPI_File_read(file, &old_num_of_processes, 1, MPI_INT, &status);

  if (num_of_processes >= old_num_of_processes){
      exold_num_of_processes = num_of_processes;
  }
  else if (old_num_of_processes % num_of_processes == 0){
    exold_num_of_processes = old_num_of_processes;
  }
  else{
    exold_num_of_processes =
      (old_num_of_processes / num_of_processes + 1) * num_of_processes;
  }

  length = new int[exold_num_of_processes];
  assert(length);
  MPI_File_read(file, length, old_num_of_processes, MPI_INT, &status);
  
  // fill rest of length with zeros
  for(int i=old_num_of_processes; i<exold_num_of_processes; i++){
    length[i] = 0;
  }

  // determine file offset
  const int index = cycle * num_of_processes + _id;
  int offset = (old_num_of_processes+1) * sizeof(int); // header
  for(int i=0; i<index; i++) offset += length[i] * sizeof(char);
  MPI_File_seek(file, offset, MPI_SEEK_SET);

  // resize the buffers
  if (receive_buffer[_id].capacity() < length[index]){
    receive_buffer[_id].resize(length[index]);
    send_buffer[_id].rcapacity = receive_buffer[_id].capacity();
  }

  // read content from file
  MPI_File_read(file, receive_buffer[_id].data+Buffer::pre_data_size, 
		length[index], MPI_CHAR, &status);
  receive_buffer[_id]._size = length[index];

  if (os){
    *os << "receive: " << id() << "<-" << filename << "   "
	<< "size: " << length[index]
	<< std::endl;
  } 


  // close file
  MPI_File_close(&file);
  delete[] length;

  // return num of cycles
  return exold_num_of_processes / num_of_processes;
}




void Communicator::write1(const char filename[])
{
  int *length = new int[num_of_processes]; // length in bytes
  assert(length);
  length[_id] = send_buffer[_id].size();
  MPI_Allgather(length+_id, 1, MPI_INT, length, 1, MPI_INT, comm);


  if (_id != master()){
    // send buffer
    MPI_Send(send_buffer[_id].data + Buffer::pre_data_size, 
	     length[_id], MPI_CHAR, master(), tag, comm);
  }
  else{
    MPI_Status status;
    std::ofstream file(filename);

    if ( !file ){
      std::cerr << __FILE__ << ", " << __LINE__ << ": "
		<< "error writing file: " 
		<< filename
		<< std::endl;
    }
    
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
	MPI_Recv(buffer, length[i], MPI_CHAR, i, tag, comm, &status);
	file.write(buffer, length[i]);
      }
      else{
	// write send_buffer
	file.write(send_buffer[_id].data + Buffer::pre_data_size, length[_id]);
      }
    }
      
    delete[] buffer;
  }

  // increment of tag
  inc_tag();


  if (os){
    *os << "send: " << id() << "->" << filename << "   "
	<< "size: " << length[_id]
	<< std::endl;
  } 

  // clear buffer
  send_buffer[_id].clear();
  delete[] length;
}






// return num_cycles
int Communicator::read1(const char filename[], int cycle)
{
  int exold_num_of_processes, *length;

  receive_buffer[_id].clear();

  if (_id == master()){ // master process
    // open file
    std::ifstream file(filename);
    if ( !file ){
      std::cerr << __FILE__ << ", " << __LINE__ << ": "
		<< "error reading file: " 
		<< filename
		<< std::endl;
      exit(-1);
    }

    // read header
    int old_num_of_processes;
    file.read( (char *)(&old_num_of_processes), sizeof(int) );

    if (num_of_processes >= old_num_of_processes){
      exold_num_of_processes = num_of_processes;
    }
    else if (old_num_of_processes % num_of_processes == 0){
      exold_num_of_processes = old_num_of_processes;
    }
    else{
      exold_num_of_processes = 
	(old_num_of_processes / num_of_processes + 1) * num_of_processes;
    }

    length = new int[exold_num_of_processes];
    assert(length);
    file.read( (char *)(length), sizeof(int)*old_num_of_processes );


    // fill rest of length with zeros
    for(int i=old_num_of_processes; i<exold_num_of_processes; i++){
      length[i] = 0;
    }    

    // determine max_length for file-read-buffer
    int max_length = 0;
    for(int i=0; i<old_num_of_processes; i++){
      if (max_length < length[i]) max_length = length[i];
    }

    // send the complete header information+extension to other processes
    for(int i=0; i<num_of_processes; i++){
      if( i != master() ){
	put<int>(i, exold_num_of_processes);
	put<int>(i, length, exold_num_of_processes);
	send_request(i);
      }
    }
    start_and_finish_communication("Communicator::read1");

    // allocate buffer memory
    char *buffer = new char[max_length];
    assert(buffer);

    // determine file offset and seek to the correct position for
    // this cycle
    int offset = (1 + old_num_of_processes) * sizeof(int); // header
    for(int k=0; k<cycle*num_of_processes; k++) offset += length[k];
    file.seekg(offset);

    // read and send
    for(int i=0; i<num_of_processes; i++){
      const int index = cycle * num_of_processes + i;

      if (i == _id){ // dest is master
	if (receive_buffer[_id].capacity() < length[index]){
	  receive_buffer[_id].resize(length[index]);
	  send_buffer[_id].rcapacity = receive_buffer[_id].capacity();
	}
	file.read(receive_buffer[_id].data+Buffer::pre_data_size,
		  length[index]);
	receive_buffer[_id]._size = length[index];
      }

      else{ // dest is other process
	file.read(buffer, length[index]);
	MPI_Send(buffer, length[index], MPI_CHAR, i, tag, comm);
      }
    }

    // free mem
    delete[] buffer;
  }
  else{ // other processes
    // receive header information + extension
    receive_request( master() );
    start_and_finish_communication();
    get<int>( master(), exold_num_of_processes );
    length = new int[exold_num_of_processes];
    assert(length);
    get<int>( master(), length, exold_num_of_processes );    

    // resize receive_buffer and receive
    const int index = cycle*num_of_processes + _id;
    if (receive_buffer[_id].capacity() < length[index]) {
      receive_buffer[_id].resize(length[index]);
      send_buffer[_id].rcapacity = receive_buffer[_id].capacity();
    }

    MPI_Status status;
    MPI_Recv(receive_buffer[_id].data + Buffer::pre_data_size, 
	     length[index], MPI_CHAR, master(), tag, comm, &status);

    receive_buffer[_id]._size = length[index];
  }

  // increment of tag
  inc_tag();


  // output
  if (os){
    *os << "receive: " << _id << "<-" << filename << "   "
	<< "size: " << receive_buffer[_id].size()
	<< std::endl;
  } 

  // free mem
  delete[] length;

  // return num of cycles
  return exold_num_of_processes / num_of_processes;
}



void Communicator::write(const char filename[])
{
  barrier();
  if (_io_version == 1) write1(filename);
  else if (_io_version == 2) write2(filename);
  else assert(0);
  barrier();
}




int Communicator::read(const char filename[], int cycle)
{
  barrier();
  if (_io_version == 1) return read1(filename, cycle);
  else if (_io_version == 2) return read2(filename, cycle);
  else{
    assert(0);
    return 0;
  }
  barrier();
}



bool Communicator::start_and_finish_communication(const char comment[])
{
  start_communication(comment);
  return finish_communication();
}




void Communicator::start_communication(const char comment[])
{
  // repermutate 
  if (counter % perm_threshold == 0) rng.permute(perm, num_of_processes);

  if (os){
    *os << std::endl
	<< "communication: " << counter << "  "
	<< "tag: " << tag << "  "
	<< comment
	<< std::endl;
  }
  counter++;

  double idle_start = time();
 
  //MPI_Barrier(comm); // ????????????????

  idle_time += time() - idle_start;
  comm_start = MPI_Wtime();

  num_send = 0;
  num_recv = 0;
  for(int i=0; i<num_of_processes; i++){
    if (sreq[i]) num_send++;
    if (rreq[i]) num_recv++;
  }
	

  // receive
  for(int source=0, i=0; source<num_of_processes; source++) {
    const int sourcep = perm[source];
    if (rreq[sourcep]){
      receive_buffer[sourcep].receive(comm, mpi_receive_request+i, tag);
      i++;
    }
  }	   


  // send
  for(int dest=0, i=0; dest<num_of_processes; dest++) {
    const int destp = perm[dest];
    if (sreq[destp]){
      send_buffer[destp].send(comm, mpi_send_request+i, tag, os);
      i++;
    }
  }

  // increment of tag
  inc_tag();
}
 

bool Communicator::finish_communication()
{
  // repermutate 
  if (counter % perm_threshold == 0) rng.permute(perm, num_of_processes);

  // wait 
  MPI_Status *send_status = new MPI_Status[num_send];
  MPI_Status *receive_status = new MPI_Status[num_recv];
  assert(send_status && receive_status);
  double idle_start = time();
  MPI_Waitall( num_send, mpi_send_request, send_status );
  MPI_Waitall( num_recv, mpi_receive_request, receive_status );
  idle_time += time() - idle_start;

  // receive_rest
  num_recv = 0;
  for(int source=0; source<num_of_processes; source++){
    const int sourcep = perm[source];
    if (rreq[sourcep]){
      if (receive_buffer[sourcep].rest_to_receive(os)){
	receive_buffer[sourcep].
	  receive_rest(comm, mpi_receive_request+num_recv, tag, os);
	num_recv++;
      }

      rreq[sourcep] = false;
    }
  }

  // send_rest
  num_send = 0;
  for(int dest=0; dest<num_of_processes; dest++){
    const int destp = perm[dest];
    if (sreq[destp]){
      if (send_buffer[destp].rest_to_send()){
	send_buffer[destp].send_rest(comm, mpi_send_request+num_send, tag, os);
	num_send++;
      }
    }
  }

  // increment of tag
  inc_tag();
  
  // wait for rest
  idle_start = time();
  MPI_Waitall( num_send, mpi_send_request, send_status );
  MPI_Waitall( num_recv, mpi_receive_request, receive_status );
  idle_time += time() - idle_start;


  // destroy arrays
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

  return true;
}

}

