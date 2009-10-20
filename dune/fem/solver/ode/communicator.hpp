#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

#include "thread.hpp"
#include <iostream>
#include <vector>
#include <cassert>

#if HAVE_MPI 
#include <mpi.h>
#include "mpicommunicator.hpp" 
#else 
#include "emptycommunicator.hpp"
#endif

#endif
