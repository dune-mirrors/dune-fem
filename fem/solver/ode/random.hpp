#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <cstdlib>

class Random
{
public:
  Random(int seed = 0);
  void seed(int seed);
  int rand();
  int rand(int hi);
  int rand(int lo, int hi);
  double rand(double hi);
  double rand(double lo, double hi);
  int rand_max();
  template<class T> void permute(T *array, int size);

private:
  int _seed;
};


inline
Random::Random(int seed) : _seed(seed)
{
  std::srand(seed);
}


inline
void Random::seed(int seed)
{
  _seed = seed;
  std::srand(_seed);
}


inline
int Random::rand()
{
  return std::rand();
}


inline
int Random::rand_max()
{
  return RAND_MAX;
}


// gives a random number in the range 0, ..., hi-1
inline
int Random::rand(int hi)
{
  return rand() % hi;
}


// gives a random number in the range lo, lo+1, ..., hi-1
inline
int Random::rand(int lo, int hi)
{
  return lo + rand() % (hi-lo);
}


// gives a random number in the range 0, ..., hi-1
inline
double Random::rand(double hi)
{
  return hi * (double)rand() / (double)rand_max();
}


// gives a random number in the range lo, lo+1, ..., hi-1
inline
double Random::rand(double lo, double hi)
{
  return lo + (hi-lo)* (double)rand() / (double)rand_max();
}


// permutates an array of size size
template<class T> 
inline
void Random::permute(T *array, int size)
{
  for(int i=0; i<size; i++){
    int i1 = rand(size);
    T tmp = array[i];
    array[i] = array[i1];
    array[i1] = tmp;      
  }
}


#endif
