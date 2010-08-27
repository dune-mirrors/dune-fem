#include <config.h>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/femtimer.hh>

namespace Dune
{

  namespace Fem
  {

    Timer< true >::Timer ()
    : timesS_(), 
      startTimesV_(),
      timesV_(0),
      timesVName_(0),
      output_(),
      stepCount_(0),
      changed_(true)
    {
      push_time();
    }


    Timer< true >::~Timer ()
    {
      double totalTime = pop_time();

      if( output_.is_open() )
      {
        output_ << "#  ******** TOTAL RUNTIME: " << totalTime
                << "   ******** " << std::endl;
        output_.close();
      }

      const MPIManager::CollectiveCommunication &comm = MPIManager::comm();
      if( comm.rank() == 0 )
      {
        double *totalTimes = new double[ comm.size() ];
        comm.gather( &totalTime, totalTimes, 1, 0 );
        double avgTime = 0.0;
        double minTime = std::numeric_limits< double >::max();
        double maxTime = std::numeric_limits< double >::min();
        for( int i = 0; i < comm.size(); ++i )
        {
          avgTime += totalTimes[ i ];
          minTime = std::min( minTime, totalTimes[ i ] );
          maxTime = std::max( maxTime, totalTimes[ i ] );
        }
        avgTime /= comm.size();
        delete[] totalTimes;

        std::cerr << "#  ******** TOTAL RUNTIME: average = " << avgTime
                  << ", minimum = " << minTime << ", maximum = " << maxTime
                  << "   ******** " << std::endl;
      }
      else comm.gather( &totalTime, (double *)0, 1, 0 );
    }


    unsigned int Timer< true >::add ( const std::string &name, int nr )
    {
      int id;
      for (id=0;id<int(timesV_.size());++id) {
        if (timesV_[id].size()==0) 
          break;
      }
      if (id==int(timesV_.size())) {
        startTimesV_.push_back(std::vector<double>(nr));
        timesV_.push_back(std::vector<double>(nr));
        timesVName_.push_back(name);
      } else {
        startTimesV_[id] = std::vector<double>(nr);
        timesV_[id] = std::vector<double>(nr);
        timesVName_[id] = name;
      }
      reset_timer(id);
      changed_=true;
      return id;
    }


    void Timer< true >::remove ( unsigned int id )
    {
      timesV_[id].clear();
      startTimesV_[id].clear();
      timesVName_[id] = "";
      changed_=true;
    }


    void Timer< true >::remove ()
    {
      timesV_.clear();
      startTimesV_.clear();
      timesVName_.clear();
      changed_=true;
    }


    void Timer< true >::print_timer ( std::ostream &out, int id )
    {
      if (timesV_.size()>1) {
        out << "(" << timesVName_[id] << ":";
      }
      out << timesV_[id][0];
      for (unsigned int i=1;i<timesV_[id].size();++i)
        out << "," << timesV_[id][i]/timesV_[id][0];
      if (timesV_.size()>1) {
        out << ") ";
      }
    }


    void Timer< true >::print_timer ( std::ostream &out, const std::string &msg )
    {
      out << msg << " : ";
      for (unsigned int i=0;i<timesV_.size();++i)
        print_timer(out,i);
      out << std::endl;
    }


    void Timer< true >::printToFile ()
    {
      for (unsigned int i=0;i<timesV_.size();++i) {
        if (timesV_[i].size()>0) {
          output_ << std::setw(6) << inMS(timesV_[i][0]);
          if (timesV_[i].size()>1) {
            output_ << " ( ";
            for (unsigned int nr=1;nr<timesV_[i].size();++nr) {
              output_ << std::setw(3) << inProz(timesV_[i][nr],timesV_[i][0])
                      << "% ";
            }
            output_ << ") ";
          }
        }
      }
      output_ << std::endl;
    }


    void Timer< true >::printToFile ( const std::string &fileName, int step )
    {
      if (!output_.is_open()) {
        output_.open(fileName.c_str());
        if( !output_ )
          DUNE_THROW( IOError, "FemTimer: Unable to open '" << fileName << "' for writing." );
        changed_=true;
      }
      if (changed_) {
        for (unsigned int i=0;i<timesV_.size();++i) {
          if (timesV_[i].size()>0) 
            output_ << std::setw(10+(timesV_[i].size()-1)*5) 
                    << timesVName_[i];
        }
        output_ << std::endl;
        stepCount_=0;
        changed_=false;
      }
      if (stepCount_%step==0) {
        printToFile();
      }
      stepCount_++;
    }


    void Timer< true >::printToFile ( const TimeProviderBase &tp, const std::string &fileName, int step )
    {
      if (!output_.is_open()) {
        output_.open(fileName.c_str());
        if( !output_ )
          DUNE_THROW( IOError, "FemTimer: Unable to open '" << fileName << "' for writing." );
        changed_=true;
      }
      if (changed_) {
        output_ << std::endl << std::endl;
        output_ << std::setw(12) << "Time" << " ";
        output_ << std::setw(12) << "dt" << " ";
        for (unsigned int i=0;i<timesV_.size();++i) {
          if (timesV_[i].size()>0) 
            output_ << std::setw(10+(timesV_[i].size()-1)*5) 
                    << timesVName_[i];
        }
        output_ << std::endl;
        stepCount_=0;
        changed_=false;
      }
      if (stepCount_%step==0) {
        output_ << std::setw(10) << std::scientific << tp.time() << " ";
        output_ << std::setw(10) << std::scientific << tp.deltaT() << " ";
        printToFile();
      }
      stepCount_++;
    }

  }

}
