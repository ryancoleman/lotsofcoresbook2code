// -*- C++ -*-
/*! @file
 * @brief Flop counters
 *
 * Flop counters
 */

#ifndef QDP_FLOPCOUNT_H
#define QDP_FLOPCOUNT_H

namespace QDP
{

  //---------------------------------------------------------------------
  /*! @defgroup qdpflops Flop Counting Mechanism
   *
   * \ingroup qdp
   *
   * This is a basic flop counter class. It is cleared on instantiation
   * but can be cleared at any time by the user. It has functions
   * to add various kinds of flopcount and to retreive the total
   * accumulated flopcount.
   * @{
   */

  //! Basic Flop Counter Clas
  class FlopCounter {
  public:
    //! Constructor - zeroes flopcount
    FlopCounter(void) : count(0), sitesOnNode((unsigned long long)Layout::sitesOnNode()) {}

    //! Destructor - kills object. No cleanup needed
    ~FlopCounter() {} 

    //! Copy Constructor
    FlopCounter(const FlopCounter& c) : count(c.count), sitesOnNode(c.sitesOnNode) {}

    //! Explicit zero method. Clears flopcounts
    inline void reset(void) { 
      count = 0;
    }

    //! Method to add raw number of flops (eg from Level 3 operators)
    inline void addFlops(unsigned long long flops) { 
      count += flops;
    }

    //! Method to add per site flop count. Count is multiplied by sitesOnNode()
    inline void addSiteFlops(unsigned long  long flops) { 
      count += (flops * sitesOnNode);
    }

    //! Method to add per site flop count for a subset of sites. Count is multiplied by the site table size of the subset (ie number of sites in a subset)
    inline void addSiteFlops(unsigned long flops, const Subset& s) {
      count += (flops * (unsigned long long)(s.numSiteTable()));
    }

    //! Method to retrieve accumulated flopcount
    inline unsigned long long getFlops(void) const { 
      return count;
    }

    //! Report floppage
    inline void report(const std::string& name, 
		       const double& time_in_seconds) {

      double mflops_per_cpu = (double)count/((double)(1000*1000)*time_in_seconds);
      double mflops_overall = mflops_per_cpu;
      QDPInternal::globalSum(mflops_overall);
      double gflops_overall = mflops_overall/(double)(1000);
      double tflops_overall = gflops_overall/(double)(1000);

      QDPIO::cout <<"QDP:FlopCount:" << name << " Performance/CPU: t=" << time_in_seconds << "(s) Flops=" << (double)count << " => " << mflops_per_cpu << " Mflops/cpu." << std::endl;
      QDPIO::cout << "QDP:FlopCount:"  << name <<" Total performance:  " << mflops_overall << " Mflops = " << gflops_overall << " Gflops = " << tflops_overall << " Tflops" << std::endl;
    }

  private:
    unsigned long long count;
    const  unsigned long long sitesOnNode;
  };

  /*! @} */  // end of group 

} // namespace QDP

#endif
