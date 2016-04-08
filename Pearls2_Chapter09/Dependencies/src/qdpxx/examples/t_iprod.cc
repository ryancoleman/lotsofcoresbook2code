#include <qdp.h>
#include <xmmintrin.h>
#define TESTING 1

namespace QDP {

#if 0
  // specialization 
  typedef PColorVector<RComplex<REAL32>,3>  CVec32;
  typedef PScalar<CVec32>  PSCVec32;
  typedef OLattice<PSCVec32 > LCVec32;
  typedef RComplex<REAL32> Cmpx32;


  void local_vcdot_kernel_contig(int lo, int hi, int my_id, VCDotKernelArgs *a); void local_vcdot_kernel_noncontig(int lo, int hi, int my_id, VCDotUnorderedKernelArgs *a);

  void
  innerProductNode( const LCVec32& s1, const LCVec32& s2, const Subset& subset,
		    REAL64 ip[2])
  {
    BinaryReturn< LCVec32, LCVec32, FnInnerProduct>::Type_t ret_val;
    
    
    REAL64 ip_re[qdpNumThreads()];
    REAL64 ip_im[qdpNumThreads()];

    if( subset.hasOrderedRep() ) {
      //    QDPIO::cout << "BJ: innerProduct s" << endl;
      
      // This BinaryReturn has Type_t
      // OScalar<OScalar<OScalar<RComplex<PScalar<REAL> > > > >
      
      unsigned long n_3vec = (subset.end() - subset.start() + 1);
      REAL* V1 =(REAL32* )&(s1.elem(subset.start()).elem().elem(0).real());
      REAL* V2 =(REAL32* )&(s2.elem(subset.start()).elem().elem(0).real());

      VCDotKernelArgs args = {V1,V2,ip_re,ip_im};
      dispatch_to_threads(3*n_3vec,args,local_vcdot_kernel_contig);


      
    
    }
    else {
      int n_3vec = subset.siteTable().size();
      const int* tab = subset.siteTable().slice();
      REAL32* V1 = (REAL32 *)&(s1.elem(0).elem().elem(0).real());
      REAL32* V2 = (REAL32 *)&(s2.elem(0).elem().elem(0).real());

      VCDotUnorderedKernelArgs args = { V1,V2, ip_re, ip_im, tab };
      dispatch_to_threads(n_3vec, args, local_vcdot_kernel_noncontig);

    }

    ip[0] = ip_re[0];
    ip[1] = ip_im[0];
    for(int i=1; i < qdpNumThreads(); i++) { 
      ip[0] += ip_re[i];
      ip[1] += ip_im[i];
    }

  }

  void innerProductMultiKernel(int lo, int hi, int my_id, InnerProdMultiArgs *a);

  //  template<> 
  inline BinaryReturn< LCVec32, LCVec32, FnInnerProductMulti >::Type_t
  innerProductMulti(const LCVec32& s1, const LCVec32& s2, const Set& ss) 
  {
    BinaryReturn<LCVec32,LCVec32, FnInnerProductMulti >::Type_t d(ss.numSubsets());

    int n_color = ss.numSubsets();
    int n_threads=qdpNumThreads();

    REAL32* s1ptr=(REAL32*)&(s1.elem(0).elem().elem(0).real());
    REAL32* s2ptr=(REAL32*)&(s2.elem(0).elem().elem(0).real());

    const int* lat_color=ss.latticeColoring().slice();
    REAL64* ip_flat = (REAL64*)QDP::Allocator::theQDPAllocator::Instance().allocate(2*n_color*n_threads*sizeof(REAL64),QDP::Allocator::DEFAULT);
    
    REAL64* part_results = (REAL64*)QDP::Allocator::theQDPAllocator::Instance().allocate(2*n_color*sizeof(REAL64),QDP::Allocator::DEFAULT);


    InnerProdMultiArgs args={s1ptr,s2ptr,n_color, lat_color, ip_flat};

    dispatch_to_threads(all.siteTable().size(), args, innerProductMultiKernel);


    for(int i=0; i < 2*ss.numSubsets(); i++) {
      part_results[i]=ip_flat[i];
    }

    for(int t=1; t < qdpNumThreads(); t++) {
      for(int i=0; i < 2*ss.numSubsets(); i++) { 
	part_results[i] += ip_flat[2*n_color*t+i];
      }
    }

    for(int i=0; i < ss.numSubsets(); i++) {
      d[i] = cmplx(Real64(part_results[2*i]),Real64(part_results[2*i+1]));
    }

    QDP::Allocator::theQDPAllocator::Instance().free(part_results);
    QDP::Allocator::theQDPAllocator::Instance().free(ip_flat);

    return d;


  }

#endif


};



using namespace QDP;

  namespace
  {
    //! Function object used for constructing the time-slice set
    class TimeSliceFunc : public SetFunc
    {
    public:
      TimeSliceFunc(int dir): dir_decay(dir) {}

      int operator() (const multi1d<int>& coordinate) const ;

      int numSubsets() const ;

    private:
      TimeSliceFunc() {}  // hide default constructor

      int dir_decay;
    };
  }

  int
  TimeSliceFunc::operator() (const multi1d<int>& coordinate) const
  {
    if ((dir_decay<0)||(dir_decay>=Nd)) {
      return 0 ;
    } else {
      return coordinate[dir_decay] ;
    }
  }

  int
  TimeSliceFunc::numSubsets() const
  {
    if ((dir_decay<0)||(dir_decay>=Nd)) {
      return 1 ;
    } else {
      return Layout::lattSize()[dir_decay] ;
    }
  }

int main(int argc, char *argv[])
{

  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the lattice size
  // NOTE: in general, the user will need/want to set the
  // lattice size at run-time
  multi1d<int> nrow(Nd);
  for(int i=0; i < Nd; ++i)
    nrow[i] = 24;         // Set the lattice size to 2^4

  nrow[Nd-1]=128;


  // Insert the lattice size into the Layout
  // There can be additional calls here like setting the number of
  // processors to use in an SMP implementation
  Layout::setLattSize(nrow);

  // Create the lattice layout
  // Afterwards, QDP is useable
  Layout::create();

  Set my_set;
  my_set.make( TimeSliceFunc(3));

  // Do some wonderful and amazing things - impress your friends
  LatticeColorVector s1, s2;
  gaussian(s1);
  gaussian(s2);

  multi1d<ComplexD> iprod1 = sumMulti(localInnerProduct(s1,s2),my_set);

#if 0
  int n_threads=qdpNumThreads();
  int n_color = my_set.numSubsets();


  multi1d<ComplexD> iprod2 = innerProductMulti(s1,s2,my_set);
  for(int i=0; i < iprod1.size(); i++) { 
    QDPIO::cout << "diff["<< i <<"] = " << norm2(iprod1[i]-iprod2[i]) << endl;
  }
  QDPIO::cout << endl;
#endif


  StopWatch swatch;

  swatch.reset();
  swatch.start();

  for(int i=0; i < 500; i++)  {
    iprod1 = sumMulti(localInnerProduct(s1,s2),my_set);
  }
  swatch.stop();
  QDPIO::cout << "Old Way Time = " << swatch.getTimeInSeconds() << endl;

#if 0
  swatch.reset();
  swatch.start();


  for(int i=0; i < 500; i++)  {
    iprod2 = innerProductMulti(s1,s2,my_set);
  }
  swatch.stop();
  QDPIO::cout << "New Way Time = " << swatch.getTimeInSeconds() << endl;
#endif

  // Possibly shutdown the machine
  QDP_finalize();

  exit(0);
}


