// $Id: t_su3.cc,v 1.11 2007-06-13 21:49:51 bjoo Exp $

#include <iostream>
#include <cstdio>

#include <time.h>

#include "qdp.h"

//#define DEBUG_BAGELQDP_LINALG 1
// #include "scalarsite_bagel_qdp/qdp_scalarsite_bagel_qdp_linalg.h"


using namespace QDP;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {6,6,6,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

#if QDP_USE_BAGEL_QDP == 1

  
  LatticeColorMatrix a, b, c, d;
  LatticeColorMatrix a2;
  LatticeColorMatrix diff;

  // -----------------------------------------------------------------
  // MM
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  unsigned long num_sites=all.end()-all.start()+1;
  
  a = b*c;
  qdp_su3_mm(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(c.elem(all.start()).elem().elem(0,0).real()),
	     num_sites, (unsigned long)0);
  
  diff = a - a2;
  
  QDPIO::cout << "MM: || diff || = " << sqrt(norm2(diff)) << endl;
  
  QDP::StopWatch swatch;
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a = b*c;
  }
  swatch.stop();
  double original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_mm(&(a2.elem(all.start()).elem().elem(0,0).real()),
	       &(b.elem(all.start()).elem().elem(0,0).real()),
	       &(c.elem(all.start()).elem().elem(0,0).real()),
	       num_sites, (unsigned long)0);
    
  }
  swatch.stop();
  
  double qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "MM: original  seconds= " << original_secs << endl;
  QDPIO::cout << "MM: bagel_qdp seconds= " << qdp_secs << endl;


  //-----------------------------------------------------------------
  // MA
  //-----------------------------------------------------------------
  gaussian(b);
  gaussian(c);

  a = b*adj(c);
  qdp_su3_ma(&(a2.elem(all.start()).elem().elem(0,0).real()),
  	     &(b.elem(all.start()).elem().elem(0,0).real()),
  	     &(c.elem(all.start()).elem().elem(0,0).real()),
  	     num_sites, (unsigned long)0);

  diff = a - a2;

  QDPIO::cout << "MA: || diff || = " << sqrt(norm2(diff)) << endl;

  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a = b*adj(c);
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();

  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_ma(&(a2.elem(all.start()).elem().elem(0,0).real()),
	       &(b.elem(all.start()).elem().elem(0,0).real()),
	       &(c.elem(all.start()).elem().elem(0,0).real()),
	       num_sites,(unsigned long)0);
  }
  swatch.stop();

  qdp_secs = swatch.getTimeInSeconds();

  QDPIO::cout << "MA: original  seconds= " << original_secs << endl;
  QDPIO::cout << "MA: bagel_qdp seconds= " << qdp_secs << endl;

  //-----------------------------------------------------------------
  // AM
  //-----------------------------------------------------------------
  gaussian(b);
  gaussian(c);

  a = adj(b)*c;
  qdp_su3_am(&(a2.elem(all.start()).elem().elem(0,0).real()),
  	     &(b.elem(all.start()).elem().elem(0,0).real()),
  	     &(c.elem(all.start()).elem().elem(0,0).real()),
  	     num_sites, (unsigned long)0);

  diff = a - a2;

  QDPIO::cout << "AM: || diff || = " << sqrt(norm2(diff)) << endl;

  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a = adj(b)*c;
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();

  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_am(&(a2.elem(all.start()).elem().elem(0,0).real()),
	       &(b.elem(all.start()).elem().elem(0,0).real()),
	       &(c.elem(all.start()).elem().elem(0,0).real()),
	       num_sites, (unsigned long)0);
  }
  swatch.stop();

  qdp_secs = swatch.getTimeInSeconds();

  QDPIO::cout << "AM: original  seconds= " << original_secs << endl;
  QDPIO::cout << "AM: bagel_qdp seconds= " << qdp_secs << endl;

  //-----------------------------------------------------------------
  // AA
  //-----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  double one_minus_i[2] QDP_ALIGN16;
  one_minus_i[0] = (double)1;
  one_minus_i[1] = (double)(-1);

  a = adj(b)*adj(c);
  qdp_su3_aa(&(a2.elem(all.start()).elem().elem(0,0).real()),
  	     &(b.elem(all.start()).elem().elem(0,0).real()),
  	     &(c.elem(all.start()).elem().elem(0,0).real()),
  	     num_sites, (unsigned long)one_minus_i);

  diff = a - a2;

  QDPIO::cout << "AA: || diff || = " << sqrt(norm2(diff)) << endl;

  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a = adj(b)*adj(c);
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();

  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_aa(&(a2.elem(all.start()).elem().elem(0,0).real()),
	       &(b.elem(all.start()).elem().elem(0,0).real()),
	       &(c.elem(all.start()).elem().elem(0,0).real()),
	       num_sites, (unsigned long)one_minus_i);
  }
  swatch.stop();

  qdp_secs = swatch.getTimeInSeconds();

  QDPIO::cout << "AA: original  seconds= " << original_secs << endl;
  QDPIO::cout << "AA: bagel_qdp seconds= " << qdp_secs << endl;


  // -----------------------------------------------------------------
  // MM += 
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);

  Complex two;
  two.elem().elem().elem().real() = 1;
  two.elem().elem().elem().imag() = 0;

  BAGELQDPFloat my_two[2] QDP_ALIGN16;
  my_two[0] =1;
  my_two[1] =0;

  a = d;
  a2 = d;


  a += b*c;
  
  
  
  qdp_su3_mm_peq(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     (BAGELQDPFloat*)my_two,
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(c.elem(all.start()).elem().elem(0,0).real()),
	     num_sites, (unsigned long)0);
  
  
  diff = a - a2;
  

  QDPIO::cout << "MM+=: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a +=  b*c;
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_mm_peq(&(a2.elem(all.start()).elem().elem(0,0).real()),
		   (BAGELQDPFloat *)my_two,
		   &(b.elem(all.start()).elem().elem(0,0).real()),
		   &(c.elem(all.start()).elem().elem(0,0).real()),
		   num_sites, (unsigned long)0);
    
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "MM: original  seconds= " << original_secs << endl;
  QDPIO::cout << "MM: bagel_qdp seconds= " << qdp_secs << endl;

  // -----------------------------------------------------------------
  // MA += 
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);
  a = d;
  a2 = d;

  a += b*adj(c);

  qdp_su3_ma_peq(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
             (BAGELQDPFloat *)my_two,
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(c.elem(all.start()).elem().elem(0,0).real()),
	     num_sites, (unsigned long)0);
  
  diff = a - a2;
  
  QDPIO::cout << "MA+=: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a += b*adj(c);
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_ma_peq(&(a2.elem(all.start()).elem().elem(0,0).real()),
		   (BAGELQDPFloat *)my_two,
		   &(b.elem(all.start()).elem().elem(0,0).real()),
		   &(c.elem(all.start()).elem().elem(0,0).real()),
		   num_sites, (unsigned long)0);
    
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "MA: original  seconds= " << original_secs << endl;
  QDPIO::cout << "MA: bagel_qdp seconds= " << qdp_secs << endl;


  // -----------------------------------------------------------------
  // AM += 
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);
  a = d;
  a2 = d;

  a += adj(b)*c;

  
  
  qdp_su3_am_peq(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     (BAGELQDPFloat *)my_two,
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(c.elem(all.start()).elem().elem(0,0).real()),
	     num_sites, (unsigned long)0);
  
  diff = a - a2;
  
  QDPIO::cout << "AM+=: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a  += adj(b)*c;
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_am_peq(&(a2.elem(all.start()).elem().elem(0,0).real()),
		   (BAGELQDPFloat *)my_two,
		   &(b.elem(all.start()).elem().elem(0,0).real()),
		   &(c.elem(all.start()).elem().elem(0,0).real()),
		   num_sites, (unsigned long)0);
    
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "AM: original  seconds= " << original_secs << endl;
  QDPIO::cout << "AM: bagel_qdp seconds= " << qdp_secs << endl;

  // -----------------------------------------------------------------
  // AA += 
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);

  a=d;
  a2 = d;
  a += adj(b)*adj(c);

  
  
  qdp_su3_aa_peq(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     (BAGELQDPFloat *)my_two,
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(c.elem(all.start()).elem().elem(0,0).real()),
	     num_sites, (unsigned long)one_minus_i);
  
  diff = a - a2;
  
  QDPIO::cout << "AA+=: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a += adj(b)*adj(c);
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_aa_peq(&(a2.elem(all.start()).elem().elem(0,0).real()),
	           (BAGELQDPFloat *)my_two,
		   &(b.elem(all.start()).elem().elem(0,0).real()),
		   &(c.elem(all.start()).elem().elem(0,0).real()),
		   num_sites, (unsigned long)one_minus_i);
    
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "AA: original  seconds= " << original_secs << endl;
  QDPIO::cout << "AA: bagel_qdp seconds= " << qdp_secs << endl;

  // -----------------------------------------------------------------
  // MM -= 
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);

  two.elem().elem().elem().real() = -1;
  two.elem().elem().elem().imag() = 0;

  my_two[0] =-1;
  my_two[1] =0;

  a = d;
  a2 = d;


  a -= b*c;
  
  
  
  qdp_su3_mm_peq(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     (BAGELQDPFloat*)my_two,
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(c.elem(all.start()).elem().elem(0,0).real()),
	     num_sites, (unsigned long)0);
  
  
  diff = a - a2;
  

  QDPIO::cout << "MM-=: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a -=  b*c;
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_mm_peq(&(a2.elem(all.start()).elem().elem(0,0).real()),
		   (BAGELQDPFloat *)my_two,
		   &(b.elem(all.start()).elem().elem(0,0).real()),
		   &(c.elem(all.start()).elem().elem(0,0).real()),
		   num_sites, (unsigned long)0);
    
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "MM: original  seconds= " << original_secs << endl;
  QDPIO::cout << "MM: bagel_qdp seconds= " << qdp_secs << endl;

  // -----------------------------------------------------------------
  // MA -= 
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);
  a = d;
  a2 = d;

  a -= b*adj(c);

  qdp_su3_ma_peq(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
             (BAGELQDPFloat *)my_two,
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(c.elem(all.start()).elem().elem(0,0).real()),
	     num_sites, (unsigned long)0);
  
  diff = a - a2;
  
  QDPIO::cout << "MA-=: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a -= b*adj(c);
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_ma_peq(&(a2.elem(all.start()).elem().elem(0,0).real()),
		   (BAGELQDPFloat *)my_two,
		   &(b.elem(all.start()).elem().elem(0,0).real()),
		   &(c.elem(all.start()).elem().elem(0,0).real()),
		   num_sites, (unsigned long)0);
    
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "MA: original  seconds= " << original_secs << endl;
  QDPIO::cout << "MA: bagel_qdp seconds= " << qdp_secs << endl;


  // -----------------------------------------------------------------
  // AM -= 
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);
  a = d;
  a2 = d;

  a -= adj(b)*c;

  
  
  qdp_su3_am_peq(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     (BAGELQDPFloat *)my_two,
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(c.elem(all.start()).elem().elem(0,0).real()),
	     num_sites, (unsigned long)0);
  
  diff = a - a2;
  
  QDPIO::cout << "AM-=: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a  -= adj(b)*c;
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_am_peq(&(a2.elem(all.start()).elem().elem(0,0).real()),
		   (BAGELQDPFloat *)my_two,
		   &(b.elem(all.start()).elem().elem(0,0).real()),
		   &(c.elem(all.start()).elem().elem(0,0).real()),
		   num_sites, (unsigned long)0);
    
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "AM: original  seconds= " << original_secs << endl;
  QDPIO::cout << "AM: bagel_qdp seconds= " << qdp_secs << endl;

  // -----------------------------------------------------------------
  // AA -= 
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);

  a=d;
  a2 = d;
  a -= adj(b)*adj(c);

  
  
  qdp_su3_aa_peq(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     (BAGELQDPFloat *)my_two,
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(c.elem(all.start()).elem().elem(0,0).real()),
	     num_sites, (unsigned long)one_minus_i);
  
  diff = a - a2;
  
  QDPIO::cout << "AA-=: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a -= adj(b)*adj(c);
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_aa_peq(&(a2.elem(all.start()).elem().elem(0,0).real()),
	           (BAGELQDPFloat *)my_two,
		   &(b.elem(all.start()).elem().elem(0,0).real()),
		   &(c.elem(all.start()).elem().elem(0,0).real()),
		   num_sites, (unsigned long)one_minus_i);
    
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "AA: original  seconds= " << original_secs << endl;
  QDPIO::cout << "AA: bagel_qdp seconds= " << qdp_secs << endl;

  // -----------------------------------------------------------------
  // MM += * 
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);

  two.elem().elem().elem().real() = 2;
  two.elem().elem().elem().imag() = 0;

  my_two[0] =2;
  my_two[1] =0;

  a = d;
  a2 = d;

  Real rtwo = 2;

  a += rtwo*b*c;

  
  
  
  qdp_su3_mm_peq(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     (BAGELQDPFloat*)my_two,
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(c.elem(all.start()).elem().elem(0,0).real()),
	     num_sites, (unsigned long)0);
  
  
  diff = a - a2;
  
  QDPIO::cout << "MM+=*: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a +=  two*b*c;
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_mm_peq(&(a2.elem(all.start()).elem().elem(0,0).real()),
		   (BAGELQDPFloat *)my_two,
		   &(b.elem(all.start()).elem().elem(0,0).real()),
		   &(c.elem(all.start()).elem().elem(0,0).real()),
		   num_sites, (unsigned long)0);
    
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "MM: original  seconds= " << original_secs << endl;
  QDPIO::cout << "MM: bagel_qdp seconds= " << qdp_secs << endl;

  // -----------------------------------------------------------------
  // MM -= * 
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);

  my_two[0] =-2;
  my_two[1] =0;

  a = d;
  a2 = d;

  rtwo = 2;

  a -= rtwo*b*c;

  
  
  
  qdp_su3_mm_peq(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     (BAGELQDPFloat*)my_two,
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(c.elem(all.start()).elem().elem(0,0).real()),
	     num_sites, (unsigned long)0);
  
  
  diff = a - a2;
  
  QDPIO::cout << "MM-=*: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a -=  two*b*c;
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_su3_mm_peq(&(a2.elem(all.start()).elem().elem(0,0).real()),
		   (BAGELQDPFloat *)my_two,
		   &(b.elem(all.start()).elem().elem(0,0).real()),
		   &(c.elem(all.start()).elem().elem(0,0).real()),
		   num_sites, (unsigned long)0);
    
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "MM: original  seconds= " << original_secs << endl;
  QDPIO::cout << "MM: bagel_qdp seconds= " << qdp_secs << endl;

  // -----------------------------------------------------------------
  // M += alpha * M
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);

  a=d;
  a2 = d;
  two.elem().elem().elem().real() = 2;
  my_two[0] = 2;
  rtwo = 2;
  a += rtwo*b;

  
  
  qdp_vaxpy3(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     (BAGELQDPFloat *)my_two,
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     3*num_sites);
  
  diff = a - a2;
  
  QDPIO::cout << "M += a*M: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a += rtwo*b;
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_vaxpy3(
	       &(a2.elem(all.start()).elem().elem(0,0).real()),
	       (BAGELQDPFloat *)my_two,
	       &(b.elem(all.start()).elem().elem(0,0).real()),
	       &(a2.elem(all.start()).elem().elem(0,0).real()),
	       3*num_sites);
    
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "M += a*M: original  seconds= " << original_secs << endl;
  QDPIO::cout << "M += a*M: bagel_qdp seconds= " << qdp_secs << endl;

  // -----------------------------------------------------------------
  // M -= alpha * M
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);

  a=d;
  a2 = d;
  two.elem().elem().elem().real() = 2;
  my_two[0] = -2;
  rtwo = 2;
  a -= rtwo*b;

  
  
  qdp_vaxpy3(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     (BAGELQDPFloat *)my_two,
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     3*num_sites);
  
  diff = a - a2;
  
  QDPIO::cout << "M -= a*M: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a -= rtwo*b;
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_vaxpy3(
	       &(a2.elem(all.start()).elem().elem(0,0).real()),
	       (BAGELQDPFloat *)my_two,
	       &(b.elem(all.start()).elem().elem(0,0).real()),
	       &(a2.elem(all.start()).elem().elem(0,0).real()),
	       3*num_sites);
    
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "M -= a*M: original  seconds= " << original_secs << endl;
  QDPIO::cout << "M -= a*M: bagel_qdp seconds= " << qdp_secs << endl;

  // -----------------------------------------------------------------
  // M += M
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);

  a=d;
  a2 = d;

  a += b;

  
  
  qdp_vadd3(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     3*num_sites);
  
  diff = a - a2;
  
  QDPIO::cout << "M += M: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a += b;
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
    qdp_vadd3(
	      &(a2.elem(all.start()).elem().elem(0,0).real()),
	      &(a2.elem(all.start()).elem().elem(0,0).real()),
	      &(b.elem(all.start()).elem().elem(0,0).real()),
	      3*num_sites);
  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "M += M: original  seconds= " << original_secs << endl;
  QDPIO::cout << "M += M: bagel_qdp seconds= " << qdp_secs << endl;

  // -----------------------------------------------------------------
  // M -= M
  // -----------------------------------------------------------------
  gaussian(b);
  gaussian(c);
  gaussian(d);

  a=d;
  a2 = d;

  a -= b;

  
  
  qdp_vsub3(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     &(a2.elem(all.start()).elem().elem(0,0).real()),	  
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     3*num_sites);
  
  diff = a - a2;
  
  QDPIO::cout << "M -= M: || diff || = " << sqrt(norm2(diff)) << endl;
  
 
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) { 
    a -= b;
  }
  swatch.stop();
  original_secs = swatch.getTimeInSeconds();
  
  swatch.reset();
  swatch.start();
  for(int i=0; i < 5000; i++) {
      qdp_vsub3(
	     &(a2.elem(all.start()).elem().elem(0,0).real()),
	     &(a2.elem(all.start()).elem().elem(0,0).real()),	  
	     &(b.elem(all.start()).elem().elem(0,0).real()),
	     3*num_sites);

  }
  swatch.stop();
  
  qdp_secs = swatch.getTimeInSeconds();
  
  QDPIO::cout << "M -= M: original  seconds= " << original_secs << endl;
  QDPIO::cout << "M -= M: bagel_qdp seconds= " << qdp_secs << endl;

#endif


  // Time to bolt
  QDP_finalize();

  exit(0);
}
