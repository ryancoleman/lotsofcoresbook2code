/*
 *  $Id: t_spectrum.cc,v 1.17 2005-03-20 18:49:16 edwards Exp $
 *
 *  This is a test program for spectroscopy using qdp++
 *
 *  It is specific to 4 D
 */

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "examples.h"
#include "qdp_util.h"

using namespace QDP;



int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);


  //  Begin by specifying the lattice sizes, of dimension Nd, defined in QDP

  multi1d<int> nsize(Nd);	// Note that the number of dimensions

  // defined in params.h

  //  For those new to C++, note that the vector class contains more 
  //  information than a simple array of elements, and in particular its size

  /*  
      for(int mu = 0; mu < nsize.size(); mu++){
      printf("Specify the size of dimension %d? ",mu);
      scanf("%d", &nsize[mu]);
      }
  */

  const int foo[] = {4, 4, 4, 8};
  nsize = foo;
  // Initialise the layout
  Layout::setLattSize(nsize);
  Layout::create();

  int j_decay = Nd-1;
  int length = Layout::lattSize()[j_decay]; // Define the temporal direction

  XMLFileWriter xml("t_spectrum.xml"); // Open file for sample output
  push(xml,"t_spectrum");


  multi1d<LatticeColorMatrix> u(Nd);

#if 1

  /*
   *  Create a Random gauge field
   */

  for(int mu = 0; mu < u.size(); mu++)
    gaussian(u[mu]);
#endif

#if 0
  /*
   *  Read in a gauge field in the usual NERSC format
   *  GTF: moved filename decl into scope to avoid compiler warning when #if 0
   */

  {
    char filename[100];

    QDPIO::cout << "Reading in NERSC file...";
    printf("Gauge file name? ");
    scanf("%s",filename);
    readArchiv(u, filename);
    QDPIO::cout << "...done\n";

    QDPIO::cout << "Gauge field is ";
    for(int mu = 0; mu < Nd; mu++){
      write(xml, "mu", mu);
      write(xml, "u_mu", u[mu]);
    }
  }
#endif

  /*
   *  Evaluate the plaquette on the gauge field
   */
 {
   Double w_plaq, s_plaq, t_plaq, link;
   QDPIO::cout << "Evaluating the plaquette...";
   MesPlq(u, w_plaq, s_plaq, t_plaq, link);
   QDPIO::cout << "...done\n";

   QDPIO::cout << "w_plaq = " << w_plaq << std::endl;
   QDPIO::cout << "link = " << link << std::endl;
 }





  /*
   *  Now the smeared lattice gauge field, which for the moment we will just
   *  set equal to the original gauge field
   */

  multi1d<LatticeColorMatrix> u_smr(Nd);

  for(int mu = 0; mu < u.size(); mu++)
    u_smr[mu] = u[mu];

  /*
   *  Read in two lattice propagators from disk
   */

  LatticePropagator quark_prop_1;
  gaussian(quark_prop_1);

  multi1d<int> t_source(Nd);	// Source coordinate of propagators
  t_source = 0;

  QDPIO::cout << "Computing simple meson spectroscopy..." << std::endl;

  {
    multi1d< multi1d<Real> > meson_prop;
    mesons(quark_prop_1, quark_prop_1, meson_prop, t_source, j_decay);

    /*
     *  Print the results
     */

    push(xml,"Point_Point_Wilson_Mesons");
    write(xml, "j_decay", j_decay);
    write(xml, "meson_prop", meson_prop);
    pop(xml);
  }

  QDPIO::cout << "...done" << std::endl;

  QDPIO::cout << "Computing simple baryon spectroscopy..." << std::endl;

  {
    multi1d< multi1d<Complex> > baryon_prop;
  
    baryon(quark_prop_1, baryon_prop, t_source, j_decay, 1);

    /*
     *  Print the results
     */

    push(xml,"Point_Point_Wilson_Baryons");
    write(xml, "j_decay", j_decay);
    write(xml, "baryon_prop", baryon_prop);
    pop(xml);
  }

  pop(xml);
  xml.close();
  QDPIO::cout << "...done" << std::endl;

  // Time to bolt
  QDP_finalize();

  exit(0);
}
