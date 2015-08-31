/***********************************************************************
 *   DL_MESO       Version 2.6                                         *
 *   Authors   :   R. S. Qin, M. A. Seaton                             *
 *   Copyright :   STFC Daresbury Laboratory                           *
 *             :   07/08/2014                                          *
 ***********************************************************************
 *   This extract used with permission of the copyright owners         *
 ***********************************************************************
 *   The full DL_MESO package is supplied to individuals under an      *
 *   academic licence, which is free of cost to academic scientists    *
 *   pursuing scientific research of a non-commercial nature           *
 *   To register please visit www.ccp5.ac.uk/DL_MESO                   *
 *   Commercial organisations interested in acquiring the package      *
 *   should approach Dr. M. A. Seaton at Daresbury Laboratory in the   *
 *   first instance. Daresbury Laboratory is the sole centre for       * 
 *   distribution of the package                                       *
 ***********************************************************************/


// Legacy VTK format output



int fsOutputLegacyVTK(const char* filename="lbout")
{
  if(lbsy.nd ==3)
    fsOutputLegacyVTK3D(filename);
  else
    fsOutputLegacyVTK2D(filename);
  return 0;
}



















int fsOutputLegacyVTK3D(const char* filename="lbout")
{
  int i, iprop, j, k;
  float ipos[3];
  int jpos;
  long ilen=0;
  char buf[80];
  sprintf(buf, "%s%.4d.vtk", filename, qVersion);
  ofstream file(buf);
  file<<"# vtk DataFile Version 3.0"<<endl;
  file<<"DL_MESO_LBE"<<endl;
  file<<"ASCII"<<endl;
  file<<"DATASET STRUCTURED_GRID"<<endl;
  file<<"DIMENSIONS  "<<lbsy.nx<<"  "<<lbsy.ny<<"  "<<lbsy.nz<<endl;
  file<<"POINTS  "<<(lbsy.nx*lbsy.ny*lbsy.nz)<<"  float"<<endl;
  // write grid points
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++)
	file<<i*float(lbdx)<<"  "<<j*float(lbdx)<<"  "<<k*float(lbdx)<<endl; 
  file<<endl;
  file<<"POINT_DATA  "<<(lbsy.nx*lbsy.ny*lbsy.nz)<<endl;
  // write densities
  for(iprop=0; iprop<lbsy.nf; iprop++) {
    file<<"SCALARS density_"<<iprop<<" float"<<endl;
    if(iprop==0)
      file<<"LOOKUP_TABLE default"<<endl;
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos[0] = float(fGetOneMassSite(&lbf[ilen*lbsitelength+iprop*lbsy.nq]));
	  file<<ipos[0]<<endl;
	}
      }
    }
  }
  // write mass fractions
  for(iprop=0; iprop<lbsy.nf; iprop++) {
    file<<"SCALARS fraction_"<<iprop<<" float"<<endl;
    if(iprop==0)
      file<<"LOOKUP_TABLE default"<<endl;
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos[0] = float(fGetFracSite(iprop, &lbf[ilen*lbsitelength]));
	  file<<ipos[0]<<endl;
	}
      }
    }
  }
  // write solute concentrations
  for(iprop=0; iprop<lbsy.nc; iprop++) {
    file<<"SCALARS concentration_"<<iprop<<" float"<<endl;
    if(iprop==0)
      file<<"LOOKUP_TABLE default"<<endl;
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos[0] = float(fGetOneConcSite(iprop, ilen));
	  file<<ipos[0]<<endl;
	}
      }
    }
  }
  // write temperatures
  if(lbsy.nt==1) {
    file<<"SCALARS temperature float"<<endl;
    if(iprop==0)
      file<<"LOOKUP_TABLE default"<<endl;
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos[0] = float(fGetTemperatureSite(ilen));
	  file<<ipos[0]<<endl;
	}
      }
    }
  }
  // write velocity
  file<<endl;
  file<<"VECTORS velocity float"<<endl;

  if(!incompress) {
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
 	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13) {
	    ipos[0] = 0.0;
            ipos[1] = 0.0;
            ipos[2] = 0.0;
          }
  	  else {
            ipos[0] = fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
            ipos[1] = fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
            ipos[2] = fGetOneDirecSpeedSite(2, &lbf[ilen*lbsitelength]);
          }
	  file<<ipos[0]<<"  "<<ipos[1]<<"  "<<ipos[2]<<endl;
	}
      }
    }
  }
  else {
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
 	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13) {
	    ipos[0] = 0.0;
            ipos[1] = 0.0;
            ipos[2] = 0.0;
          }
  	  else {
            ipos[0] = fGetOneDirecSpeedIncomSite(0, &lbf[ilen*lbsitelength]);
            ipos[1] = fGetOneDirecSpeedIncomSite(1, &lbf[ilen*lbsitelength]);
            ipos[2] = fGetOneDirecSpeedIncomSite(2, &lbf[ilen*lbsitelength]);
          }
	  file<<ipos[0]<<"  "<<ipos[1]<<"  "<<ipos[2]<<endl;
	}
      }
    }
  }
  file<<endl;

  // write phase field/space property
  file<<"SCALARS phase_field int"<<endl;
  file<<"LOOKUP_TABLE default"<<endl;
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++) 
	{
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  jpos=lbphi[ilen];
	  file<<jpos<<endl;
	}
  file.close();
  return 0;
}




int fsOutputLegacyVTK2D(const char* filename="lbout")
{
  int i, iprop, j;
  float ipos[3];
  int jpos;
  long ilen=0;
  char buf[80];
  sprintf(buf, "%s%.4d.vtk", filename, qVersion);
  ofstream file(buf);
  file<<"# vtk DataFile Version 3.0"<<endl;
  file<<"DL_MESO_LBE"<<endl;
  file<<"ASCII"<<endl;
  file<<"DATASET STRUCTURED_GRID"<<endl;
  file<<"DIMENSIONS  "<<lbsy.nx<<"  "<<lbsy.ny<<"  "<<1<<endl;
  file<<"POINTS  "<<(lbsy.nx*lbsy.ny)<<"  float"<<endl;
  // write grid points
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++)
      file<<i*float(lbdx)<<"  "<<j*float(lbdx)<<"  "<<0.0<<endl; 
  file<<endl;
  file<<"POINT_DATA  "<<(lbsy.nx*lbsy.ny)<<endl;
  // write densities
  for(iprop=0; iprop<lbsy.nf; iprop++) {
    file<<"SCALARS density_"<<iprop<<" float"<<endl;
    if(iprop==0)
      file<<"LOOKUP_TABLE default"<<endl;
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	ipos[0] = float(fGetOneMassSite(&lbf[ilen*lbsitelength+iprop*lbsy.nq]));
	file<<ipos[0]<<endl;
      }
    }
  }
  // write mass fractions
  for(iprop=0; iprop<lbsy.nf; iprop++) {
    file<<"SCALARS fraction_"<<iprop<<" float"<<endl;
    if(iprop==0)
      file<<"LOOKUP_TABLE default"<<endl;
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	ipos[0] = float(fGetFracSite(iprop, &lbf[ilen*lbsitelength]));
	file<<ipos[0]<<endl;
      }
    }
  }
  // write solute concentrations
  for(iprop=0; iprop<lbsy.nc; iprop++) {
    file<<"SCALARS concentration_"<<iprop<<" float"<<endl;
    if(iprop==0)
      file<<"LOOKUP_TABLE default"<<endl;
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	ipos[0] = float(fGetOneConcSite(iprop, ilen));
	file<<ipos[0]<<endl;
      }
    }
  }
  // write temperatures
  if(lbsy.nt==1) {
    file<<"SCALARS temperature float"<<endl;
    if(iprop==0)
      file<<"LOOKUP_TABLE default"<<endl;
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	ipos[0] = float(fGetTemperatureSite(ilen));
	file<<ipos[0]<<endl;
      }
    }
  }

  // write velocity
  file<<endl;
  file<<"VECTORS velocity float"<<endl;
  if(!incompress) {
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid) * lbdm.youter + j + lbdm.bwid;
 	if(lbphi[ilen] == 12 || lbphi[ilen] == 13) {
	  ipos[0] = 0.0;
          ipos[1] = 0.0;
        }
	else {
          ipos[0] = fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
          ipos[1] = fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
        }
	file<<ipos[0]<<"  "<<ipos[1]<<"  "<<0.0<<endl;
      }
    }
  }
  else {
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid) * lbdm.youter + j + lbdm.bwid;
 	if(lbphi[ilen] == 12 || lbphi[ilen] == 13) {
	  ipos[0] = 0.0;
          ipos[1] = 0.0;
        }
	else {
          ipos[0] = fGetOneDirecSpeedIncomSite(0, &lbf[ilen*lbsitelength]);
          ipos[1] = fGetOneDirecSpeedIncomSite(1, &lbf[ilen*lbsitelength]);
        }
	file<<ipos[0]<<"  "<<ipos[1]<<"  "<<0.0<<endl;
      }
    }
  }

  // write phase field/space property
  file<<endl;
  file<<"SCALARS phase_field int"<<endl;
  file<<"LOOKUP_TABLE default"<<endl;
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++) 
      {
	ilen = (i + lbdm.bwid) * lbdm.youter + j + lbdm.bwid;
	jpos=lbphi[ilen];
	file<<jpos<<endl;
	}
  file.close();
  return 0;
}
































