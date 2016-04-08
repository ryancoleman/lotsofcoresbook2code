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


// XML VTK (structured grid) format output



int fsOutputVTK(const char* filename="lbout")
{
  if(lbsy.nd ==3)
    fsOutputVTK3D(filename);
  else
    fsOutputVTK2D(filename);
  return 0;
}



















int fsOutputVTK3D(const char* filename="lbout")
{
  int i, j, k, iprop;
  float ipos[3];
  int jpos;
  long ilen=0;
  char buf[80];
  sprintf(buf, "%s%.4d.vts", filename, qVersion);
  ofstream file(buf);
  file<<"<?xml version=\"1.0\"?>"<<endl;
  file<<"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">"<<endl;
  file<<"<StructuredGrid WholeExtent=\"";
  file<<"0 "<<(lbsy.nx-1)<<" 0 "<<(lbsy.ny-1)<<" 0 "<<(lbsy.nz-1)<<"\">"<<endl;
  file<<"<Piece Extent=\"0 "<<(lbsy.nx-1)<<" 0 "<<(lbsy.ny-1)<<" 0 "<<(lbsy.nz-1)<<"\">"<<endl;
  // write grid points
  file<<"<Points>"<<endl;
  file<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++)
	file<<i*float(lbdx)<<" "<<j*float(lbdx)<<" "<<k*float(lbdx)<<" ";
  file<<endl;
  file<<"</DataArray>"<<endl;
  file<<"</Points>"<<endl;
  file<<"<PointData Scalars=\"density\" Vectors=\"velocity\">"<<endl;
  // write densities
  for(iprop=0; iprop<lbsy.nf; iprop++) {
    file<<"<DataArray Name=\"density_"<<iprop<<"\" type=\"Float32\" format=\"ascii\">"<<endl;
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos[0] = float(fGetOneMassSite(&lbf[ilen*lbsitelength+iprop*lbsy.nq]));
	  file<<ipos[0]<<" ";
	}
      }
    }
    file<<endl;
    file<<"</DataArray>"<<endl;
  }
  // write mass fractions
  for(iprop=0; iprop<lbsy.nf; iprop++) {
    file<<"<DataArray Name=\"fraction_"<<iprop<<"\" type=\"Float32\" format=\"ascii\">"<<endl;
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos[0] = float(fGetFracSite(iprop, &lbf[ilen*lbsitelength]));
	  file<<ipos[0]<<" ";
	}
      }
    }
    file<<endl;
    file<<"</DataArray>"<<endl;
  }
  // write solute concentrations
  for(iprop=0; iprop<lbsy.nc; iprop++) {
    file<<"<DataArray Name=\"concentration_"<<iprop<<"\" type=\"Float32\" format=\"ascii\">"<<endl;
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos[0] = float(fGetOneConcSite(iprop, ilen));
	  file<<ipos[0]<<" ";
	}
      }
    }
    file<<endl;
    file<<"</DataArray>"<<endl;
  }
  // write temperatures
  if(lbsy.nt==1) {
    file<<"<DataArray Name=\"temperature\" type=\"Float32\" format=\"ascii\">"<<endl;
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos[0] = float(fGetTemperatureSite(ilen));
	  file<<ipos[0]<<" ";
	}
      }
    }
    file<<endl;
    file<<"</DataArray>"<<endl;
  }

  // write velocity
  file<<"<DataArray Name=\"velocity\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
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
  	  file<<ipos[0]<<" "<<ipos[1]<<" "<<ipos[2]<<" ";
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
	  file<<ipos[0]<<" "<<ipos[1]<<" "<<ipos[2]<<" ";
	}
      }
    }
  }
  file<<endl;
  file<<"</DataArray>"<<endl;
  // write phase field/space property
  file<<"<DataArray Name=\"phase_field\" type=\"Int32\" format=\"ascii\">"<<endl;
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++) 
	{
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  jpos=lbphi[ilen];
	  file<<jpos<<" ";
	}

  file<<endl;
  file<<"</DataArray>"<<endl;
  file<<"</PointData>"<<endl;
  file<<"</Piece>"<<endl;
  file<<"</StructuredGrid>"<<endl;
  file<<"</VTKFile>"<<endl;
  file.close();
  return 0;
}




int fsOutputVTK2D(const char* filename="lbout")
{
  int i, j, iprop;
  float ipos[3];
  int jpos;
  long ilen=0;
  char buf[80];
  sprintf(buf, "%s%.4d.vts", filename, qVersion);
  ofstream file(buf);
  file<<"<?xml version=\"1.0\"?>"<<endl;
  file<<"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">"<<endl;
  file<<"<StructuredGrid WholeExtent=\"";
  file<<"0 "<<(lbsy.nx-1)<<" 0 "<<(lbsy.ny-1)<<" 0 0\">"<<endl;
  file<<"<Piece Extent=\"0 "<<(lbsy.nx-1)<<" 0 "<<(lbsy.ny-1)<<" 0 0\">"<<endl;
  // write grid points
  file<<"<Points>"<<endl;
  file<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++)
      file<<i*float(lbdx)<<" "<<j*float(lbdx)<<" "<<0.0<<" ";
  file<<endl;
  file<<"</DataArray>"<<endl;
  file<<"</Points>"<<endl;
  file<<"<PointData Scalars=\"density\" Vectors=\"velocity\">"<<endl;
  // write densities
  for(iprop=0; iprop<lbsy.nf; iprop++) {
    file<<"<DataArray Name=\"density_"<<iprop<<"\" type=\"Float32\" format=\"ascii\">"<<endl;
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	ipos[0] = float(fGetOneMassSite(&lbf[ilen*lbsitelength+iprop*lbsy.nq]));
	file<<ipos[0]<<" ";
      }
    }
    file<<endl;
    file<<"</DataArray>"<<endl;
  }
  // write mass fractions
  for(iprop=0; iprop<lbsy.nf; iprop++) {
    file<<"<DataArray Name=\"fraction_"<<iprop<<"\" type=\"Float32\" format=\"ascii\">"<<endl;
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	ipos[0] = float(fGetFracSite(iprop, &lbf[ilen*lbsitelength]));
	file<<ipos[0]<<" ";
      }
    }
    file<<endl;
    file<<"</DataArray>"<<endl;
  }
  // write solute concentrations
  for(iprop=0; iprop<lbsy.nc; iprop++) {
    file<<"<DataArray Name=\"concentration_"<<iprop<<"\" type=\"Float32\" format=\"ascii\">"<<endl;
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	ipos[0] = float(fGetOneConcSite(iprop, ilen));
	file<<ipos[0]<<" ";
      }
    }
    file<<endl;
    file<<"</DataArray>"<<endl;
  }
  // write temperatures
  if(lbsy.nt==1) {
    file<<"<DataArray Name=\"temperature\" type=\"Float32\" format=\"ascii\">"<<endl;
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	ipos[0] = float(fGetTemperatureSite(ilen));
	file<<ipos[0]<<" ";
      }
    }
    file<<endl;
    file<<"</DataArray>"<<endl;
  }

  // write velocity
  file<<"<DataArray Name=\"velocity\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
  if (!incompress) {
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
	file<<ipos[0]<<" "<<ipos[1]<<" "<<0.0<<" ";
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
	file<<ipos[0]<<" "<<ipos[1]<<" "<<0.0<<" ";
      }
    }
  }

  file<<endl;
  file<<"</DataArray>"<<endl;
  // write phase field/space property
  file<<"<DataArray Name=\"phase_field\" type=\"Int32\" format=\"ascii\">"<<endl;
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++) 
      {
	ilen = (i + lbdm.bwid) * lbdm.youter + j + lbdm.bwid;
	jpos=lbphi[ilen];
	file<<jpos<<" ";
      }
  file<<endl;
  file<<"</DataArray>"<<endl;
  file<<"</PointData>"<<endl;
  file<<"</Piece>"<<endl;
  file<<"</StructuredGrid>"<<endl;
  file<<"</VTKFile>"<<endl;
  file.close();
  return 0;
}


































