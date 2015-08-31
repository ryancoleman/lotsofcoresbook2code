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


// Plot3D format output


int fsOutputGrid(const char* filename="lbout")
{
  if(lbsy.nd == 3)
    fsOutputGrid3D(filename);
  else
    fsOutputGrid2D(filename);
  return 0;
}


int fsOutputQ(const char* filename="lbout")
{
  char namebuf[80];
  if(lbsy.nd==3) {
    for(int iprop=0; iprop<lbsy.nf; iprop++) {
      sprintf(namebuf, "%s%.2ddens", filename, iprop);
      fsOutputQP3D(namebuf, iprop);
      sprintf(namebuf, "%s%.2dfrac", filename, iprop);
      fsOutputQCA3D(namebuf, iprop);
    }
    for(int iprop=0; iprop<lbsy.nc; iprop++) {
      sprintf(namebuf, "%s%.2dconc", filename, iprop);
      fsOutputQCB3D(namebuf, iprop);
    }
    if(lbsy.nt==1) {
      sprintf(namebuf, "%stemp", filename);
      fsOutputQT3D(namebuf);
    }
  }
  else {
    for(int iprop=0; iprop<lbsy.nf; iprop++) {
      sprintf(namebuf, "%s%.2ddens", filename, iprop);
      fsOutputQP2D(namebuf, iprop);
      sprintf(namebuf, "%s%.2dfrac", filename, iprop);
      fsOutputQCA2D(namebuf, iprop);
    }
    for(int iprop=0; iprop<lbsy.nc; iprop++) {
      sprintf(namebuf, "%s%.2dconc", filename, iprop);
      fsOutputQCB2D(namebuf, iprop);
    }
    if(lbsy.nt==1) {
      sprintf(namebuf, "%stemp", filename);
      fsOutputQT2D(namebuf);
    }
  }
  return 0;
}





















int fsOutputGrid3D(const char* filename="lbout")
{
  int i, j, k;
  char buf[80];
  sprintf(buf, "%s.xyz", filename);
  ofstream file(buf);
  file<<lbsy.nx<<" "<<lbsy.ny<<" "<<lbsy.nz<<" "<<endl;
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++)
	file<<i*float(lbdx)<<" "; 
  file<<endl;
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++)
	file<<j*float(lbdx)<<" "; 
  file<<endl;
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++)
	file<<k*float(lbdx)<<" ";  
  file<<endl;
  file.close();
  return 0;
}


int fsOutputGrid2D(const char* filename="lbout")
{
  int i, j;
  char buf[80];
  sprintf(buf, "%s.xy", filename);
  ofstream file(buf);
  file<<lbsy.nx<<" "<<lbsy.ny<<" "<<endl;
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++)
      file<<i*float(lbdx)<<" "; 
  file<<endl;
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++)
      file<<j*float(lbdx)<<" "; 
  file<<endl;
  file.close();
  return 0;
}




int fsOutputQP3D(const char* filename="lbout", int iprop=0)
{
  int i, j, k;
  float ipos=0.0;
  long ilen=0;
  char buf[80];
  sprintf(buf, "%s%.4d.q", filename, qVersion);
  ofstream file(buf);
  file<<lbsy.nx<<" "<<lbsy.ny<<" "<<lbsy.nz<<" "<<endl;
  file<<lbsoundv<<" "<<1.0<<" "<<lbreynolds<<" "<<qVersion<<" "<<endl;
  // write density
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++) 
	{
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
          ipos = float(fGetOneMassSite(&lbf[ilen*lbsitelength+iprop*lbsy.nq]));
	  file<<ipos<<" ";
	}

  if(!incompress) {
  // write x-component of velocity (v_x)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write y-component of velocity (v_y)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write z-component of velocity (v_z)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(2, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  }
  else {
  // write x-component of velocity (v_x)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write y-component of velocity (v_y)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write z-component of velocity (v_z)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(2, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  }

  // write phase field/space property
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++) 
	{
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos=(float)lbphi[ilen];
	  file<<ipos;
          if(ilen<(lbsy.nx*lbsy.ny*lbsy.nz-1))
            file<<" ";
	}
  file<<endl;
  file.close();
  return 0;
}




int fsOutputQP2D(const char* filename="lbout", int iprop=0)
{
  int i, j;
  float ipos=0.0;
  long ilen=0;
  char buf[80];
  sprintf(buf, "%s%.4d.q", filename, qVersion);
  ofstream file(buf);
  file<<lbsy.nx<<" "<<lbsy.ny<<endl;
  file<<lbsoundv<<" "<<1.0<<" "<<lbreynolds<<" "<<qVersion<<" "<<endl;
  // write density
  for(j=lbdm.bwid; j<lbdm.youter-lbdm.bwid; j++)
    for(i=lbdm.bwid; i<lbdm.xouter-lbdm.bwid; i++) 
      {
	ilen = i* lbdm.youter + j;
	ipos = float(fGetOneMassSite(&lbf[ilen*lbsitelength+iprop*lbsy.nq]));
	file<<ipos<<" ";
      }

  if(!incompress) {
  // write x-component of velocity (v_x)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  // write y-component of velocity (v_y)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  }
  else {
  // write x-component of velocity (v_x)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  // write y-component of velocity (v_y)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  }

  // write phase field/space property
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++) 
      {
	ilen = (i + lbdm.bwid) * lbdm.youter + j + lbdm.bwid;
	ipos=(float)lbphi[ilen];
	file<<ipos;
        if(ilen<(lbsy.nx*lbsy.ny-1))
          file<<" ";
      }
  file<<endl;
  file.close();
  return 0;
}




int fsOutputQCA3D(const char* filename="lbout", int iprop=0)
{
  int i, j, k;
  float ipos=0.0;
  long ilen=0;
  char buf[80];
  sprintf(buf, "%s%.4d.q", filename, qVersion);
  ofstream file(buf);
  file<<lbsy.nx<<" "<<lbsy.ny<<" "<<lbsy.nz<<" "<<endl;
  file<<lbsoundv<<" "<<1.0<<" "<<lbreynolds<<" "<<qVersion<<" "<<endl;
  // write mass fraction
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++) 
	{
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos = float(fGetFracSite(iprop, &lbf[ilen*lbsitelength]));
	  file<<ipos<<" ";
	}

  if(!incompress) {
  // write x-component of velocity (v_x)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write y-component of velocity (v_y)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write z-component of velocity (v_z)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(2, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  }
  else {
  // write x-component of velocity (v_x)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write y-component of velocity (v_y)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write z-component of velocity (v_z)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(2, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  }

  // write phase-field/space property
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++) 
	{
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos=(float)lbphi[ilen];
	  file<<ipos;
          if(ilen<(lbsy.nx*lbsy.ny*lbsy.nz-1))
            file<<" ";
	}
  file<<endl;
  file.close();
  return 0;
}




int fsOutputQCA2D(const char* filename="lbout", int iprop=0)
{
  int i, j;
  float ipos=0.0;
  long ilen=0;
  char buf[80];
  sprintf(buf, "%s%.4d.q", filename, qVersion);
  ofstream file(buf);
  file<<lbsy.nx<<" "<<lbsy.ny<<endl;
  file<<lbsoundv<<" "<<1.0<<" "<<lbreynolds<<" "<<qVersion<<" "<<endl;
  // write mass fraction
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++) 
      {
	ilen = (i + lbdm.bwid) * lbdm.youter + j + lbdm.bwid;
	ipos = float(fGetFracSite(iprop, &lbf[ilen*lbsitelength]));
	file<<ipos<<" ";
	}

  if(!incompress) {
  // write x-component of velocity (v_x)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  // write y-component of velocity (v_y)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  }
  else {
  // write x-component of velocity (v_x)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  // write y-component of velocity (v_y)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  }

  // write phase field/space property
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++) 
      {
	ilen = (i + lbdm.bwid) * lbdm.youter + j + lbdm.bwid;
	ipos=(float)lbphi[ilen];
	file<<ipos;
        if(ilen<(lbsy.nx*lbsy.ny-1))
          file<<" ";
      }
  file<<endl;
  file.close();
  return 0;
}




int fsOutputQCB3D(const char* filename="lbout", int iprop=0)
{
  int i, j, k;
  float ipos=0.0;
  long ilen=0;
  char buf[80];
  sprintf(buf, "%s%.4d.q", filename, qVersion);
  ofstream file(buf);
  file<<lbsy.nx<<" "<<lbsy.ny<<" "<<lbsy.nz<<" "<<endl;
  file<<lbsoundv<<" "<<1.0<<" "<<lbreynolds<<" "<<qVersion<<" "<<endl;
  // write solute concentration
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++) 
	{
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
          ipos = float(fGetOneConcSite(iprop, ilen));
	  file<<ipos<<" ";
	}

  if(!incompress) {
  // write x-component of velocity (v_x)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write y-component of velocity (v_y)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write z-component of velocity (v_z)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(2, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  }
  else {
  // write x-component of velocity (v_x)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write y-component of velocity (v_y)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write z-component of velocity (v_z)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(2, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  }

  // write phase field/space property
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++) 
	{
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos=(float)lbphi[ilen];
	  file<<ipos;
          if(ilen<(lbsy.nx*lbsy.ny*lbsy.nz-1))
            file<<" ";
	}
  file<<endl;
  file.close();
  return 0;
}




int fsOutputQCB2D(const char* filename="lbout", int iprop=0)
{
  int i, j;
  float ipos=0.0;
  long ilen=0;
  char buf[80];
  sprintf(buf, "%s%.4d.q", filename, qVersion);
  ofstream file(buf);
  file<<lbsy.nx<<" "<<lbsy.ny<<endl;
  file<<lbsoundv<<" "<<1.0<<" "<<lbreynolds<<" "<<qVersion<<" "<<endl;
  // write solute concentration
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++) 
      {
	ilen = (i + lbdm.bwid) * lbdm.youter + j + lbdm.bwid;
	ipos = float(fGetOneConcSite(iprop, ilen));
	file<<ipos<<" ";
	}

  if(!incompress) {
  // write x-component of velocity (v_x)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  // write y-component of velocity (v_y)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  }
  else {
  // write x-component of velocity (v_x)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  // write y-component of velocity (v_y)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  }

  // write phase field/space property
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++) 
      {
	ilen = (i + lbdm.bwid) * lbdm.youter + j + lbdm.bwid;
	ipos=(float)lbphi[ilen];
	file<<ipos;
        if(ilen<(lbsy.nx*lbsy.ny-1))
          file<<" ";
      }
  file<<endl;
  file.close();
  return 0;
}




int fsOutputQT3D(const char* filename="lbout")
{
  int i, j, k;
  float ipos=0.0;
  long ilen=0;
  char buf[80];
  sprintf(buf, "%s%.4d.q", filename, qVersion);
  ofstream file(buf);
  file<<lbsy.nx<<" "<<lbsy.ny<<" "<<lbsy.nz<<" "<<endl;
  file<<lbsoundv<<" "<<1.0<<" "<<lbreynolds<<" "<<qVersion<<" "<<endl;
  // write temperature
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++) 
	{
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos = float(fGetTemperatureSite(ilen));
	  file<<ipos<<" ";
	}

  if(!incompress) {
  // write x-component of velocity (v_x)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write y-component of velocity (v_y)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write z-component of velocity (v_z)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(2, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  }
  else {
  // write x-component of velocity (v_x)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write y-component of velocity (v_y)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  // write z-component of velocity (v_z)
    for(k=0; k<lbsy.nz; k++) {
      for(j=0; j<lbsy.ny; j++) {
        for(i=0; i<lbsy.nx; i++) {
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	    ipos = 0;
 	  else
	    ipos=fGetOneDirecSpeedSite(2, &lbf[ilen*lbsitelength]);
	  file<<ipos<<" ";
	}
      }
    }
  }

  // write phase field/space property
  for(k=0; k<lbsy.nz; k++)
    for(j=0; j<lbsy.ny; j++)
      for(i=0; i<lbsy.nx; i++) 
	{
	  ilen = ((i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid) * lbdm.zouter + k + lbdm.bwid;
	  ipos=(float)lbphi[ilen];
	  file<<ipos;
          if(ilen<(lbsy.nx*lbsy.ny*lbsy.nz-1))
            file<<" ";
	}
  file<<endl;
  file.close();
  return 0;
}




int fsOutputQT2D(const char* filename="lbout")
{
  int i, j;
  float ipos=0.0;
  long ilen=0;
  char buf[80];
  sprintf(buf, "%s%d.q", filename, qVersion);
  ofstream file(buf);
  file<<lbsy.nx<<" "<<lbsy.ny<<endl;
  file<<lbsoundv<<" "<<1.0<<" "<<lbreynolds<<" "<<qVersion<<" "<<endl;
  // write temperature
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++) 
      {
	ilen = (i + lbdm.bwid) * lbdm.youter + j + lbdm.bwid;
	ipos = float(fGetTemperatureSite(ilen));
	file<<ipos<<" ";
	}

  if(!incompress) {
  // write x-component of velocity (v_x)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  // write y-component of velocity (v_y)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  }
  else {
  // write x-component of velocity (v_x)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(0, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  // write y-component of velocity (v_y)
    for(j=0; j<lbsy.ny; j++) {
      for(i=0; i<lbsy.nx; i++) {
	ilen = (i + lbdm.bwid)* lbdm.youter + j + lbdm.bwid;
	if(lbphi[ilen] == 12 || lbphi[ilen] == 13)
	  ipos = 0;
 	else
	  ipos=fGetOneDirecSpeedSite(1, &lbf[ilen*lbsitelength]);
	file<<ipos<<" ";
      }
    }
  }

  // write phase field/space property
  for(j=0; j<lbsy.ny; j++)
    for(i=0; i<lbsy.nx; i++) 
      {
	ilen = (i + lbdm.bwid) * lbdm.youter + j + lbdm.bwid;
	ipos=(float)lbphi[ilen];
	file<<ipos;
        if(ilen<(lbsy.nx*lbsy.ny-1))
          file<<" ";
      }
  file<<endl;
  file.close();
  return 0;
}



