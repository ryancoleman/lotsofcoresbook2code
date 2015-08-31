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


int fDefineSystem(const char* filename = "lbin.sys")
{
  
  // read calculation parameters from system data file (default: lbin.sys)

  int numb=0;
  string line,issue,value;
  ifstream myfile;

  myfile.open(filename, ios::in);

  if(!myfile)
    {
      cout << "error opening system file" << endl;
      exit(1);
    }

  lbrestart = 0;
  incompress = 0;
  collide = 0;
  interact = 0;
  outformat = 0;

  while (!myfile.eof())
    {

      getline(myfile,line);
      issue = fReadString(line, 0);
      value = fReadString(line, 1);

      if(issue.compare(0,15,"space_dimension")==0) { 
	    lbsy.nd = fStringToNumber <int> (value); numb++; }

      else if(issue.compare(0,14,"discrete_speed")==0) { 
	    lbsy.nq = fStringToNumber <int> (value); numb++; }

      else if(issue.compare(0,15,"number_of_fluid")==0) { 
	    lbsy.nf = fStringToNumber <int> (value); numb++; }

      else if(issue.compare(0,16,"number_of_solute")==0) {
	    lbsy.nc = fStringToNumber <int> (value); numb++; }

      else if(issue.compare(0,18,"temperature_scalar")==0) {
	    lbsy.nt = fStringToNumber <int> (value); numb++; }

      else if(issue.compare(0,11,"phase_field")==0) {
	    lbsy.pf = fStringToNumber <int> (value); numb++; }

      else if(issue.compare(0,13,"grid_number_x")==0) {
	    lbsy.nx = fStringToNumber <int> (value); numb++; }

      else if(issue.compare(0,13,"grid_number_y")==0) {
	    lbsy.ny = fStringToNumber <int> (value); numb++; }

      else if(issue.compare(0,13,"grid_number_z")==0) {
	    lbsy.nz = fStringToNumber <int> (value); numb++; }

      else if(issue.compare(0,21,"domain_boundary_width")==0) { 
	    lbdm.bwid = fStringToNumber <int> (value); numb++; }

      else if(issue.compare(0,21,"incompressible_fluids")==0) { 
	    incompress = fStringToNumber <int> (value); }

      else if(issue.compare(0,18,"restart_simulation")==0) { 
	    lbrestart = fStringToNumber <int> (value); }

      else if(issue.compare(0,14,"collision_type")==0) { 
        if(value.compare("BGK")==0) {
          collide = 0;
        }
        else if(value.compare("BGKEDM")==0) {
          collide = 1;
        }
        else if(value.compare("BGKGuo")==0) {
          collide = 2;
        }
        else if(value.compare("MRT")==0) {
          collide = 3;
        }
        else if(value.compare("MRTEDM")==0) {
          collide = 4;
        }
        else if(value.compare("MRTGuo")==0) {
          collide = 5;
        }
        else {
          collide = fStringToNumber <int> (value);
        }
      }

      else if(issue.compare(0,16,"interaction_type")==0) { 
        if(value.compare("ShanChen")==0) {
          interact = 1;
        }
        else if(value.compare("ShanChenWetting")==0) {
          interact = 1;
        }
        else if(value.compare("ShanChenQuadratic")==0) {
          interact = 2;
        }
        else if(value.compare("Lishchuk")==0) {
          interact = 3;
        }
        else {
          interact = fStringToNumber <int> (value);
        }
      }

      else if(issue.compare(0,13,"output_format")==0) {
	    if(value.compare("VTK")==0) {
          outformat = 0;
        }
        else if(value.compare("LegacyVTK")==0) {
          outformat = 1;
        }
        else if(value.compare("Plot3D")==0) {
          outformat = 2;
        }
        else {
          outformat = fStringToNumber <int> (value);
        }
      }

    }

  if(numb != 10) {
    cout << "error: insufficient system information" << endl;
    exit(1); }
  
  if((lbsy.nd == 2) && (lbsy.nz != 1)) {
    lbsy.nz = 1;
    cout << "warning: lbsy.nz is reset to 1" << endl; }
  
  lbsitelength = (lbsy.nf + lbsy.nc + lbsy.nt) * lbsy.nq + lbsy.pf;

  myfile.close();
  
  return 0;
  
}


int fInputParameters(const char* filename="lbin.sys")
{
  
  // reads system parameters from system data file (default: lbin.sys)
  
  string line,issue,value,str1;
  stringstream sstm;
  ifstream myfile;
  
  myfile.open(filename, ios::in);
  
  if(!myfile)
    {
      cout << "error: cannot open system file" << endl;
      exit(1);
    }

  // default values for equilibration steps and evaporation limit

  lbequstep = 0;
  lbevaplim = 1e-8;
  lbdump = 10000;

  while (!myfile.eof())
    {
      getline(myfile,line);
      issue = fReadString(line, 0);
      value = fReadString(line, 1);
      if(issue.compare(0,10,"total_step")==0)
 	    lbtotstep = fStringToNumber <int> (value);
      else if(issue.compare(0,18,"equilibration_step")==0)
        lbequstep = fStringToNumber <int> (value);
      else if(issue.compare(0,9,"save_span")==0)
	    lbsave = fStringToNumber <int> (value);
      else if(issue.compare(0,9,"dump_span")==0)
	    lbdump = fStringToNumber <int> (value);
      else if(issue.compare(0,16,"calculation_time")==0)
        lbcalctime = fStringToNumber <double> (value);
      else if(issue.compare(0,16,"closedown_time")==0)
        lbendtime = fStringToNumber <double> (value);
      else if(issue.compare(0,15,"noise_intensity")==0)
	    lbnoise = fStringToNumber <double> (value);
      else if(issue.compare(0,11,"sound_speed")==0)
	    lbsoundv = fStringToNumber <double> (value);
      else if(issue.compare(0,17,"kinetic_viscosity")==0)
	    lbkinetic = fStringToNumber <double> (value);
      else if(issue.compare(0,17,"evaporation_limit")==0)
	    lbevaplim = fStringToNumber <double> (value);

      
      else if(issue.compare(0,15,"temperature_ini")==0)
	    lbinit = fStringToNumber <double> (value);
      else if(issue.compare(0,15,"temperature_top")==0)
	    lbtopt = fStringToNumber <double> (value);
      else if(issue.compare(0,18,"temperature_bottom")==0)
	    lbbott = fStringToNumber <double> (value);
      else if(issue.compare(0,17,"temperature_front")==0)
	    lbfrot = fStringToNumber <double> (value);
      else if(issue.compare(0,16,"temperature_back")==0)
	    lbbact = fStringToNumber <double> (value);
      else if(issue.compare(0,16,"temperature_left")==0)
	    lbleft = fStringToNumber <double> (value);
      else if(issue.compare(0,17,"temperature_right")==0)
	    lbrigt = fStringToNumber <double> (value);
      else if(issue.compare(0,27,"temperature_boussinesq_high")==0)
	    lbbousth = fStringToNumber <double> (value);
      else if(issue.compare(0,26,"temperature_boussinesq_low")==0)
	    lbboustl = fStringToNumber <double> (value);
      
      else if(issue.compare(0,16,"heating_rate_sys")==0)
	    lbsysdt = fStringToNumber <double> (value);
      else if(issue.compare(0,16,"heating_rate_top")==0)
	    lbtopdt = fStringToNumber <double> (value);
      else if(issue.compare(0,19,"heating_rate_bottom")==0)
	    lbbotdt = fStringToNumber <double> (value);
      else if(issue.compare(0,18,"heating_rate_front")==0)
	    lbfrodt = fStringToNumber <double> (value);
      else if(issue.compare(0,17,"heating_rate_back")==0)
	    lbbacdt = fStringToNumber <double> (value);
      else if(issue.compare(0,17,"heating_rate_left")==0)
	    lblefdt = fStringToNumber <double> (value);
      else if(issue.compare(0,18,"heating_rate_right")==0)
	    lbrigdt = fStringToNumber <double> (value);

      else if(issue.compare(0,11,"segregation")==0) {
        for(int i=0; i<lbsy.nf*lbsy.nf; i++) {
	      lbseg[i]=fStringToNumber <double> (value);
        }
      }

	  else if(issue.compare(0,14,"potential_type")==0) {
        if(value.compare("ShanChen1993")==0) {
          for(int i=0; i<lbsy.nf; i++) {
            lbscpot[i] = 0;
          }
        }
        else if(value.compare("ShanChen1994")==0) {
          for(int i=0; i<lbsy.nf; i++) {
            lbscpot[i] = 1;
          }
        }
        else if(value.compare("Rho")==0) {
          for(int i=0; i<lbsy.nf; i++) {
            lbscpot[i] = 2;
          }
        }
        else {
          for(int i=0; i<lbsy.nf; i++) {
            lbscpot[i] = fStringToNumber <double> (value);
          }
        }
	  }

      else {
	    for(int i=0; i<3; i++) {
          sstm.str(string());
          sstm << "speed_ini_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbiniv[i] = fStringToNumber <double> (value); goto Found;}
	  
          sstm.str(string());
          sstm << "speed_top_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbtopv[i] = fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "speed_bot_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbbotv[i] = fStringToNumber <double> (value); goto Found;}
	  
          sstm.str(string());
          sstm << "speed_fro_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbfrov[i] = fStringToNumber <double> (value); goto Found;}
	  
          sstm.str(string());
          sstm << "speed_bac_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbbacv[i] = fStringToNumber <double> (value); goto Found;}
	  
          sstm.str(string());
          sstm << "speed_lef_" << i;
	      str1 = sstm.str();
    	  if(issue.compare(str1)==0) {
	        lblefv[i] = fStringToNumber <double> (value); goto Found;}
	  
          sstm.str(string());
          sstm << "speed_rig_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbrigv[i] = fStringToNumber <double> (value); goto Found;}
	  
	}
	for(int i=0; i<lbsy.nf; i++) {

/*
          sstm.str(string());
          sstm << "embryo_radius_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbemradius[i] = fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "embryo_number_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbemnumber[i] = fStringToNumber <int> (value); goto Found;}
*/

          sstm.str(string());
          sstm << "density_inc_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbincp[i] = fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "density_ini_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbinip[i] = fStringToNumber <double> (value);
            if(!incompress && lbincp[i]==0.0)
              lbincp[i] = fStringToNumber <double> (value);
            goto Found;}
	  
          sstm.str(string());
          sstm << "density_top_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbtopp[i] = fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "density_bot_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbbotp[i] = fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "density_fro_" << i;
  	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbfrop[i] = fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "density_bac_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbbacp[i] = fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "density_lef_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lblefp[i] = fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "density_rig_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbrigp[i] = fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "relaxation_fluid_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbtf[i] = 1.0 / fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "relax_freq_fluid_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbtf[i] = fStringToNumber <double> (value);  goto Found;}

          sstm.str(string());
          sstm << "bulk_relaxation_fluid_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbtfbulk[i] = 1.0 / fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "bulk_relax_freq_fluid_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbtfbulk[i] = fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "shanchen_psi0_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbpsi0[i] = fStringToNumber <double> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "potential_type_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
            if(value.compare("ShanChen1993")==0) {
              lbscpot[i] = 0;
            }
            else if(value.compare("ShanChen1994")==0) {
              lbscpot[i] = 1;
            }
            else if(value.compare("Rho")==0) {
              lbscpot[i] = 2;
            }
            else {
	          lbscpot[i] = fStringToNumber <double> (value);
            }
            goto Found;}
	  
	      for(int j=0; j<lbsy.nf; j++) {
            sstm.str(string());
            sstm << "interaction_" << (i*lbsy.nf+j);
	        str1 = sstm.str();
	        if(issue.compare(str1)==0) {
	          lbg[i*lbsy.nf+j]=fStringToNumber <double> (value);  goto Found;}

            sstm.str(string());
            sstm << "interaction_" << i << "_" << j;
	        str1 = sstm.str();
	        if(issue.compare(str1)==0) {
	          lbg[i*lbsy.nf+j]=fStringToNumber <double> (value);
	          lbg[j*lbsy.nf+i]=fStringToNumber <double> (value);  goto Found;}

            sstm.str(string());
            sstm << "segregation_" << (i*lbsy.nf+j);
	        str1 = sstm.str();
	        if(issue.compare(str1)==0) {
	          lbseg[i*lbsy.nf+j]=fStringToNumber <double> (value);  goto Found;}

            sstm.str(string());
            sstm << "segregation_" << i << "_" << j;
	        str1 = sstm.str();
	        if(issue.compare(str1)==0) {
	          lbseg[i*lbsy.nf+j]=fStringToNumber <double> (value);
	          lbseg[j*lbsy.nf+i]=fStringToNumber <double> (value);  goto Found;}
	  }

          sstm.str(string());
          sstm << "wall_interaction_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbgwall[i] = fStringToNumber <double> (value);  goto Found;}


          for(int j=0; j<3; j++) {
            sstm.str(string());
            sstm << "body_force_" << (i*3+j);
	        str1 = sstm.str();
	        if(issue.compare(str1)==0) {
	          lbbdforce[i*3+j] = fStringToNumber <double> (value);  goto Found;}

            sstm.str(string());
            sstm << "boussinesq_force_" << (i*3+j);
	        str1 = sstm.str();
	        if(issue.compare(str1)==0) {
	          lbbousforce[i*3+j] = fStringToNumber <double> (value);  goto Found;}
	  }

          sstm.str(string());
          sstm << "body_force_x_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbbdforce[i*3] = fStringToNumber <double> (value);  goto Found;}
          sstm.str(string());
          sstm << "body_force_y_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbbdforce[i*3+1] = fStringToNumber <double> (value);  goto Found;}
          sstm.str(string());
          sstm << "body_force_z_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbbdforce[i*3+2] = fStringToNumber <double> (value);  goto Found;}
          sstm.str(string());
          sstm << "boussinesq_force_x_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbbousforce[i*3] = fStringToNumber <double> (value);  goto Found;}
          sstm.str(string());
          sstm << "boussinesq_force_y_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbbousforce[i*3+1] = fStringToNumber <double> (value);  goto Found;}
          sstm.str(string());
          sstm << "boussinesq_force_z_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbbousforce[i*3+2] = fStringToNumber <double> (value);  goto Found;}
	  
	}

	for(int j=0; j<lbsy.nc; j++) {
          sstm.str(string());
          sstm << "relax_solute_" << j;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbtc[j]= 1.0 / fStringToNumber <double> (value);  goto Found;}
          sstm.str(string());
          sstm << "relax_freq_solute_" << j;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbtc[j]= fStringToNumber <double> (value);  goto Found;}
	}
	  
	for(int j=0; j<lbsy.nt; j++) {
	  if(issue.compare(0,13,"relax_thermal")==0) {
	    lbtt[j]= 1.0 / fStringToNumber <double> (value); goto Found;}
	  if(issue.compare(0,18,"relax_freq_thermal")==0) {
	    lbtt[j]= fStringToNumber <double> (value); goto Found;}
	}

/*
	for(int i=0; i<(lbsy.nf * (lbsy.nf-1))/2; i++) {
          sstm.str(string());
          sstm << "interface_fold_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbanifold[i] = fStringToNumber <int> (value);  goto Found;}
	  
          sstm.str(string());
          sstm << "anisotropy_intensity_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbanitens[i] = fStringToNumber <double> (value);  goto Found;}
	}
*/	
	
	for(int i=0; i<lbsy.nc; i++) {
          sstm.str(string());
          sstm << "solute_ini_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbinic[i] = fStringToNumber <double> (value); goto Found;}
	  
          sstm.str(string());
          sstm << "solute_top_" << i;
  	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbtopc[i] = fStringToNumber <double> (value); goto Found;}
	  
          sstm.str(string());
          sstm << "solute_bot_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbbotc[i] = fStringToNumber <double> (value); goto Found;}
	  
          sstm.str(string());
          sstm << "solute_fro_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbfroc[i] = fStringToNumber <double> (value); goto Found;}
	  
          sstm.str(string());
          sstm << "solute_bac_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbbacc[i] = fStringToNumber <double> (value); goto Found;}
	  
          sstm.str(string());
          sstm << "solute_lef_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lblefc[i] = fStringToNumber <double> (value); goto Found;}
	  
          sstm.str(string());
          sstm << "solute_rig_" << i;
	      str1 = sstm.str();
	      if(issue.compare(str1)==0) {
	        lbrigc[i] = fStringToNumber <double> (value); goto Found;}
	    }
      }
    Found: ;
    }

  lbdt = lbkinetic / (lbsoundv * lbsoundv * (1.0/lbtf[0] -0.5));

  if(lbdt <= 0) {
    cout << "The given relaxation time is not suitable for LBE calculation!" << endl;
    cout << "error: program terminated" << endl;
    exit(1);
  }

  postequil = (lbequstep==0);

  myfile.close();
  
  return 0;
}

int fReadSpace2D(const char* filename="lbin.spa")
{

  // read 2D space parameters from data file (default: lbin.spa)
   
  int xpos, ypos, zpos, value;
  unsigned long rpos;

  ifstream myfile;
  myfile.open(filename, ios::in);
  if(!myfile)
    {
      cout << "error: cannot open input file" << endl;
      exit(1);
    }
  while (!myfile.eof())
    {
      myfile >> xpos >> ypos >> zpos >> value;
      if(!(postequil) && value>110 && value<890)
        value = 0;
      if((xpos >= lbdm.xs) && (xpos < lbdm.xe) && 
	 (ypos >= lbdm.ys) && (ypos < lbdm.ye) &&
	 (zpos == 0))
	{
	  rpos = (xpos - lbdm.xs + lbdm.bwid) * lbdm.youter 
	    + ypos - lbdm.ys + lbdm.bwid; 
	  lbphi[rpos] = value;
	}
    }
  myfile.close();
  return 0;
}


int fReadSpace3D(const char* filename="lbin.spa")
{
  
  // read 3D space parameters from data file (default: lbin.spa)
  
  int xpos, ypos, zpos, value;
  unsigned long rpos;
  
  ifstream myfile;
  myfile.open(filename, ios::in);
  if(!myfile)
    {
      cout << "Error opening input file" << endl;
      exit(1);
    }
  while (!myfile.eof())
    {
      myfile >> xpos >> ypos >> zpos >> value;
      if((xpos >= lbdm.xs) && (xpos <lbdm.xe) && 
	 (ypos >= lbdm.ys) && (ypos <lbdm.ye) &&
	 (zpos >= lbdm.zs) && (zpos <lbdm.ze))
	{
	  rpos = ((xpos - lbdm.xs + lbdm.bwid) * lbdm.youter 
		  + (ypos - lbdm.ys + lbdm.bwid)) * lbdm.zouter
	    + zpos - lbdm.zs + lbdm.bwid; 
	  lbphi[rpos] = value;
	}
    }
  myfile.close();
  return 0;
}


int fReadSpaceParameter(const char* filename="lbin.spa")
{

  if(lbsy.nd == 2)
    fReadSpace2D(filename);
  else
    fReadSpace3D(filename);
  return 0;
}


int fReadInitialState2D(const char* filename="lbin.init")
{
  
  // read initial velocities, densities, concentrations, temperature
  // from data file to replace default values (default file: lbin.init)

  int xpos, ypos, zpos;
  unsigned long rpos;
  double vel[3];
  double dens[lbsy.nf];
  double conc[lbsy.nc];
  double temp;
  double * pt2;

  ifstream myfile;
  myfile.open(filename, ios::in);
  if(myfile)
    {
    while (!myfile.eof())
      {
        myfile >> xpos >> ypos >> zpos;
        myfile >> vel[0] >> vel[1] >> vel[2];
        for (int m=0; m<lbsy.nf; m++) {
          myfile >> dens[m];
        }
        if (lbsy.nc>0) {
          for (int m=0; m<lbsy.nc; m++) {
            myfile >> conc[m];
          }
        }
        if (lbsy.nt==1) {
          myfile >> temp;
        }
  //  replace populations for points inside domain

        if((xpos >= lbdm.xs) && (xpos <lbdm.xe) && 
	   (ypos >= lbdm.ys) && (ypos <lbdm.ye) &&
           (zpos == 0))
  	  {
	    rpos = (xpos - lbdm.xs + lbdm.bwid) * lbdm.youter 
		  + (ypos - lbdm.ys + lbdm.bwid); 
	    pt2 = &lbf[lbsitelength*rpos];
            if(!incompress) {
              for (int m=0; m<lbsy.nf; m++) {
                fGetEquilibriumF(lbfeq, vel, dens[m]);
                for (int n=0; n<lbsy.nq; n++) {
                  *pt2 = lbfeq[n];
                  pt2++;
                }
              }
            }
            else {
              for (int m=0; m<lbsy.nf; m++) {
                fGetEquilibriumFIncom(lbfeq, vel, dens[m], lbinip[m]);
                for (int n=0; n<lbsy.nq; n++) {
                  *pt2 = lbfeq[n];
                  pt2++;
                }
              }
            }
            for (int m=0; m<lbsy.nc; m++) {
              fGetEquilibriumC(lbfeq, vel, conc[m]);
              for (int n=0; n<lbsy.nq; n++) {
                *pt2 = lbfeq[n];
                pt2++;
              }
            }
            if (lbsy.nt==1) {
              fGetEquilibriumT(lbfeq, vel, temp);
              for (int n=0; n<lbsy.nq; n++) {
                *pt2 = lbfeq[n];
                pt2++;
              }            
            }
  	  }

      }
    }
  myfile.close();
  pt2=NULL;

  return 0;
}


int fReadInitialState3D(const char* filename="lbin.init")
{
  
  // read initial velocities, densities, concentrations, temperature
  // from data file to replace default values (default file: lbin.init)

  int xpos, ypos, zpos;
  unsigned long rpos;
  double vel[3];
  double dens[lbsy.nf];
  double conc[lbsy.nc];
  double temp;
  double * pt2;

  ifstream myfile;
  myfile.open(filename, ios::in);
  if(myfile)
    {
    while (!myfile.eof())
      {
        myfile >> xpos >> ypos >> zpos;
        myfile >> vel[0] >> vel[1] >> vel[2];
        for (int m=0; m<lbsy.nf; m++) {
          myfile >> dens[m];
        }
        if (lbsy.nc>0) {
          for (int m=0; m<lbsy.nc; m++) {
            myfile >> conc[m];
          }
        }
        if (lbsy.nt==1) {
          myfile >> temp;
        }
  //  replace populations for points inside domain

        if((xpos >= lbdm.xs) && (xpos <lbdm.xe) && 
	   (ypos >= lbdm.ys) && (ypos <lbdm.ye) &&
  	   (zpos >= lbdm.zs) && (zpos <lbdm.ze))
  	  {
	    rpos = ((xpos - lbdm.xs + lbdm.bwid) * lbdm.youter 
		  + (ypos - lbdm.ys + lbdm.bwid)) * lbdm.zouter
	          + zpos - lbdm.zs + lbdm.bwid; 
	    pt2 = &lbf[lbsitelength*rpos];
            if(!incompress) {
              for (int m=0; m<lbsy.nf; m++) {
                fGetEquilibriumF(lbfeq, vel, dens[m]);
                for (int n=0; n<lbsy.nq; n++) {
                  *pt2 = lbfeq[n];
                  pt2++;
                }
              }
            }
            else {
              for (int m=0; m<lbsy.nf; m++) {
                fGetEquilibriumFIncom(lbfeq, vel, dens[m], lbinip[m]);
                for (int n=0; n<lbsy.nq; n++) {
                  *pt2 = lbfeq[n];
                  pt2++;
                }
              }
            }
            for (int m=0; m<lbsy.nc; m++) {
              fGetEquilibriumC(lbfeq, vel, conc[m]);
              for (int n=0; n<lbsy.nq; n++) {
                *pt2 = lbfeq[n];
                pt2++;
              }
            }
            if (lbsy.nt==1) {
              fGetEquilibriumT(lbfeq, vel, temp);
              for (int n=0; n<lbsy.nq; n++) {
                *pt2 = lbfeq[n];
                pt2++;
              }            
            }
  	  }

      }
    }
  myfile.close();
  pt2=NULL;

  return 0;
}

int fReadInitialState(const char* filename="lbin.init")
{

  if(lbsy.nd == 2)
    fReadInitialState2D(filename);
  else
    fReadInitialState3D(filename);
  return 0;
}






int fPrintSystemInfo()
{

  // print system information

  string buff, bufc, buft, bufcoll, bufforc, bufmeso, bufwet;
//  char buff[30], bufc[20], buft[30], bufcoll[20], bufforc[20], bufmeso[20], bufwet[20];
  int coll=collide/3;
  int forc=collide%3;

  if(!incompress)
    buff.assign((lbsy.nf>1)?"compressible phases":"compressible phase");
  else
    buff.assign((lbsy.nf>1)?"incompressible phases":"incompressible phase");

  bufc.assign((lbsy.nc>1)?"solutes":"solute");
  buft.assign((lbsy.nt==1)?"temperature scalar":"no temperature scalar");
  switch (coll) {
    case 0:
      bufcoll.assign("BGK");
      break;
    case 1:
      bufcoll.assign("MRT");
      break;
  }
  switch (forc) {
    case 0:
      bufforc.assign("standard");
      break;
    case 1:
      bufforc.assign("EDM");
      break;
    case 2:
      bufforc.assign("Guo");
      break;
  }
  switch (interact) {
    case 0:
      bufmeso.assign("no");
      bufwet.assign(" ");
      break;
    case 1:
      bufmeso.assign("Shan-Chen");
      bufwet.assign(" ");
      break;
    case 2:
      bufmeso.assign("Shan-Chen");
      bufwet.assign("with wetting");
      break;
    case 3:
      bufmeso.assign("Shan-Chen");
      bufwet.assign(" ");
      break;
    case 4:
      bufmeso.assign("Shan-Chen");
      bufwet.assign("with wetting");
      break;
    case 5:
      bufmeso.assign("Lishchuk");
      bufwet.assign(" ");
      break;
  }
  if(lbdm.rank == 0) {
    fPrintDoubleLine();
    cout << lbsy.nd << "-dimensional system with grid size: nx=" << lbsy.nx << ", ny=" << lbsy.ny;
    if(lbsy.nd == 3)
      cout << ", nz=" << lbsy.nz; 
    cout << endl;
    cout << "boundary width = " << lbdm.bwid << endl;
    cout << "system includes " << lbsy.nf << " " << buff << ", " << lbsy.nc 
         << " " << bufc << " and " << buft << endl;
    cout << "D" << lbsy.nd << "Q" << lbsy.nq << " lattice Boltzmann model is defined using " 
         << bufcoll << " collisions" << endl;
    cout << "with " << bufforc << " forcing and " << bufmeso << " interactions " << bufwet << endl;
    fPrintLine();
  }
  
  return 0;
}


int fPrintEndEquilibration()
{
  if(lbdm.rank == 0) {
    fPrintDoubleLine();
    cout << "EQUILIBRATION PERIOD NOW COMPLETE" << endl;
    fPrintDoubleLine();   
  }
  return 0;
}


int fPrintDomainMass()
{

  // display masses of individual and all fluids in domain

  cout<<"MASS: "<<"total = "<<fGetTotMassDomain();
  for(int i=0; i<lbsy.nf; i++)
    cout<<", fluid "<<i+1<<" = "<<fGetOneMassDomain(i);
  cout<<endl;
  return 0;
}


int fPrintDomainMomentum()
{

  // display total fluid momentum in domain

  double mome[3];
  fGetTotMomentDomain(mome);
  cout<<"MOMENTUM: "<<"x: "<<mome[0]<<", y: "<<mome[1];
  if(lbsy.nd == 3)
    cout<<", z: "<<mome[2];
  cout<<endl;
  return 0;  
}


int fsOutput(const char* filename="lbout")
{

  // output all properties at current timestep according to required format

  switch (outformat) {
    case 0:
      fsOutputVTK(filename);
      break;
    case 1:
      fsOutputLegacyVTK(filename);
      break;
    case 2:
      fsOutputQ(filename);
      break;
  }
  qVersion++;
  return 0;
}


