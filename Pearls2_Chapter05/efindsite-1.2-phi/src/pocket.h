/*
===============================================================================
         ___________.__            .____________.__  __          
     ____\_   _____/|__| ____    __| _/   _____/|__|/  |_  ____  
   _/ __ \|    __)  |  |/    \  / __ |\_____  \ |  \   __\/ __ \ 
   \  ___/|     \   |  |   |  \/ /_/ |/        \|  ||  | \  ___/ 
    \___  >___  /   |__|___|  /\____ /_______  /|__||__|  \___  >
        \/    \/            \/      \/       \/               \/ 

                                                  
   eFindSite - ligand-binding site prediction from meta-threading

   Computational Systems Biology Group
   Department of Biological Sciences
   Center for Computation & Technology
   Louisiana State University
   407 Choppin Hall, Baton Rouge, LA 70803, USA

   http://www.brylinski.org

   Report bugs to michal@brylinski.org wfeinstein@lsu.edu

   Copyright 2013 Michal Brylinski, Wei Feinstein

   This file is part of eFindSite.

   eFindSite is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   eFindSite is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with eFindSite. If not, see <http://www.gnu.org/licenses/>.

===============================================================================
*/


#ifndef __POCKET_H_
#define __POCKET_H_

#include<list>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<cstring>
#include<bitset>
#include<map>

#include "size.h"
#include "template.h"
#include "target.h"
#include "tanimoto.h"
#include "runsvm.h"
#include "data.h"
#include "cluster.h"
#include "cmps.h"

using namespace std;

struct binding_residue
{
 int         residue_number;
 std::string residue_code;
 double      residue_comp_svm;
 double      residue_bind_svm;
 double      residue_score;
 double      residue_distance;
 double      residue_profile1[20];
 double      residue_profile2[20];
 int         template_incontact;
 double      template_fraction;
 double      template_fraction_bkg;
 double      residue_coords[3];
};

struct ligand_fpt_smi
{
 int                 ftp_number;
 double              fpt_fraction;
 std::string         fpt_pdbid;
 std::bitset<MAXSMI> fpt_fingerprint;
};

struct ligand_fpt_mac
{
 int                 ftp_number;
 double              fpt_fraction;
 std::string         fpt_pdbid;
 std::bitset<MAXMAC> fpt_fingerprint;
};

class Pocket {
        
  private:
    
    list< Template * >       _complexes;                // template complexes
    
    int                      _pocket_number;            // pocket number
    int                      _ligand_tot;               // number of ligands
    int                      _ligclust_tot_smi;         // number of ligand clusters SMILES
    int                      _ligclust_tot_mac;         // number of ligand clusters MACCS
    
    int                      _binres_tot;               // number of predicted binding residues
    
    double                   _pocket_center[4];         // ligand geometric center
    
    double                   _pocket_fraction;          // fraction of template-bound ligands
    
    list<std::string>        _list_prot;                // list of proteins
    list<std::string>        _list_lig;                 // list of ligands
    
    map<int,binding_residue> _binding_res;              // ligand-binding residues
    
    list<ligand_fpt_smi>     _ligand_fpts_smi;          // fingerprints SMILES
    list<ligand_fpt_mac>     _ligand_fpts_mac;          // fingerprints MACCS
    
    double                   _profile_fpts_smi[MAXSMI]; // fingerprint profile
    double                   _profile_fpts_mac[MAXMAC]; // fingerprint profile
    
    double                   _ave_mw;                   // weighted average mw
    double                   _ave_logp;                 // weighted average logp
    double                   _ave_psa;                  // weighted average psa
    double                   _ave_mr;                   // weighted average mr
    double                   _ave_hbd;                  // weighted average hb donors
    double                   _ave_hba;                  // weighted average hb acceptors
    
    double                   _std_mw;                   // weighted stdev mw
    double                   _std_logp;                 // weighted stdev logp
    double                   _std_psa;                  // weighted stdev psa
    double                   _std_mr;                   // weighted stdev mr
    double                   _std_hbd;                  // weighted stdev hb donors
    double                   _std_hba;                  // weighted stdev hb acceptors
    
    double                   _score_cmps[MAXSV4];       // auxiliary cmps scores
    
    double                   _svm_cmps;                 // svm score for auxiliary cmps
    
    double                   _plb;                      // protein-ligand binding index
    
    double                   _ave_bres;                 // average score for binding residues
    
    double                   _svm_confidence;           // svm confidence index
    
  public:
    
    Pocket( int );
    
    Pocket( void );
    
    ~Pocket();
    
    int getProteinsTotal( void );
    
    int getLigandsTotal( void );
    
    void addTemplate( Template * );
    
    void calculatePocketCenter( void );
    
    int calculateBindingResidues( Target *, ModelSVM *, double );
    
    void calculateFingerprintsSMILES( double, std::string );
    
    void calculateFingerprintsMACCS( double, std::string );
    
    void setCenter( double, double );
    
    void dumpProteinsAlignments( std::string, map<string,bool> &, Target * );
    
    void dumpLigands( std::string, map<string,bool> &, int );
    
    void dumpPocket( std::string, Target *, double, int );
    
    void setPocketFraction( double );
    
    void calculateCmpsScores( Cmps *, ModelSVM * );
    
    double calculateConfidence( bool, ModelSVM * );
    
    double getConfidence( void );
};

#endif
