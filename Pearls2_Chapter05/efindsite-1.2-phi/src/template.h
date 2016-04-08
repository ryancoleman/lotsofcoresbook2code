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


#ifndef __TEMPLATE_H_
#define __TEMPLATE_H_

#include<algorithm>
#include<cmath>
#include<cstring>
#include<vector>
#include<bitset>
#include<list>
#include<map>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<gzstream.h>

#include "coords.h"
#include "data.h"
#include "nwalign.h"
#include "size.h"
#include "target.h"

using namespace std;

struct tpl_alignment
{
 double tpl_tms;
 double tpl_rms;
 double tpl_sid;
 int    tpl_len;
 
 double tpl_t[3];
 double tpl_u[3][3];
 
 std::string tpl_aln1;
 std::string tpl_aln2;
};

struct lig_contact
{
 int         residue_number;
 char        residue_name;
 std::string residue_atom;
 int         atom_number;
 std::string atom_name;
 double      contact_distance;
};

struct lig_binding
{
 int    residue_number;
 char   residue_name;
 bool   residue_sidechain;
 double contact_distance;
};

class Template {
        
  private:
    
    vector<CoordsProtein>     _protein_xyz;             // protein coords
    vector<CoordsLigand>      _ligand_xyz[MAXLIG];      // ligand heavy atom coords
    
    double                    _ligand_cen[MAXLIG][3];   // ligand geometric centers
    
    list<string>              _ligand_sdf[MAXLIG];      // ligands in SDF
    
    bitset<MAXSMI>            _ligand_fpt_smi[MAXLIG];  // ligand fingerprint SMILES
    bitset<MAXMAC>            _ligand_fpt_mac[MAXLIG];  // ligand fingerprint MACCS
    
    string                    _protein_id;              // protein PDB-ID
    string                    _ligand_id[MAXLIG];       // ligand PDB-IDs
    
    int                       _na;                      // number of protein atoms
    int                       _nr;                      // number of protein residues
    int                       _nl;                      // number of ligands
    
    int                       _nla[MAXLIG];             // number of ligand heavy atoms
    int                       _nlb[MAXLIG];             // number of ligand bonds
    
    string                    _protein_seq1;            // aa sequence
    char                      _protein_seq2[MAXPRO];    // aa sequence
    int                       _protein_seq3[MAXPRO];    // aa sequence numbering
    
    int                       _protein_align1[MAXPRO];  // alignment to target from NW
    int                       _protein_align2[MAXPRO];  // alignment to target from TMalign
    
    double                    _score_nw;                // sequence identity to target
    double                    _score_tm;                // TM-score to target
    double                    _score_rms;               // RMSD to target
    double                    _score_id;                // sequence identity to target (from structure alignment)
    
    int                       _score_len_tm;            // structure alignment length
    int                       _score_len_nw;            // sequence alignment length
    
    double                    _tm_t[3];                 // TM-align transformation vector to target
    double                    _tm_u[3][3];              // TM-align rotation matrix to target
    
    int                       _pocket_num[MAXLIG];      // pocket number
    int                       _pocket_clu_smi[MAXLIG];  // ligand cluster number SMILES
    int                       _pocket_clu_mac[MAXLIG];  // ligand cluster number MACCS
    
    double                    _ligand_mw[MAXLIG];       // ligand molecular weight
    double                    _ligand_logp[MAXLIG];     // ligand logp
    double                    _ligand_psa[MAXLIG];      // ligand polar surface area
    double                    _ligand_mr[MAXLIG];       // molar refarctivity
    int                       _ligand_hbd[MAXLIG];      // number of HB donors
    int                       _ligand_hba[MAXLIG];      // number of HB acceptors
    
    map<string,tpl_alignment> _aligment_pair;           // template pairwise alignments
    
    double                    _prob_pkt;                // probability for pocket location
    double                    _prob_lig;                // probability for ligand similarity
    
    list<lig_contact>         _contacts_ligand[MAXLIG]; // protein-ligand contacts
    list<lig_binding>         _bindings_ligand[MAXLIG]; // binding residues
    
  public:
    
    Template( int, int, int, double, double );
    
    Template( void );
    
    ~Template();
    
    bool loadTemplate( std::string, list<string> & );
    
    double alignNW( std::string );
    
    int getLigandsTotal( void );
    
    int getProteinResiduesTotal( void );
    
    std::string getProteinSequence( void );
    
    int getProteinAtomsTotal( void );
    
    int getProteinCoordsCA( double [][3] );
    
    void getProteinCoords1D( double [] );
    
    void setPocketNumber( int, int );
    
    int getPocketNumber( int );
    
    void purgeAlignments( list<string> & );
    
    void calculateContacts( void );
    
    void getLigandCenter( int, double [], bool );
    
    std::string getProteinID( void );
    
    std::string getLigandID( int );
    
    void setProteinID( std::string );
    
    void setLigandID( int, std::string );
    
    bool getMatrix( std::string, double [], double [][3] );
    
    void getNWalignment( int [], int );
    
    void getTMalignment( int [], int );
    
    void setTMalignment( int [], int );
    
    void setMatrix( double [], double [][3] );
    
    void getBindingResidues( int, list<lig_binding> & );
    
    double getProbPkt( void );
    
    double getProbLig( void );
    
    double getLigandProp( int, int );
    
    void getLigandFingerprintSMILES ( int, bitset<MAXSMI> & );
    
    void getLigandFingerprintMACCS ( int, bitset<MAXMAC> & );
    
    int getPocketClusterNumberSMILES( int );
    
    void setPocketClusterNumberSMILES( int, int );
    
    int getPocketClusterNumberMACCS( int );
    
    void setPocketClusterNumberMACCS( int, int );
    
    void dumpProtein( std::string, bool );
    
    void dumpAlignment( std::string, int, std::string, double [][3] );
    
    void dumpLigand( std::string, int, bool, int );
    
    double getProteinSeqID1( void );
    
    double getProteinSeqID2( void );
    
    double getProteinTMscore( void );
    
    double getProteinRMSD( void );
    
    int getProteinLengthNW( void );
    
    int getProteinLengthTM( void );
    
    void setProteinSeqID1( double );
    
    void setProteinSeqID2( double );
    
    void setProteinTMscore( double );
    
    void setProteinRMSD( double );
    
    void setProteinLengthNW( int );
    
    void setProteinLengthTM( int );
    
    int getLigandAtomsTotal( int );
    
    int getLigandBondsTotal( int );
};

#endif
