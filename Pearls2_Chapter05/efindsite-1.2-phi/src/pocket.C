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


#include "pocket.h"

using namespace std;

Pocket::Pocket( int an )
{
 _pocket_number = an;
 
 _pocket_fraction = 0.0;
 
 _ligand_tot = 0;
 _ligclust_tot_smi = 0;
 _ligclust_tot_mac = 0;
 _binres_tot = 0;
 
 _ave_mw = 0.0;
 _ave_logp = 0.0;
 _ave_psa = 0.0;
 _ave_mr = 0.0;
 _ave_hba = 0.0;
 _ave_hbd = 0.0;
 
 _std_mw = 0.0;
 _std_logp = 0.0;
 _std_psa = 0.0;
 _std_mr = 0.0;
 _std_hba = 0.0;
 _std_hbd = 0.0;
 
 _svm_cmps = 0.0;
 
 _plb = 0.0;
 
 _ave_bres = 0.0;
 
 _svm_confidence = 0.0;
 
 _complexes.clear();
 _binding_res.clear();
 
 for ( int i1 = 0; i1 < 4; i1++ )
  _pocket_center[i1] = 0.0;
 
 for ( int i1 = 0; i1 < MAXSMI; i1++ )
  _profile_fpts_smi[i1] = 0.0;
  
 for ( int i1 = 0; i1 < MAXMAC; i1++ )
  _profile_fpts_mac[i1] = 0.0;
 
 for ( int i1 = 0; i1 < MAXSV4; i1++ )
  _score_cmps[i1] = 0.0;
}

Pocket::Pocket( void )
{
 _pocket_number = 0;
 
 _pocket_fraction = 0.0;
 
 _ligand_tot = 0;
 _ligclust_tot_smi = 0;
 _ligclust_tot_mac = 0;
 _binres_tot = 0;
 
 _ave_mw = 0.0;
 _ave_logp = 0.0;
 _ave_psa = 0.0;
 _ave_mr = 0.0;
 _ave_hba = 0.0;
 _ave_hbd = 0.0;
 
 _std_mw = 0.0;
 _std_logp = 0.0;
 _std_psa = 0.0;
 _std_mr = 0.0;
 _std_hba = 0.0;
 _std_hbd = 0.0;
 
 _svm_cmps = 0.0;
 
 _plb = 0.0;
 
 _ave_bres = 0.0;
 
 _svm_confidence = 0.0;
 
 _complexes.clear();
 _binding_res.clear();
 
 for ( int i1 = 0; i1 < 4; i1++ )
  _pocket_center[i1] = 0.0;
 
 for ( int i1 = 0; i1 < MAXSMI; i1++ )
  _profile_fpts_smi[i1] = 0.0;
  
 for ( int i1 = 0; i1 < MAXMAC; i1++ )
  _profile_fpts_mac[i1] = 0.0;
 
 for ( int i1 = 0; i1 < MAXSV4; i1++ )
  _score_cmps[i1] = 0.0;
}

Pocket::~Pocket() {}


// ==================================================================================   getProteinsTotal

int Pocket::getProteinsTotal( void )
{
 return _list_prot.size();
}


// ==================================================================================   getLigandsTotal

int Pocket::getLigandsTotal( void )
{
 return _list_lig.size();
}


// ==================================================================================   addTemplate

void Pocket::addTemplate( Template * template_new1 )
{
 _complexes.push_back( template_new1 );
 
 _list_prot.clear();
 _list_lig.clear();
 
 for ( list< Template * >::iterator tp1 = _complexes.begin(); tp1 != _complexes.end(); tp1++ )
 {
  _list_prot.push_back( (*tp1)->getProteinID() );
  
  for ( int tl1 = 0; tl1 < (*tp1)->getLigandsTotal(); tl1++ )
   if ( (*tp1)->getPocketNumber(tl1) == _pocket_number ) 
    _list_lig.push_back( (*tp1)->getLigandID(tl1) );
 }
 
 _list_prot.sort();
 _list_prot.unique();
 
 _list_lig.sort();
 _list_lig.unique();
}


// ==================================================================================   calculatePocketCenter

void Pocket::calculatePocketCenter( void )
{
 list< Template * >::iterator ipt1;
 
 for ( int ipt3 = 0; ipt3 < 4; ipt3++ )
  _pocket_center[ipt3] = 0.0;
 
 _ligand_tot = 0;
 
 double wg1 = 0.0;
 double wg2 = 0.0;
 
 for ( ipt1 = _complexes.begin(); ipt1 != _complexes.end(); ipt1++ )
  for ( int ipt2 = 0; ipt2 < (*ipt1)->getLigandsTotal(); ipt2++ )
   if ( (*ipt1)->getPocketNumber(ipt2) == _pocket_number )
   {
    wg1 += (*ipt1)->getProbPkt();
    
    wg2 += (*ipt1)->getProbLig();
   }
 
 for ( ipt1 = _complexes.begin(); ipt1 != _complexes.end(); ipt1++ )
 {
  double cpt1[3];
  
  for ( int ipt2 = 0; ipt2 < (*ipt1)->getLigandsTotal(); ipt2++ )
   if ( (*ipt1)->getPocketNumber(ipt2) == _pocket_number )
   {
    (*ipt1)->getLigandCenter(ipt2, cpt1, true);
    
    for ( int ipt3 = 0; ipt3 < 3; ipt3++ )
     _pocket_center[ipt3] += ( (*ipt1)->getProbPkt() / wg1 ) * cpt1[ipt3];
    
    _ave_mw += ( (*ipt1)->getProbLig() / wg2 ) * (*ipt1)->getLigandProp( ipt2, 1 );
    _ave_logp += ( (*ipt1)->getProbLig() / wg2 ) * (*ipt1)->getLigandProp( ipt2, 2 );
    _ave_psa += ( (*ipt1)->getProbLig() / wg2 ) * (*ipt1)->getLigandProp( ipt2, 3 );
    _ave_mr += ( (*ipt1)->getProbLig() / wg2 ) * (*ipt1)->getLigandProp( ipt2, 4 );
    _ave_hbd += ( (*ipt1)->getProbLig() / wg2 ) * (*ipt1)->getLigandProp( ipt2, 5 );
    _ave_hba += ( (*ipt1)->getProbLig() / wg2 ) * (*ipt1)->getLigandProp( ipt2, 6 );
    
    _ligand_tot++;
   }
 }
 
 for ( ipt1 = _complexes.begin(); ipt1 != _complexes.end(); ipt1++ )
  for ( int ipt2 = 0; ipt2 < (*ipt1)->getLigandsTotal(); ipt2++ )
   if ( (*ipt1)->getPocketNumber(ipt2) == _pocket_number )
   {
    _std_mw += (*ipt1)->getProbLig() * pow( (*ipt1)->getLigandProp( ipt2, 1 ) - _ave_mw, 2.0 );
    _std_logp += (*ipt1)->getProbLig() * pow( (*ipt1)->getLigandProp( ipt2, 2 ) - _ave_logp, 2.0 );
    _std_psa += (*ipt1)->getProbLig() * pow( (*ipt1)->getLigandProp( ipt2, 3 ) - _ave_psa, 2.0 );
    _std_mr += (*ipt1)->getProbLig() * pow( (*ipt1)->getLigandProp( ipt2, 4 ) - _ave_mr, 2.0 );
    _std_hbd += (*ipt1)->getProbLig() * pow( (*ipt1)->getLigandProp( ipt2, 5 ) - _ave_hbd, 2.0 );
    _std_hba += (*ipt1)->getProbLig() * pow( (*ipt1)->getLigandProp( ipt2, 6 ) - _ave_hba, 2.0 );
   }
 
 _std_mw = sqrt( _std_mw / wg2 );
 _std_logp = sqrt( _std_logp / wg2 );
 _std_psa = sqrt( _std_psa / wg2 );
 _std_mr = sqrt( _std_mr / wg2 );
 _std_hbd = sqrt( _std_hbd / wg2 );
 _std_hba = sqrt( _std_hba / wg2 );
 
 for ( ipt1 = _complexes.begin(); ipt1 != _complexes.end(); ipt1++ )
 {
  double cpt1[3];
  
  for ( int ipt2 = 0; ipt2 < (*ipt1)->getLigandsTotal(); ipt2++ )
   if ( (*ipt1)->getPocketNumber(ipt2) == _pocket_number )
   {
    (*ipt1)->getLigandCenter(ipt2, cpt1, true);
    
    _pocket_center[3] += ( (*ipt1)->getProbPkt() / wg1 ) * sqrt(pow(cpt1[0]-_pocket_center[0],2.0) + pow(cpt1[1]-_pocket_center[1],2.0) + pow(cpt1[2]-_pocket_center[2],2.0));
   }
 }
}


// ==================================================================================   calculateBindingResidues

int Pocket::calculateBindingResidues( Target * target1, ModelSVM * model1, double thres1 )
{
 _binding_res.clear();
 
 string res4 = target1->getProteinSequence();
 
 list< Template * >::iterator ipt1;
 
 map<int,binding_residue>::iterator ipt2;
 
 int res1 = target1->getProteinResiduesTotal();
 
 for ( ipt1 = _complexes.begin(); ipt1 != _complexes.end(); ipt1++ )
 {
  int aln1[MAXPRO];
  
  (*ipt1)->getTMalignment(aln1, res1);
  
  for ( int ipt2 = 0; ipt2 < res1; ipt2++ )
   if ( aln1[ipt2] > -1 )
    for ( int ilt1 = 0; ilt1 < (*ipt1)->getLigandsTotal(); ilt1++ )
     if ( (*ipt1)->getPocketNumber(ilt1) == _pocket_number )
     {
      list<lig_binding> res3;
      
      (*ipt1)->getBindingResidues(ilt1, res3);
      
      list<lig_binding>::iterator irt1;
      
      for ( irt1 = res3.begin(); irt1 != res3.end(); irt1++ )
       if ( aln1[ipt2] == (*irt1).residue_number )
       {
        map<int,binding_residue>::iterator ipt3;
        
        ipt3 = _binding_res.find( ipt2 );
        
        if ( ipt3 != _binding_res.end() )
        {
         ((*ipt3).second).template_incontact++;
         
         ((*ipt3).second).residue_profile1[one2num((*irt1).residue_name)] += 1.0;
         ((*ipt3).second).residue_profile2[one2num((*irt1).residue_name)] += 1.0;
        }
        else
        {
         binding_residue res2;
         
         res2.residue_number = ipt2;
         res2.residue_code = ( res4.substr( ipt2, 1 ) ).c_str();
         res2.residue_comp_svm = target1->getCompositionScoreSVM(ipt2);
         res2.residue_bind_svm = 0.0;
         res2.residue_distance = 0.0;
         res2.residue_score = 0.0;
         
         res2.template_incontact = 1;
         res2.template_fraction = 0.0;
         res2.template_fraction_bkg = 0.0;
         
         for ( int irt2 = 0; irt2 < 20; irt2++ )
         {
          res2.residue_profile1[irt2] = 0.0;
          res2.residue_profile2[irt2] = 0.0;
         }
         
         for ( int irt2 = 0; irt2 < 3; irt2++ )
          res2.residue_coords[irt2] = 0.0;
         
         res2.residue_profile1[one2num((*irt1).residue_name)] = 1.0;
         res2.residue_profile2[one2num((*irt1).residue_name)] = 1.0;
         
         _binding_res.insert( pair<int,binding_residue>(ipt2,res2) );
        }
       }
     }
 }
 
 double target2[MAXPRO][3];
 
 target1->getProteinCoordsCA(target2);
 
 ipt2 = _binding_res.begin();
 
 while ( ipt2 != _binding_res.end() )
 {
  if ( _ligand_tot > 0 && ((*ipt2).second).template_incontact > 0 )
  {
   ((*ipt2).second).template_fraction = ( (double) ( ((*ipt2).second).template_incontact ) ) / ( (double) _ligand_tot );
   
   ((*ipt2).second).template_fraction_bkg = ( ( (double) ( ((*ipt2).second).template_incontact ) ) + getBkgFreq1( res4.substr( ((*ipt2).second).residue_number, 1) ) * sqrt( (double) _ligand_tot ) ) / ( ( (double) _ligand_tot ) + sqrt( (double) _ligand_tot ) );
   
   ((*ipt2).second).residue_distance = sqrt( pow( target2[((*ipt2).second).residue_number][0] - _pocket_center[0], 2 ) + pow( target2[((*ipt2).second).residue_number][1] - _pocket_center[1], 2 ) + pow( target2[((*ipt2).second).residue_number][2] - _pocket_center[2], 2 ) );
   
   for ( int irt2 = 0; irt2 < 3; irt2++ )
    ((*ipt2).second).residue_coords[irt2] = target2[((*ipt2).second).residue_number][irt2];
   
   double psvm4[MAXSV2];
   
   for ( int irt3 = 0; irt3 < MAXSV2; irt3++ )
   {
    ((*ipt2).second).residue_profile1[irt3] /= ( (double) ( ((*ipt2).second).template_incontact ) );
    
    ((*ipt2).second).residue_profile2[irt3] = ( ((*ipt2).second).residue_profile2[irt3] + getBkgFreq1( num2one(irt3) ) * sqrt( (double) ( ((*ipt2).second).template_incontact ) ) ) / ( ( (double) ( ((*ipt2).second).template_incontact ) ) + sqrt( (double) ( ((*ipt2).second).template_incontact ) ) );
    
    psvm4[irt3] = ((*ipt2).second).residue_profile2[irt3];
   }
   
   ((*ipt2).second).residue_bind_svm = model1->SVMpredict( 2, psvm4 );
   
   ++ipt2;
  }
  else
   _binding_res.erase(ipt2++);
 }
 
 _binres_tot = 0;
 
 for ( ipt2 = _binding_res.begin(); ipt2 != _binding_res.end(); ipt2++ )
 {
  double rsvm4[MAXSV3];
  
  rsvm4[0] = ((*ipt2).second).residue_comp_svm;
  rsvm4[1] = ((*ipt2).second).residue_bind_svm;
  rsvm4[2] = 1.0 / ((*ipt2).second).residue_distance;
  rsvm4[3] = ((*ipt2).second).template_fraction_bkg;
  rsvm4[4] = _pocket_center[3];
  rsvm4[5] = _ave_mw;
  
  if ( rsvm4[2] > 1.0 )
   rsvm4[2] = 1.0;
  
  if ( rsvm4[4] < 1.0 )
   rsvm4[4] = 1.0;
  
  rsvm4[4] = 1.0 / rsvm4[4];
  
  rsvm4[5] = log(rsvm4[4]);
  
  ((*ipt2).second).residue_score = model1->SVMpredict( 3, rsvm4 );
  
  if ( ((*ipt2).second).residue_score >= thres1 )
  {
   _plb += one2plb(((*ipt2).second).residue_code);
   
   _ave_bres += ((*ipt2).second).residue_score;
   
   _binres_tot++;
  }
 }
 
 if ( _binres_tot > 0 )
 {
  _plb /= (double) _binres_tot;
  
  _ave_bres /= (double) _binres_tot;
 }
 
 return _binres_tot;
}


// ==================================================================================   calculateFingerprintsSMILES

void Pocket::calculateFingerprintsSMILES( double fthr1, std::string fthr2 )
{
 _ligand_fpts_smi.clear();
 
 list< Template * >::iterator fipt1;
 list< Template * >::iterator fipt2;
 
 int fsim2 = 0;
 
 for ( fipt1 = _complexes.begin(); fipt1 != _complexes.end(); fipt1++ )
  for ( int ilt1 = 0; ilt1 < (*fipt1)->getLigandsTotal(); ilt1++ )
   if ( (*fipt1)->getPocketNumber(ilt1) == _pocket_number )
    fsim2++;
 
 double *fsim1 = new double[fsim2*fsim2];
 
 int ftpl1n = 0;
 
 for ( fipt1 = _complexes.begin(); fipt1 != _complexes.end(); fipt1++ )
  for ( int ilt1 = 0; ilt1 < (*fipt1)->getLigandsTotal(); ilt1++ )
   if ( (*fipt1)->getPocketNumber(ilt1) == _pocket_number )
   {
    std::bitset<MAXSMI> ffpt1;
    
    (*fipt1)->getLigandFingerprintSMILES( ilt1, ffpt1 );
    
    int ftpl2n = 0;
    
    for ( fipt2 = _complexes.begin(); fipt2 != _complexes.end(); fipt2++ )
     for ( int ilt2 = 0; ilt2 < (*fipt2)->getLigandsTotal(); ilt2++ )
      if ( (*fipt2)->getPocketNumber(ilt2) == _pocket_number )
      {
       if ( ftpl1n < ftpl2n )
       {
        std::bitset<MAXSMI> ffpt2;
        
        (*fipt2)->getLigandFingerprintSMILES( ilt2, ffpt2 );
        
        if ( fthr2 == "T" )
         fsim1[ftpl1n*fsim2+ftpl2n] = getTanimoto1024(ffpt1, ffpt2);
        else
         fsim1[ftpl1n*fsim2+ftpl2n] = getTanimotoAve1024(ffpt1, ffpt2);
        
        fsim1[ftpl2n*fsim2+ftpl1n] = fsim1[ftpl1n*fsim2+ftpl2n];
       }
       
       ftpl2n++;
      }
    
    ftpl1n++;
   }
 
 int * fclu1 = new int [fsim2];
 
 int fclu2 = cluster_avelink( fsim1, fclu1, fsim2, fthr1, "max" );
 
 _ligclust_tot_smi = fclu2;
 
 for ( int fclu3 = 0; fclu3 < fclu2; fclu3++ )
 {
  int fclu4 = 0;
  double fbst1 = -1.0;
  std::bitset<MAXSMI> fbst2;
  std::string fbst3;
  
  ftpl1n = 0;
  
  for ( fipt1 = _complexes.begin(); fipt1 != _complexes.end(); fipt1++ )
   for ( int ilt1 = 0; ilt1 < (*fipt1)->getLigandsTotal(); ilt1++ )
    if ( (*fipt1)->getPocketNumber(ilt1) == _pocket_number )
    {
     if ( fclu1[ftpl1n] == fclu3 )
     {
      (*fipt1)->setPocketClusterNumberSMILES( ilt1, fclu3 + 1 );
      
      fclu4++;
      
      double ave1 = 0.0;
      double ave2 = 0.0;
      
      std::bitset<MAXSMI> ffpt1;
      
      (*fipt1)->getLigandFingerprintSMILES( ilt1, ffpt1 );
      
      int ftpl2n = 0;
      
      for ( fipt2 = _complexes.begin(); fipt2 != _complexes.end(); fipt2++ )
       for ( int ilt2 = 0; ilt2 < (*fipt2)->getLigandsTotal(); ilt2++ )
        if ( (*fipt2)->getPocketNumber(ilt2) == _pocket_number )
        {
         if ( fclu1[ftpl2n] == fclu3 )
         {
          if ( ftpl1n == ftpl2n )
           ave1 += 1.0;
          else
           ave1 += fsim1[ftpl1n*fsim2+ftpl2n];
          
          ave2 += 1.0;
         }
         
         ave1 -= 1.0;
         ave2 -= 1.0;
         
         if ( ave2 > 0.1 )
          ave1 /= ave2;
         else
          ave1 = 1.0;
         
         if ( ave1 > fbst1 )
         {
          fbst1 = ave1;
          fbst2 = ffpt1;
          fbst3 = (*fipt1)->getLigandID(ilt1);
         }
         
         ftpl2n++;
        }
     }
     
     ftpl1n++;
    }
  
  ligand_fpt_smi fclu5 = { fclu3 + 1, ( (double) fclu4 / (double) _complexes.size() ), fbst3, fbst2 };
  
  _ligand_fpts_smi.push_back( fclu5 );
 }
 
 delete [] fsim1;
 delete [] fclu1;
 
 for ( fipt1 = _complexes.begin(); fipt1 != _complexes.end(); fipt1++ )
  for ( int ilt1 = 0; ilt1 < (*fipt1)->getLigandsTotal(); ilt1++ )
   if ( (*fipt1)->getPocketNumber(ilt1) == _pocket_number )
   {
    std::bitset<MAXSMI> ffpt1;
    
    (*fipt1)->getLigandFingerprintSMILES( ilt1, ffpt1 );
    
    for ( int fipt3 = 0; fipt3 < MAXSMI; fipt3++ )
     _profile_fpts_smi[fipt3] += ( (double) ffpt1[MAXSMI-1-fipt3] / (double) _ligand_tot );
   }
}


// ==================================================================================   calculateFingerprintsMACCS

void Pocket::calculateFingerprintsMACCS( double fthr1, std::string fthr2 )
{
 _ligand_fpts_mac.clear();
 
 list< Template * >::iterator fipt1;
 list< Template * >::iterator fipt2;
 
 int fsim2 = 0;
 
 for ( fipt1 = _complexes.begin(); fipt1 != _complexes.end(); fipt1++ )
  for ( int ilt1 = 0; ilt1 < (*fipt1)->getLigandsTotal(); ilt1++ )
   if ( (*fipt1)->getPocketNumber(ilt1) == _pocket_number )
    fsim2++;
 
 double *fsim1 = new double[fsim2*fsim2];
 
 int ftpl1n = 0;
 
 for ( fipt1 = _complexes.begin(); fipt1 != _complexes.end(); fipt1++ )
  for ( int ilt1 = 0; ilt1 < (*fipt1)->getLigandsTotal(); ilt1++ )
   if ( (*fipt1)->getPocketNumber(ilt1) == _pocket_number )
   {
    std::bitset<MAXMAC> ffpt1;
    
    (*fipt1)->getLigandFingerprintMACCS( ilt1, ffpt1 );
    
    int ftpl2n = 0;
    
    for ( fipt2 = _complexes.begin(); fipt2 != _complexes.end(); fipt2++ )
     for ( int ilt2 = 0; ilt2 < (*fipt2)->getLigandsTotal(); ilt2++ )
      if ( (*fipt2)->getPocketNumber(ilt2) == _pocket_number )
      {
       if ( ftpl1n < ftpl2n )
       {
        std::bitset<MAXMAC> ffpt2;
        
        (*fipt2)->getLigandFingerprintMACCS( ilt2, ffpt2 );
        
        if ( fthr2 == "T" )
         fsim1[ftpl1n*fsim2+ftpl2n] = getTanimoto166(ffpt1, ffpt2);
        else
         fsim1[ftpl1n*fsim2+ftpl2n] = getTanimotoAve166(ffpt1, ffpt2);
        
        fsim1[ftpl2n*fsim2+ftpl1n] = fsim1[ftpl1n*fsim2+ftpl2n];
       }
       
       ftpl2n++;
      }
    
    ftpl1n++;
   }
 
 int * fclu1 = new int [fsim2];
 
 int fclu2 = cluster_avelink( fsim1, fclu1, fsim2, fthr1, "max" );
 
 _ligclust_tot_mac = fclu2;
 
 for ( int fclu3 = 0; fclu3 < fclu2; fclu3++ )
 {
  int fclu4 = 0;
  double fbst1 = -1.0;
  std::bitset<MAXMAC> fbst2;
  std::string fbst3;
  
  ftpl1n = 0;
  
  for ( fipt1 = _complexes.begin(); fipt1 != _complexes.end(); fipt1++ )
   for ( int ilt1 = 0; ilt1 < (*fipt1)->getLigandsTotal(); ilt1++ )
    if ( (*fipt1)->getPocketNumber(ilt1) == _pocket_number )
    {
     if ( fclu1[ftpl1n] == fclu3 )
     {
      (*fipt1)->setPocketClusterNumberMACCS( ilt1, fclu3 + 1 );
      
      fclu4++;
      
      double ave1 = 0.0;
      double ave2 = 0.0;
      
      std::bitset<MAXMAC> ffpt1;
      
      (*fipt1)->getLigandFingerprintMACCS( ilt1, ffpt1 );
      
      int ftpl2n = 0;
      
      for ( fipt2 = _complexes.begin(); fipt2 != _complexes.end(); fipt2++ )
       for ( int ilt2 = 0; ilt2 < (*fipt2)->getLigandsTotal(); ilt2++ )
        if ( (*fipt2)->getPocketNumber(ilt2) == _pocket_number )
        {
         if ( fclu1[ftpl2n] == fclu3 )
         {
          if ( ftpl1n == ftpl2n )
           ave1 += 1.0;
          else
           ave1 += fsim1[ftpl1n*fsim2+ftpl2n];
          
          ave2 += 1.0;
         }
         
         ave1 -= 1.0;
         ave2 -= 1.0;
         
         if ( ave2 > 0.1 )
          ave1 /= ave2;
         else
          ave1 = 1.0;
         
         if ( ave1 > fbst1 )
         {
          fbst1 = ave1;
          fbst2 = ffpt1;
          fbst3 = (*fipt1)->getLigandID(ilt1);
         }
         
         ftpl2n++;
        }
     }
     
     ftpl1n++;
    }
  
  ligand_fpt_mac fclu5 = { fclu3 + 1, ( (double) fclu4 / (double) _complexes.size() ), fbst3, fbst2 };
  
  _ligand_fpts_mac.push_back( fclu5 );
 }
 
 delete [] fsim1;
 delete [] fclu1;
 
 for ( fipt1 = _complexes.begin(); fipt1 != _complexes.end(); fipt1++ )
  for ( int ilt1 = 0; ilt1 < (*fipt1)->getLigandsTotal(); ilt1++ )
   if ( (*fipt1)->getPocketNumber(ilt1) == _pocket_number )
   {
    std::bitset<MAXMAC> ffpt1;
    
    (*fipt1)->getLigandFingerprintMACCS( ilt1, ffpt1 );
    
    for ( int fipt3 = 0; fipt3 < MAXMAC; fipt3++ )
     _profile_fpts_mac[fipt3] += ( (double) ffpt1[MAXMAC-1-fipt3] / (double) _ligand_tot );
   }
}


// ==================================================================================   setCenter

void Pocket::setCenter( double tres1, double tdis1 )
{
 double tcen[3];
 
 for ( int itr1 = 0; itr1 < 3; itr1++ )
  tcen[itr1] = 0.0;
 
 double itr2 = 0.0;
 
 map<int,binding_residue>::iterator itr3;
 
 for ( itr3 = _binding_res.begin() ; itr3 != _binding_res.end(); itr3++ )
  if ( ((*itr3).second).residue_score >= tres1 )
  {
   itr2 += ((*itr3).second).residue_score;
   
   for ( int itr1 = 0; itr1 < 3; itr1++ )
    tcen[itr1] += ((*itr3).second).residue_score * ((*itr3).second).residue_coords[itr1];
  }
 
 for ( int itr1 = 0; itr1 < 3; itr1++ )
  tcen[itr1] /= itr2;
 
 if ( sqrt( pow( tcen[0] - _pocket_center[0], 2 ) + pow( tcen[1] - _pocket_center[1], 2 ) + pow( tcen[2] - _pocket_center[2], 2 ) ) > tdis1 )
 {
  for ( int itr1 = 0; itr1 < 3; itr1++ )
   _pocket_center[itr1] = tcen[itr1];
  
  for ( itr3 = _binding_res.begin() ; itr3 != _binding_res.end(); itr3++ )
   ((*itr3).second).residue_distance = sqrt( pow( ((*itr3).second).residue_coords[0] - _pocket_center[0], 2 ) + pow( ((*itr3).second).residue_coords[1] - _pocket_center[1], 2 ) + pow( ((*itr3).second).residue_coords[2] - _pocket_center[2], 2 ) );
 }
}


// ==================================================================================   dumpProteinsAlignments

void Pocket::dumpProteinsAlignments( std::string c1_name, map<string,bool> &chk2, Target * target2 )
{
 list< Template * >::iterator cipt1;
 
 for ( cipt1 = _complexes.begin(); cipt1 != _complexes.end(); cipt1++ )
  for ( int ilt1 = 0; ilt1 < (*cipt1)->getLigandsTotal(); ilt1++ )
   if ( (*cipt1)->getPocketNumber(ilt1) == _pocket_number )
    if ( chk2.count( (*cipt1)->getProteinID() ) < 1 )
    {
     chk2[(*cipt1)->getProteinID()] = true;
     
     (*cipt1)->dumpProtein( c1_name, true );
     
     double tarca1[MAXPRO][3];
     
     target2->getProteinCoordsCA(tarca1);
     
     (*cipt1)->dumpAlignment( c1_name, target2->getProteinResiduesTotal(), target2->getProteinSequence(), tarca1 );
    }
}


// ==================================================================================   dumpLigands

void Pocket::dumpLigands( std::string s1_name, map<string,bool> &chk3, int pktnum1 )
{
 list< Template * >::iterator fipt1;
 
 for ( int cipt1 = 0; cipt1 < _ligclust_tot_smi; cipt1++ )
  for ( fipt1 = _complexes.begin(); fipt1 != _complexes.end(); fipt1++ )
   for ( int ilt1 = 0; ilt1 < (*fipt1)->getLigandsTotal(); ilt1++ )
    if ( (*fipt1)->getPocketNumber(ilt1) == _pocket_number && (*fipt1)->getPocketClusterNumberSMILES(ilt1) == cipt1 + 1 )
     if ( chk3.count( (*fipt1)->getLigandID(ilt1) ) < 1 )
     {
      chk3[(*fipt1)->getLigandID(ilt1)] = true;
      
      (*fipt1)->dumpLigand( s1_name, ilt1, true, pktnum1 );
     }
}


// ==================================================================================   dumpPocket

void Pocket::dumpPocket( std::string p1_name, Target * target3, double tres1, int pktnum1 )
{
 ofstream outpkt1( (p1_name+".pockets.pdb").c_str(), ios_base::out|ios_base::app );
 ofstream outpkt2( (p1_name+".pockets.dat").c_str(), ios_base::out|ios_base::app );
 
 outpkt1 << "HETATM" << setw(5) << pktnum1
                     << "  EFS PKT M"
                     << setw(4) << _ligand_tot
                     << fixed << setw(12) << setprecision(3) << _pocket_center[0]
                     << fixed << setw(8) << setprecision(3) << _pocket_center[1]
                     << fixed << setw(8) << setprecision(3) << _pocket_center[2] << endl;
 
 outpkt2 << "POCKET" << setw(5) << pktnum1
                     << fixed << setw(5) << _complexes.size()
                     << fixed << setw(8) << setprecision(4) << _pocket_fraction
                     << fixed << setw(8) << setprecision(4) << _svm_confidence << endl;
 
 map<string,bool> chk2;
 
 list< Template * >::iterator itpl1;
 
 int itpl3 = 0;
 
 for ( itpl1 = _complexes.begin(); itpl1 != _complexes.end(); itpl1++ )
  if ( chk2.count( (*itpl1)->getProteinID() ) < 1 )
  {
   chk2[(*itpl1)->getProteinID()] = true;
   
   outpkt2 << "TEMPLATE" << setw(5) << ++itpl3 << setw(6) << (*itpl1)->getProteinID() 
                         << fixed << setw(5) << setprecision(0) << (*itpl1)->getProteinResiduesTotal()
                         << fixed << setw(5) << setprecision(0) << (*itpl1)->getProteinLengthTM()
                         << fixed << setw(8) << setprecision(3) << (*itpl1)->getProteinSeqID1()
                         << fixed << setw(8) << setprecision(3) << (*itpl1)->getProteinSeqID2()
                         << fixed << setw(8) << setprecision(3) << (*itpl1)->getProteinTMscore()
                         << fixed << setw(8) << setprecision(3) << (*itpl1)->getProteinRMSD()
                         << fixed << setw(8) << setprecision(3) << (*itpl1)->getProbPkt() << endl;
  }
 
 map<string,bool> chk3;
 
 for ( itpl1 = _complexes.begin(), itpl3 = 0; itpl1 != _complexes.end(); itpl1++, itpl3++ )
  for ( int itpl2 = 0; itpl2 < (*itpl1)->getLigandsTotal(); itpl2++ )
   if ( (*itpl1)->getPocketNumber(itpl2) == _pocket_number )
    if ( chk3.count( (*itpl1)->getLigandID(itpl2) ) < 1 )
    {
     chk3[(*itpl1)->getLigandID(itpl2)] = true;
     
     outpkt2 << "LIGAND" << setw(5) << itpl3 + 1 << setw(8) << (*itpl1)->getLigandID(itpl2) 
                       << fixed << setw(5) << setprecision(0) << (*itpl1)->getLigandAtomsTotal(itpl2)
                       << fixed << setw(10) << setprecision(3) << (*itpl1)->getLigandProp(itpl2,1)
                       << fixed << setw(8) << setprecision(3) << (*itpl1)->getLigandProp(itpl2,2)
                       << fixed << setw(9) << setprecision(2) << (*itpl1)->getLigandProp(itpl2,3)
                       << fixed << setw(8) << setprecision(2) << (*itpl1)->getLigandProp(itpl2,4)
                       << fixed << setw(4) << (int) (*itpl1)->getLigandProp(itpl2,5)
                       << fixed << setw(4) << (int) (*itpl1)->getLigandProp(itpl2,6)
                       << fixed << setw(8) << setprecision(3) << (*itpl1)->getProbLig() << endl;
    }
 
 outpkt2 << "PROP_MW   " << fixed << setw(10) << setprecision(3) << _ave_mw << fixed << setw(10) << setprecision(3) << _std_mw << endl
         << "PROP_LOGP " << fixed << setw(10) << setprecision(3) << _ave_logp << fixed << setw(10) << setprecision(3) << _std_logp << endl
         << "PROP_PSA  " << fixed << setw(10) << setprecision(3) << _ave_psa << fixed << setw(10) << setprecision(3) << _std_psa << endl
         << "PROP_MR   " << fixed << setw(10) << setprecision(3) << _ave_mr << fixed << setw(10) << setprecision(3) << _std_mr << endl
         << "PROP_HBD  " << fixed << setw(10) << setprecision(3) << _ave_hbd << fixed << setw(10) << setprecision(3) << _std_hbd << endl
         << "PROP_HBA  " << fixed << setw(10) << setprecision(3) << _ave_hba << fixed << setw(10) << setprecision(3) << _std_hba << endl;
 
 outpkt2 << "CENTER" << fixed << setw(9) << setprecision(3) << _pocket_center[0]
                     << fixed << setw(9) << setprecision(3) << _pocket_center[1]
                     << fixed << setw(9) << setprecision(3) << _pocket_center[2]
                     << fixed << setw(9) << setprecision(3) << _pocket_center[3] << endl;
 
 list<ligand_fpt_smi>::iterator itpl4;
 
 outpkt2 << "CLUSTERS" << fixed << setw(5) << _ligclust_tot_smi << fixed << setw(5) << _ligclust_tot_mac << endl;
 
 for ( itpl4 = _ligand_fpts_smi.begin(); itpl4 != _ligand_fpts_smi.end(); itpl4++ )
  outpkt2 << "SMILES" << setw(5) << (*itpl4).ftp_number 
                      << setw(8) << (*itpl4).fpt_pdbid
                      << fixed << setw(8) << setprecision(5) << (*itpl4).fpt_fraction
                      << setw(1025) << ((*itpl4).fpt_fingerprint).to_string<char,char_traits<char>,allocator<char> >() << endl;
 
 list<ligand_fpt_mac>::iterator itpl7;
 
 for ( itpl7 = _ligand_fpts_mac.begin(); itpl7 != _ligand_fpts_mac.end(); itpl7++ )
  outpkt2 << "MACCS " << setw(5) << (*itpl7).ftp_number 
                      << setw(8) << (*itpl7).fpt_pdbid
                      << fixed << setw(8) << setprecision(5) << (*itpl7).fpt_fraction
                      << setw(169) << ((*itpl7).fpt_fingerprint).to_string<char,char_traits<char>,allocator<char> >() << endl;
 
 outpkt2 << "PROFSMI" << flush;
 
 for ( int itpl6 = 0; itpl6 < MAXSMI; itpl6++ )
  outpkt2 << setw(5) << setprecision(2) << _profile_fpts_smi[itpl6] << flush;
 
 outpkt2 << endl;
 
 outpkt2 << "PROFMAC" << flush;
 
 for ( int itpl8 = 0; itpl8 < MAXMAC; itpl8++ )
  outpkt2 << setw(5) << setprecision(2) << _profile_fpts_mac[itpl8] << flush;
 
 outpkt2 << endl;
 
 map<int,binding_residue>::iterator itpl5;
 
 for ( itpl5 = _binding_res.begin() ; itpl5 != _binding_res.end(); itpl5++ )
 {
  string bind1 = " ";
  
  if ( ((*itpl5).second).residue_score >= tres1 )
   bind1 = "*";
  
  outpkt2 << "RESIDUE" << setw(6) << target3->getProteinNumbering(((*itpl5).second).residue_number)
                       << setw(2) << ((*itpl5).second).residue_code
                       << setw(2) << bind1
                       << setw(5) << ((*itpl5).second).template_incontact
                       << fixed << setw(9) << setprecision(3) << ((*itpl5).second).residue_distance
                       << fixed << setw(8) << setprecision(5) << ((*itpl5).second).template_fraction_bkg
                       << fixed << setw(8) << setprecision(5) << ((*itpl5).second).residue_score << flush;
  
  outpkt2 << endl;
 }
 
 outpkt2 << "PLBINDEX" << fixed << setw(6) << setprecision(3) << _plb << endl;
 
 outpkt2 << "AUXLIG" << fixed << setw(8) << setprecision(4) << _svm_cmps << endl;
 
 outpkt2 << "TER" << endl;
 
 outpkt1.close();
 outpkt2.close();
}


// ==================================================================================   setPocketFraction

void Pocket::setPocketFraction( double p1_fra )
{
 _pocket_fraction = p1_fra;
}


// ==================================================================================   calculateCmpsScores

void Pocket::calculateCmpsScores( Cmps * compounds1, ModelSVM * model1 )
{
 for ( int sc1 = 0; sc1 < compounds1->getCmpsTotal(); sc1++ )
 {
  double t_sco1[MAXSV4];
  
  for ( int sc2 = 0; sc2 < MAXSV4; sc2++ )
   t_sco1[sc2] = 0.0;
  
  bitset<MAXSMI> fpt1;
  
  compounds1->getSMILES(sc1, fpt1);
  
  list<ligand_fpt_smi>::iterator sc3;
  
  for ( sc3 = _ligand_fpts_smi.begin(); sc3 != _ligand_fpts_smi.end(); sc3++ )
  {
   t_sco1[0] += (*sc3).fpt_fraction * getTanimoto1024(fpt1, (*sc3).fpt_fingerprint);
   t_sco1[1] += (*sc3).fpt_fraction * getTanimotoAve1024(fpt1, (*sc3).fpt_fingerprint);
  }
  
  t_sco1[2] = getTanimotoCnt1024( fpt1, _profile_fpts_smi );
  
  bitset<MAXMAC> fpt2;
  
  compounds1->getMACCS(sc1, fpt2);
  
  list<ligand_fpt_mac>::iterator sc4;
  
  for ( sc4 = _ligand_fpts_mac.begin(); sc4 != _ligand_fpts_mac.end(); sc4++ )
  {
   t_sco1[3] += (*sc4).fpt_fraction * getTanimoto166(fpt2, (*sc4).fpt_fingerprint);
   t_sco1[4] += (*sc4).fpt_fraction * getTanimotoAve166(fpt2, (*sc4).fpt_fingerprint);
  }
  
  t_sco1[5] = getTanimotoCnt166( fpt2, _profile_fpts_mac );
  
  double t_mw = _std_mw;
  double t_logp = _std_logp;
  double t_psa = _std_psa;
  double t_hbd = _std_hbd;
  double t_hba = _std_hba;
  
  if ( t_mw < 0.1 ) t_mw = 0.1;
  if ( t_logp < 0.1 ) t_logp = 0.1;
  if ( t_psa < 0.1 ) t_psa = 0.1;
  if ( t_hbd < 0.1 ) t_hbd = 0.1;
  if ( t_hba < 0.1 ) t_hba = 0.1;
  
  t_sco1[6] = 0.5 * pow( ( _ave_mw - compounds1->getMW(sc1) ) / t_mw, 2.0) - log( 1.0 / ( t_mw * sqrt( 2.0 * PI ) ) );
  t_sco1[7] = 0.5 * pow( ( _ave_logp - compounds1->getLOGP(sc1) ) / t_logp, 2.0) - log( 1.0 / ( t_logp * sqrt( 2.0 * PI ) ) );
  t_sco1[8] = 0.5 * pow( ( _ave_psa - compounds1->getPSA(sc1) ) / t_psa, 2.0) - log( 1.0 / ( t_psa * sqrt( 2.0 * PI ) ) );
  t_sco1[9] = 0.5 * pow( ( _ave_hbd - (double) compounds1->getHBD(sc1) ) / t_hbd, 2.0) - log( 1.0 / ( t_hbd * sqrt( 2.0 * PI ) ) );
  t_sco1[10] = 0.5 * pow( ( _ave_hba - (double) compounds1->getHBA(sc1) ) / t_hba, 2.0) - log( 1.0 / ( t_hba * sqrt( 2.0 * PI ) ) );
  
  for ( int sc2 = 0; sc2 < MAXSV4; sc2++ )
   _score_cmps[sc2] += t_sco1[sc2];
 }
 
 for ( int sc2 = 0; sc2 < MAXSV4; sc2++ )
  _score_cmps[sc2] /= (double) compounds1->getCmpsTotal();
 
 _svm_cmps = model1->SVMpredict( 4, _score_cmps );
}


// ==================================================================================   calculateConfidence

double Pocket::calculateConfidence( bool opt1, ModelSVM * model1 )
{
 double a_tms1 = 0.0;
 double a_tms2 = 0.0;
 double a_pkt1 = 0.0;
 
 map<string,bool> chk2;
 
 list< Template * >::iterator itpl1;
 
 for ( itpl1 = _complexes.begin(); itpl1 != _complexes.end(); itpl1++ )
  if ( chk2.count( (*itpl1)->getProteinID() ) < 1 )
  {
   chk2[(*itpl1)->getProteinID()] = true;
   
   a_tms1 += (*itpl1)->getProteinTMscore();
   a_pkt1 += (*itpl1)->getProbPkt();
   
   a_tms2 += 1.0;
  }
 
 if ( a_tms2 > 0.0 )
 {
  a_tms1 /= a_tms2;
  a_pkt1 /= a_tms2;
 }
 
 if ( opt1 )
 {
  double t_sco1[MAXSV6];
  
  t_sco1[0] = _pocket_fraction;
  t_sco1[1] = log(_complexes.size());
  t_sco1[2] = a_tms1;
  t_sco1[3] = _ave_bres;
  t_sco1[4] = log(_binres_tot);
  t_sco1[5] = _plb;
  t_sco1[6] = a_pkt1;
  t_sco1[7] = _svm_cmps;
  
  _svm_confidence = model1->SVMpredict( 6, t_sco1 );
 }
 else
 {
  double t_sco1[MAXSV5];
  
  t_sco1[0] = _pocket_fraction;
  t_sco1[1] = log(_complexes.size());
  t_sco1[2] = a_tms1;
  t_sco1[3] = _ave_bres;
  t_sco1[4] = log(_binres_tot);
  t_sco1[5] = _plb;
  t_sco1[6] = a_pkt1;
  
  _svm_confidence = model1->SVMpredict( 5, t_sco1 );
 }
 
 return _svm_confidence;
}


// ==================================================================================   getConfidence

double Pocket::getConfidence( void )
{
 return _svm_confidence;
}
