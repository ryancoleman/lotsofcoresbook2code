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


#include "template.h"

using namespace std;

Template::Template( int aa, int ar, int al, double p1, double p2 )
{
 _na      = aa;
 _nr      = ar;
 _nl      = al;
 
 _prob_pkt = p1;
 _prob_lig = p2;
 
 _score_nw  = 0.0;
 _score_tm  = 0.0;
 _score_rms = 0.0;
 _score_id  = 0.0;
 
 _score_len_tm = 0;
 _score_len_nw = 0;
 
 for ( int i1 = 0; i1 < 3; i1++ )
 {
  _tm_t[i1] = 0.0;
  
  for ( int i2 = 0; i2 < 3; i2++ )
   _tm_u[i1][i2] = 0.0;
 }
 
 for ( int i1 = 0; i1 < MAXLIG; i1++ )
 {
  for ( int i2 = 0; i2 < 3; i2++ )
   _ligand_cen[i1][i2] = 0.0;
  
  _ligand_fpt_smi[i1].reset();
  _ligand_fpt_mac[i1].reset();
  
  _nla[i1] = 0;
  _nlb[i1] = 0;
  
  _ligand_mw[i1]   = 0.0;
  _ligand_logp[i1] = 0.0;
  _ligand_psa[i1]  = 0.0;
  _ligand_mr[i1]   = 0.0;
  _ligand_hbd[i1]  = 0;
  _ligand_hba[i1]  = 0;
 }
}

Template::Template( void )
{
 _na      = 0;
 _nr      = 0;
 _nl      = 0;
 
 _prob_pkt = 0.0;
 _prob_lig = 0.0;
 
 _score_nw  = 0.0;
 _score_tm  = 0.0;
 _score_rms = 0.0;
 _score_id  = 0.0;
 
 _score_len_tm = 0;
 _score_len_nw = 0;
 
 for ( int i1 = 0; i1 < 3; i1++ )
 {
  _tm_t[i1] = 0.0;
  
  for ( int i2 = 0; i2 < 3; i2++ )
   _tm_u[i1][i2] = 0.0;
 }
 
 for ( int i1 = 0; i1 < MAXLIG; i1++ )
 {
  for ( int i2 = 0; i2 < 3; i2++ )
   _ligand_cen[i1][i2] = 0.0;
  
  _ligand_fpt_smi[i1].reset();
  _ligand_fpt_mac[i1].reset();
  
  _nla[i1] = 0;
  _nlb[i1] = 0;
  
  _ligand_mw[i1]   = 0.0;
  _ligand_logp[i1] = 0.0;
  _ligand_psa[i1]  = 0.0;
  _ligand_mr[i1]   = 0.0;
  _ligand_hbd[i1]  = 0;
  _ligand_hba[i1]  = 0;
 }
}

Template::~Template() {}


// ==================================================================================   getLigandsTotal

int Template::getLigandsTotal( void )
{
 return _nl;
}


// ==================================================================================   getProteinResiduesTotal

int Template::getProteinResiduesTotal( void )
{
 return _nr;
}


// ==================================================================================   getProteinSequence

std::string Template::getProteinSequence( void )
{
 return _protein_seq1;
}


// ==================================================================================   getProteinAtomsTotal

int Template::getProteinAtomsTotal( void )
{
 return _na;
}


// ==================================================================================   getProteinID

string Template::getProteinID( void )
{
 return _protein_id;
}


// ==================================================================================   getLigandID

string Template::getLigandID( int lid1 )
{
 return _ligand_id[lid1];
}


// ==================================================================================   setProteinID

void Template::setProteinID( string pid1 )
{
 _protein_id = pid1;
}


// ==================================================================================   setLigandID

void Template::setLigandID( int lid1, string lid2 )
{
 _ligand_id[lid1] = lid2;
}


// ==================================================================================   getProteinCoordsCA

int Template::getProteinCoordsCA( double tabCA1[][3] )
{
 vector<CoordsProtein>::iterator ipca1;
 
 int ipca2 = 0;
 
 for ( ipca1 = _protein_xyz.begin(); ipca1 < _protein_xyz.end(); ipca1++ )
  if ( (*ipca1).getAtomName() == " CA " )
  {
   for ( int ipca3 = 0; ipca3 < 3; ipca3++ )
    tabCA1[ipca2][ipca3] = (*ipca1).getCoords(ipca3+1);
   
   ipca2++;
  }
 
 return ipca2;
}


// ==================================================================================   getProteinCoords1D

void Template::getProteinCoords1D( double tabCA1[] )
{
 vector<CoordsProtein>::iterator ipca1;
 
 int ipca2 = 0;
 
 for ( ipca1 = _protein_xyz.begin(); ipca1 < _protein_xyz.end(); ipca1++ )
  if ( (*ipca1).getAtomName() == " CA " )
   for ( int ipca3 = 0; ipca3 < 3; ipca3++ )
    tabCA1[ipca2++] = (*ipca1).getCoords(ipca3+1);
}


// ==================================================================================   setPocketNumber

void Template::setPocketNumber( int lnum1, int pnum1 )
{
 _pocket_num[lnum1] = pnum1;
}


// ==================================================================================   getPocketNumber

int Template::getPocketNumber( int lnum1 )
{
 return _pocket_num[lnum1];
}


// ==================================================================================   loadTemplate

bool Template::loadTemplate( std::string tpl_name, list<string> &template_list )
{
 string line1;
 
 igzstream tpl_file( tpl_name.c_str() );
 
 if ( ! tpl_file.good() )
 {
  std::cerr << "ERROR: Opening file `" << tpl_name << "' failed\n";
  return EXIT_FAILURE;
 }
 
 std::vector<string> tmp_sdf_t;
 
 while ( getline(tpl_file,line1) )
  if ( line1.length() > 3 )
  {
   if ( line1.compare(0, 3, "TAR") == 0 )
    _protein_id = line1.substr(4,5);
   
   else if ( line1.compare(0, 3, "PRT") == 0 )
   {
    string residue1 = line1.substr(21,3);
    string atom1    = line1.substr(16,4);
    
    int residue2 = atoi(line1.substr(26,4).c_str());
    int atom2    = atoi(line1.substr(10,5).c_str());
    
    double x1 = atof(line1.substr(34,8).c_str());
    double y1 = atof(line1.substr(42,8).c_str());
    double z1 = atof(line1.substr(50,8).c_str());
    
    _protein_xyz.push_back( CoordsProtein( atom2, residue2, x1, y1, z1, residue1, atom1 ) );
    
    _na++;
    
    if ( atom1 == " CA " )
    {
     _protein_seq1.append( three2oneS( residue1 ) );
     
     _protein_seq3[_nr] = residue2;
     
     _nr++;
    }
   }
   else if ( line1.compare(0, 3, "LGN") == 0 )
   {
    _ligand_sdf[_nl].push_back( line1.substr(4,line1.length()-4) );
    
    tmp_sdf_t.push_back(line1.substr(4,line1.length()-4));
    
    if ( line1.compare("LGN $$$$") == 0 )
    {
     _nla[_nl] = atoi(tmp_sdf_t[3].substr(0,3).c_str());
     _nlb[_nl] = atoi(tmp_sdf_t[3].substr(3,3).c_str());
     
     for ( int i1 = 0; i1 < 3; i1++ )
      _ligand_cen[_nl][i1] = 0.0;
     
     for ( int i1 = 4; i1 < _nla[_nl] + 4; i1++ )
     {
      double l1_x = atof(tmp_sdf_t[i1].substr(0,10).c_str());
      double l1_y = atof(tmp_sdf_t[i1].substr(10,10).c_str());
      double l1_z = atof(tmp_sdf_t[i1].substr(20,10).c_str());
      
      _ligand_cen[_nl][0] += l1_x;
      _ligand_cen[_nl][1] += l1_y;
      _ligand_cen[_nl][2] += l1_z;
      
      _ligand_xyz[_nl].push_back( CoordsLigand( i1-3, l1_x, l1_y, l1_z, tmp_sdf_t[i1].substr(30,3) ) );
     }
     
     for ( int i1 = 0; i1 < 3; i1++ )
      _ligand_cen[_nl][i1] /= _nla[_nl];
     
     for ( int i1 = 4 + _nla[_nl] + _nlb[_nl]; i1 < (int) tmp_sdf_t.size() - 1; i1++ )
     {
      if ( tmp_sdf_t[i1].find("FINGERPRINT") != string::npos )
      {
       string::iterator i2;
       int i3;
       
       for ( i2 = tmp_sdf_t[i1+1].begin(), i3 = 0; i2 < tmp_sdf_t[i1+1].end(); i2++, i3++ )
        switch(*i2)
        {
         case '0': _ligand_fpt_smi[_nl][MAXSMI-1-i3] = 0; break;
         case '1': _ligand_fpt_smi[_nl][MAXSMI-1-i3] = 1; break;
        }
      }
      
      else if ( tmp_sdf_t[i1].find("MACCS166") != string::npos )
      {
       string::iterator i2;
       int i3;
       
       for ( i2 = tmp_sdf_t[i1+1].begin(), i3 = 0; i2 < tmp_sdf_t[i1+1].end(); i2++, i3++ )
        switch(*i2)
        {
         case '0': _ligand_fpt_mac[_nl][MAXMAC-1-i3] = 0; break;
         case '1': _ligand_fpt_mac[_nl][MAXMAC-1-i3] = 1; break;
        }
      }
      
      else if ( tmp_sdf_t[i1].find("MOLID") != string::npos )
       _ligand_id[_nl] = tmp_sdf_t[i1+1];
      
      else if ( tmp_sdf_t[i1].find("OB_MW") != string::npos )
       _ligand_mw[_nl] = atof(tmp_sdf_t[i1+1].c_str());
      
      else if ( tmp_sdf_t[i1].find("OB_logP") != string::npos )
       _ligand_logp[_nl] = atof(tmp_sdf_t[i1+1].c_str());
      
      else if ( tmp_sdf_t[i1].find("OB_PSA") != string::npos )
       _ligand_psa[_nl] = atof(tmp_sdf_t[i1+1].c_str());
      
      else if ( tmp_sdf_t[i1].find("OB_MR") != string::npos )
       _ligand_mr[_nl] = atof(tmp_sdf_t[i1+1].c_str());
       
      else if ( tmp_sdf_t[i1].find("MCT_HBD") != string::npos )
       _ligand_hbd[_nl] = atoi(tmp_sdf_t[i1+1].c_str());
       
      else if ( tmp_sdf_t[i1].find("MCT_HBA") != string::npos )
       _ligand_hba[_nl] = atoi(tmp_sdf_t[i1+1].c_str());
     }
     
     tmp_sdf_t.clear();
     
     _nl++;
    }
   }
   else if ( line1.compare(0, 3, "ALN") == 0 )
   {
    istringstream iss1(line1);
    
    std::string iss2;
    
    int iss3 = 0;
    int iss4 = 0;
    
    std::string tmp_nam;
    
    tpl_alignment tmp_pair;
    
    while (iss1)
    {
     iss1 >> iss2;
     
     if ( iss3 == 1 )
     {
      list<string>::iterator it2;
      
      it2 = find(template_list.begin(), template_list.end(), iss2);
      
      if ( it2 != template_list.end() )
       iss4 = 1;
      
      tmp_nam = iss2;
     }
     else if ( iss3 > 1 && iss3 < 20 && iss4 > 0 )
      switch ( iss3 )
      {
       case 2  : tmp_pair.tpl_tms = atof(iss2.c_str()); break;
       case 3  : tmp_pair.tpl_rms = atof(iss2.c_str()); break;
       case 4  : tmp_pair.tpl_sid = atof(iss2.c_str()); break;
       case 5  : tmp_pair.tpl_len = atoi(iss2.c_str()); break;
       case 6  : tmp_pair.tpl_t[0] = atof(iss2.c_str()); break;
       case 7  : tmp_pair.tpl_u[0][0] = atof(iss2.c_str()); break;
       case 8  : tmp_pair.tpl_u[0][1] = atof(iss2.c_str()); break;
       case 9  : tmp_pair.tpl_u[0][2] = atof(iss2.c_str()); break;
       case 10 : tmp_pair.tpl_t[1] = atof(iss2.c_str()); break;
       case 11 : tmp_pair.tpl_u[1][0] = atof(iss2.c_str()); break;
       case 12 : tmp_pair.tpl_u[1][1] = atof(iss2.c_str()); break;
       case 13 : tmp_pair.tpl_u[1][2] = atof(iss2.c_str()); break;
       case 14 : tmp_pair.tpl_t[2] = atof(iss2.c_str()); break;
       case 15 : tmp_pair.tpl_u[2][0] = atof(iss2.c_str()); break;
       case 16 : tmp_pair.tpl_u[2][1] = atof(iss2.c_str()); break;
       case 17 : tmp_pair.tpl_u[2][2] = atof(iss2.c_str()); break;
       case 18 : tmp_pair.tpl_aln1 = iss2; break;
       case 19 : tmp_pair.tpl_aln1 = iss2; break;
      }
     
     iss3++;
    }
    
    if ( iss4 > 0 )
     _aligment_pair.insert( pair<string,tpl_alignment>(tmp_nam, tmp_pair) );
   }
  }
 
 tpl_file.close();
 
 strcpy(_protein_seq2, _protein_seq1.c_str());
 
 return EXIT_SUCCESS;
}


// ==================================================================================   alignNW

double Template::alignNW( std::string target_seq1 )
{
 double target_sco1, target_sco2;
 
 int target_len1 = _protein_seq1.size();
 int target_len2 = target_seq1.size();
 
 char target_seq2[MAXPRO];
 
 strcpy(target_seq2, target_seq1.c_str());
 
 int target_ali1[MAXPRO];
 
 nwalign_( &target_sco1, &target_sco2, &target_len1, &target_len2, &_protein_seq2, &target_seq2, &target_ali1 );
 
 for ( int ia1 = 0; ia1 < target_len2; ia1++ )
  _protein_align1[ia1] = target_ali1[ia1];
 
 _score_nw     = target_sco1;
 _score_len_nw = target_len2;
 
 return target_sco1;
}


// ==================================================================================   purgeAlignments

void Template::purgeAlignments( list<string> &template_list )
{
 map<string,tpl_alignment>::iterator it1 = _aligment_pair.begin();
 
 list<string>::iterator it2;
 
 while ( it1 != _aligment_pair.end() )
 {
  it2 = find (template_list.begin(), template_list.end(), (*it1).first);
  
  if ( it2 != template_list.end() )
   ++it1;
  else
   _aligment_pair.erase(it1++);
 }
}


// ==================================================================================   calculateContacts

void Template::calculateContacts( void )
{
 for ( int it1 = 0 ; it1 < _nl; it1++ )
 {
  map<int,double> bres1;
  map<int,char> bres2;
  multimap<int,string> bres3;
  
  vector<CoordsLigand>::iterator it2;
  
  for ( it2 = _ligand_xyz[it1].begin(); it2 < _ligand_xyz[it1].end(); it2++ )
  {
   vector<CoordsProtein>::iterator it3;
   
   for ( it3 = _protein_xyz.begin(); it3 < _protein_xyz.end(); it3++ )
   {
    double rr1 = sqrt( pow((*it2).getCoords(1)-(*it3).getCoords(1),2.0) + 
                       pow((*it2).getCoords(2)-(*it3).getCoords(2),2.0) + 
                       pow((*it2).getCoords(3)-(*it3).getCoords(3),2.0) );
    
    if ( rr1 <= 4.5 )
    {
     lig_contact tmp_cnt;
     
     tmp_cnt.residue_number   = (*it3).getResidueNumber();
     tmp_cnt.residue_name     = three2oneC( (*it3).getResidueName() );
     tmp_cnt.residue_atom     = (*it3).getAtomName();
     
     tmp_cnt.atom_number      = (*it2).getAtomNumber();
     tmp_cnt.atom_name        = (*it2).getAtomName();
     tmp_cnt.contact_distance = rr1;
     
     _contacts_ligand[it1].push_back( tmp_cnt );
     
     map<int,double>::iterator it4;
     
     it4 = bres1.find( (*it3).getResidueNumber() );
     
     if ( it4 == bres1.end() )
     {
      bres1.insert( pair<int,double>((*it3).getResidueNumber(), rr1) );
      
      bres2.insert( pair<int,char>((*it3).getResidueNumber(), three2oneC( (*it3).getResidueName() )) );
     }
     else
      if ( rr1 < (*it4).second )
       (*it4).second = rr1;
     
     bres3.insert( pair<int,string>((*it3).getResidueNumber(), (*it3).getAtomName()) );
    }
   }
  }
  
  map<int,double>::iterator it5;
  
  for ( it5 = bres1.begin(); it5 != bres1.end(); it5++ )
  {
   map<int,char>::iterator it6;
   
   it6 = bres2.find( (*it5).first );
   
   lig_binding tmp_res;
   
   tmp_res.residue_number = (*it5).first - 1;
   tmp_res.residue_name = (*it6).second;
   tmp_res.residue_sidechain = false;
   tmp_res.contact_distance = (*it5).second;
   
   multimap<int,string>::iterator it7;
   pair<multimap<int,string>::iterator,multimap<int,string>::iterator> it8;
   
   it8 = bres3.equal_range( (*it5).first );
   
   for ( it7 = it8.first; it7 != it8.second; it7++ )
    if ( (*it7).second.compare(" N  ") != 0 && (*it7).second.compare(" CA ") != 0 && (*it7).second.compare(" C  ") != 0 && (*it7).second.compare(" O  ") != 0 )
     tmp_res.residue_sidechain = true;
   
   _bindings_ligand[it1].push_back( tmp_res );
  }
 }
}


// ==================================================================================   getLigandCenter

void Template::getLigandCenter( int ipc1, double tab1[], bool rot1 )
{
 if ( rot1 )
  for ( int ipc2 = 0; ipc2 < 3; ipc2++ )
   tab1[ipc2] = _tm_t[ipc2] + _tm_u[ipc2][0] * _ligand_cen[ipc1][0] + _tm_u[ipc2][1] * _ligand_cen[ipc1][1] + _tm_u[ipc2][2] * _ligand_cen[ipc1][2];
 
 else
  for ( int ipc2 = 0; ipc2 < 3; ipc2++ )
   tab1[ipc2] = _ligand_cen[ipc1][ipc2];
}


// ==================================================================================   getMatrix

bool Template::getMatrix( std::string pid1, double tt1[], double uu1[][3] )
{
 map<string,tpl_alignment>::iterator it1;
 
 it1 = _aligment_pair.find( pid1 );
 
 if ( it1 == _aligment_pair.end() )
  return EXIT_FAILURE;
 
 for ( int i1 = 0; i1 < 3; i1++ )
 {
  tt1[i1] = (*it1).second.tpl_t[i1];
  
  for ( int i2 = 0; i2 < 3; i2++ )
   uu1[i1][i2] = (*it1).second.tpl_u[i1][i2];
 }
 
 return EXIT_SUCCESS;
}


// ==================================================================================   getNWalignment

void Template::getNWalignment( int tab1[], int tab2 )
{
 for ( int ial1 = 0; ial1 < tab2; ial1++ )
  tab1[ial1] = _protein_align1[ial1];
}


// ==================================================================================   getTMalignment

void Template::getTMalignment( int tab1[], int tab2 )
{
 for ( int ial1 = 0; ial1 < tab2; ial1++ )
  tab1[ial1] = _protein_align2[ial1];
}


// ==================================================================================   setTMalignment

void Template::setTMalignment( int tab1[], int tab2 )
{
 for ( int ial1 = 0; ial1 < tab2; ial1++ )
  _protein_align2[ial1] = tab1[ial1];
}


// ==================================================================================   setMatrix

void Template::setMatrix( double tt1[], double uu1[][3] )
{
 for ( int i1 = 0; i1 < 3; i1++ )
 {
  _tm_t[i1] = tt1[i1];
  
  for ( int i2 = 0; i2 < 3; i2++ )
   _tm_u[i1][i2] = uu1[i1][i2];
 }
}


// ==================================================================================   getBindingResidues

void Template::getBindingResidues( int lig1, list<lig_binding> &res1 )
{
 list<lig_binding>::iterator it1;
 
 for ( it1 = _bindings_ligand[lig1].begin(); it1 != _bindings_ligand[lig1].end(); it1++ )
  res1.push_back( *it1 );
}


// ==================================================================================   getProbPkt

double Template::getProbPkt( void )
{
 return _prob_pkt;
}


// ==================================================================================   getProbLig

double Template::getProbLig( void )
{
 return _prob_lig;
}


// ==================================================================================   getLigandProp

double Template::getLigandProp( int al, int ap )
{
 switch (ap)
 {
  case  1  :  return _ligand_mw[al];
  case  2  :  return _ligand_logp[al];
  case  3  :  return _ligand_psa[al];
  case  4  :  return _ligand_mr[al];
  case  5  :  return (double) _ligand_hbd[al];
  case  6  :  return (double) _ligand_hba[al];
  
  default  :  return 0;
 }
}


// ==================================================================================   getLigandFingerprintSMILES

void Template::getLigandFingerprintSMILES ( int tfpt1, bitset<MAXSMI> &tfpt2 )
{
 tfpt2 = _ligand_fpt_smi[tfpt1];
}


// ==================================================================================   getLigandFingerprintMACCS

void Template::getLigandFingerprintMACCS ( int tfpt1, bitset<MAXMAC> &tfpt2 )
{
 tfpt2 = _ligand_fpt_mac[tfpt1];
}


// ==================================================================================   getPocketClusterNumberSMILES

int Template::getPocketClusterNumberSMILES( int pclu1 )
{
 return _pocket_clu_smi[pclu1];
}


// ==================================================================================   setPocketClusterNumberSMILES

void Template::setPocketClusterNumberSMILES( int pclu1, int pclu2 )
{
 _pocket_clu_smi[pclu1] = pclu2;
}


// ==================================================================================   getPocketClusterNumberMACCS

int Template::getPocketClusterNumberMACCS( int pclu1 )
{
 return _pocket_clu_mac[pclu1];
}


// ==================================================================================   setPocketClusterNumberMACCS

void Template::setPocketClusterNumberMACCS( int pclu1, int pclu2 )
{
 _pocket_clu_mac[pclu1] = pclu2;
}


// ==================================================================================   dumpProtein

void Template::dumpProtein( std::string p1_name, bool p1_align )
{
 vector<CoordsProtein>::iterator p1_i;
 
 ofstream outprot( (p1_name+".templates.pdb").c_str(), ios_base::out|ios_base::app );
 
 outprot << "REMARK   PDB-ID:  " << setw(8) << _protein_id << endl
         << "REMARK   SEQID1:  " << fixed << setw(8) << setprecision(3) << _score_nw << endl
         << "REMARK   SEQID2:  " << fixed << setw(8) << setprecision(3) << _score_id << endl
         << "REMARK   TM-SCORE:" << fixed << setw(8) << setprecision(3) << _score_tm << endl
         << "REMARK   RMSD:    " << fixed << setw(8) << setprecision(3) << _score_rms << endl;
 
 
 
 for ( p1_i = _protein_xyz.begin(); p1_i < _protein_xyz.end(); p1_i++ )
 {
  double xyz1[3];
  
  if ( p1_align )
   for ( int ipc1 = 0; ipc1 < 3; ipc1++ )
    xyz1[ipc1] = _tm_t[ipc1] + _tm_u[ipc1][0] * (*p1_i).getCoords(1) + _tm_u[ipc1][1] * (*p1_i).getCoords(2) + _tm_u[ipc1][2] * (*p1_i).getCoords(3);
  
  else
   for ( int ipc1 = 0; ipc1 < 3; ipc1++ )
    xyz1[ipc1] = (*p1_i).getCoords(ipc1+1);
  
  outprot << "ATOM  " << setw(5) << (*p1_i).getAtomNumber()
                      << setw(5) << (*p1_i).getAtomName()
                      << setw(4) << (*p1_i).getResidueName()
                      << setw(6) << (*p1_i).getResidueNumber()
                      << fixed << setw(12) << setprecision(3) << xyz1[0]
                      << fixed << setw(8)  << setprecision(3) << xyz1[1]
                      << fixed << setw(8)  << setprecision(3) << xyz1[2] << endl;
 }
 
 outprot << "TER" << endl;
 
 outprot.close();
}


// ==================================================================================   dumpAlignment

void Template::dumpAlignment( std::string p1_name, int tarlen1, std::string tarseq1, double tarca2[][3] )
{
 list<int> frsq1;
 list<int> frsq2;
 
 list<int>::iterator isq1;
 list<int>::iterator isq2;
 
 std::string asq1 = "";
 std::string asq2 = "";
 std::string asq3 = "";
 
 for ( int iaf1 = 0; iaf1 < tarlen1; iaf1++ )
 {
  if ( iaf1 == 0 )
  {
   if ( _protein_align2[iaf1] > 0 )
   {
    for ( int iaf2 = 0; iaf2 < _protein_align2[iaf1]; iaf2++ )
    {
     frsq1.push_back(-1);
     frsq2.push_back(iaf2);
    }
   }
   else if ( _protein_align2[iaf1] < 0 )
   {
    int iaf3 = -1;
    
    for ( int iaf2 = 0; iaf2 < tarlen1; iaf2++ )
    {
     if ( _protein_align2[iaf2] > -1 )
     {
      iaf3 = _protein_align2[iaf2];
      
      break;
     }
    }
    
    for ( int iaf2 = 0; iaf2 < iaf3; iaf2++ )
    {
     frsq1.push_back(-1);
     frsq2.push_back(iaf2);
    }
   }
   
   frsq1.push_back(iaf1);
   frsq2.push_back(_protein_align2[iaf1]);
  }
  else if ( iaf1 == tarlen1 -1 )
  {
   frsq1.push_back(iaf1);
   frsq2.push_back(_protein_align2[iaf1]);
   
   if ( _protein_align2[iaf1] >= 0 )
   {
    if ( _protein_align2[iaf1] + 1 < _nr )
    {
     for ( int iaf2 = _protein_align2[iaf1] + 1; iaf2 < _nr; iaf2++ )
     {
      frsq1.push_back(-1);
      frsq2.push_back(iaf2);
     }
    }
   }
   else
   {
    int iaf3 = -1;
    
    for ( int iaf2 = iaf1; iaf2 >= 0 ; iaf2-- )
    {
     if ( _protein_align2[iaf2] > -1 )
     {
      iaf3 = _protein_align2[iaf2];
      
      break;
     }
    }
    
    for ( int iaf2 = 1; iaf2 < _nr - iaf3; iaf2++ )
    {
     frsq1.push_back(-1);
     frsq2.push_back(iaf3+iaf2);
    }
   }
  }
  else
  {
   frsq1.push_back(iaf1);
   frsq2.push_back(_protein_align2[iaf1]);
   
   if ( _protein_align2[iaf1+1] < 0 && _protein_align2[iaf1] > -1 )
   {
    int iaf3 = -1;
    
    for ( int iaf2 = iaf1 + 1; iaf2 < tarlen1; iaf2++ )
    {
     if ( _protein_align2[iaf2] > -1 )
     {
      iaf3 = _protein_align2[iaf2];
      
      break;
     }
    }
    
    for ( int iaf2 = 1; iaf2 < iaf3 - _protein_align2[iaf1]; iaf2++ )
    {
     frsq1.push_back(-1);
     frsq2.push_back(_protein_align2[iaf1] + iaf2);
    }
   }
   else if ( _protein_align2[iaf1] + 1 < _protein_align2[iaf1+1] && _protein_align2[iaf1] > -1 )
   {
    for ( int iaf2 = 1; iaf2 < ( _protein_align2[iaf1+1] - _protein_align2[iaf1] ); iaf2++ )
    {
     frsq1.push_back(-1);
     frsq2.push_back(_protein_align2[iaf1] + iaf2);
    }
   }
  }
 }
 
 double tplca2[MAXPRO][3];
 
 vector<CoordsProtein>::iterator itca1;
 
 int itca2 = 0;
 
 for ( itca1 = _protein_xyz.begin(); itca1 < _protein_xyz.end(); itca1++ )
  if ( (*itca1).getAtomName() == " CA " )
  {
   for ( int itca3 = 0; itca3 < 3; itca3++ )
    tplca2[itca2][itca3] = _tm_t[itca3] + _tm_u[itca3][0] * (*itca1).getCoords(1) + _tm_u[itca3][1] * (*itca1).getCoords(2) + _tm_u[itca3][2] * (*itca1).getCoords(3);
   
   itca2++;
  }
 
 for ( isq1 = frsq1.begin(), isq2 = frsq2.begin(); isq1 != frsq1.end(); isq1++, isq2++ )
 {
  if ( *isq1 < 0 )
   asq1 += "-";
  else
   asq1 += tarseq1.substr(*isq1,1).c_str();
  
  if ( *isq2 < 0 )
   asq2 += "-";
  else
   asq2 += _protein_seq1.substr(*isq2,1).c_str();
  
  if ( *isq1 > -1 && *isq2 > -1 )
  {
   if ( sqrt( pow(tarca2[*isq1][0]-tplca2[*isq2][0], 2) + pow(tarca2[*isq1][1]-tplca2[*isq2][1], 2) + pow(tarca2[*isq1][2]-tplca2[*isq2][2], 2)) < 5.0 )
    asq3 += ":";
   else
    asq3 += ".";
  }
  else
   asq3 += " ";
 }
 
 ofstream outali( (p1_name+".alignments.dat").c_str(), ios_base::out|ios_base::app );
 
 outali << ">"  << _protein_id << " " << _nr << " " << _score_len_tm << " " << fixed << setprecision(3) << _score_tm << " " << _score_rms << " " << _score_id << endl << asq1 << endl << asq3 << endl << asq2 << endl << "*" << endl;
 
 outali.close();
}


// ==================================================================================   dumpLigand

void Template::dumpLigand( std::string l1_name, int l1_num, bool l1_align, int l1_pocket )
{
 vector<CoordsLigand>::iterator il1;
 
 std::vector<double> l1_xyz;
 
 for ( il1 = _ligand_xyz[l1_num].begin(); il1 < _ligand_xyz[l1_num].end(); il1++ )
 {
  if ( l1_align )
   for ( int ipc1 = 0; ipc1 < 3; ipc1++ )
    l1_xyz.push_back(_tm_t[ipc1] + _tm_u[ipc1][0] * (*il1).getCoords(1) + _tm_u[ipc1][1] * (*il1).getCoords(2) + _tm_u[ipc1][2] * (*il1).getCoords(3));
   
  else
   for ( int ipc1 = 0; ipc1 < 3; ipc1++ )
    l1_xyz.push_back((*il1).getCoords(ipc1+1));
 }
 
 list<std::string>::iterator il2;
 
 std::string pat1;
 
 int il3;
 
 for ( il2 = _ligand_sdf[l1_num].begin(), il3 = 0; il2 != _ligand_sdf[l1_num].end(); il2++, il3++ )
 {
  if ( il3 > 3 && il3 < _nla[l1_num] + 4 )
  {
   stringstream il4;
   
   il4 << fixed << setw(10) << setprecision(4) << l1_xyz[(il3-4)*3+0]
       << fixed << setw(10) << setprecision(4) << l1_xyz[(il3-4)*3+1]
       << fixed << setw(10) << setprecision(4) << l1_xyz[(il3-4)*3+2];
   
   (*il2).replace( 0, 30, il4.str() );
  }
  
  else if ( (*il2).find("MOLID") != string::npos )
   pat1 = (*il2);
 }
 
 l1_xyz.clear();
 
 std::stringstream pkt_n1;
 std::stringstream pkt_s1;
 std::stringstream pkt_m1;
 
 pkt_n1 << l1_pocket;
 pkt_s1 << _pocket_clu_smi[l1_num];
 pkt_m1 << _pocket_clu_mac[l1_num];
 
 ofstream outlig( (l1_name+".ligands.sdf").c_str(), ios_base::out|ios_base::app );
 
 for ( il2 = _ligand_sdf[l1_num].begin(); il2 != _ligand_sdf[l1_num].end(); il2++ )
 {
  if ( (*il2).find("$$$$") != string::npos )
  {
   std::string pat2 = pat1;
   
   pat2.replace(pat2.find("MOLID"),5,"EFINDSITE_POCKET");
   
   std::string pat3 = pat1;
   
   pat3.replace(pat3.find("MOLID"),5,"FINDSITE_CLUSTER_SMILES");
   
   std::string pat4 = pat1;
   
   pat4.replace(pat4.find("MOLID"),5,"FINDSITE_CLUSTER_MACCS");
   
   outlig << pat2 << endl
          << pkt_n1.str() << endl << endl
          << pat3 << endl
          << pkt_s1.str() << endl << endl
          << pat4 << endl
          << pkt_m1.str() << endl << endl;
  }
  
  outlig << (*il2) << endl;
 }
 
 outlig.close();
}


// ==================================================================================   getProteinSeqID1

double Template::getProteinSeqID1( void )
{
 return _score_nw;
}


// ==================================================================================   getProteinSeqID2

double Template::getProteinSeqID2( void )
{
 return _score_id;
}


// ==================================================================================   getProteinTMscore

double Template::getProteinTMscore( void )
{
 return _score_tm;
}


// ==================================================================================   getProteinRMSD

double Template::getProteinRMSD( void )
{
 return _score_rms;
}


// ==================================================================================   getProteinLengthNW

int Template::getProteinLengthNW( void )
{
 return _score_len_nw;
}


// ==================================================================================   getProteinLengthTM

int Template::getProteinLengthTM( void )
{
 return _score_len_tm;
}


// ==================================================================================   setProteinSeqID1

void Template::setProteinSeqID1( double nw1 )
{
 _score_nw = nw1;
}


// ==================================================================================   setProteinSeqID2

void Template::setProteinSeqID2( double id1 )
{
 _score_id = id1;
}


// ==================================================================================   setProteinTMscore

void Template::setProteinTMscore( double tm1 )
{
 _score_tm = tm1;
}


// ==================================================================================   setProteinRMSD

void Template::setProteinRMSD( double rms1 )
{
 _score_rms = rms1;
}


// ==================================================================================   setProteinLengthNW

void Template::setProteinLengthNW( int lnw1 )
{
 _score_len_nw = lnw1;
}


// ==================================================================================   setProteinLengthTM

void Template::setProteinLengthTM( int ltm1 )
{
 _score_len_tm = ltm1;
}


// ==================================================================================   getLigandAtomsTotal

int Template::getLigandAtomsTotal( int l1_num )
{
 return _nla[l1_num];
}


// ==================================================================================   getLigandBondsTotal

int Template::getLigandBondsTotal( int l1_num )
{
 return _nlb[l1_num];
}
