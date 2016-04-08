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


#include "efindsite_scr.h"

using namespace std;

int main(int argc, char *argv[])
{
 time_t t_start, t_end, t_bench1, t_bench2;
 time(&t_start);
 
 std::string smethod = "sum";
 int pktnum = 1;
 double ctres1 = 0.571;
 double ctres2 = 0.592;
 
 cout << "------------------------------------------------------------" << endl
      << "                       efindsite_scr" << endl
      << "                        version 1.2" << endl
      << "               ligand-based virtual screening" << endl << endl
      << "       report bugs and issues to michal@brylinski.org" << endl
      << "------------------------------------------------------------" << endl << endl;
 
 if ( argc < 7 )
 {
  cout << " efindsite_scr -p <eFindSite pockets>" << endl
       << "               -s <screening library, gzipped>" << endl
       << "               -o <output filename>" << endl << endl
       << " additional options:" << endl << endl
       << "               -n <pocket number (default 1)>" << endl
       << "               -m <scoring function (default sum)>" << endl << endl
       << " available single scoring functions:" << endl << endl
       << "               tst - classical Tanimoto coeff for Daylight" << endl
       << "               tsa - average Tanimoto coeff for Daylight" << endl
       << "               tsc - continuous Tanimoto coeff for Daylight" << endl
       << "               tmt - classical Tanimoto coeff for MACCS" << endl
       << "               tma - average Tanimoto coeff for MACCS" << endl
       << "               tmc - continuous Tanimoto coeff for MACCS" << endl << endl
       << " available combined scoring functions:" << endl << endl
       << "               sum - data fusion using sum rule" << endl
       << "               max - data fusion using max rule" << endl
       << "               min - data fusion using min rule" << endl
       << "               svm - machine learning using SVM" << endl << endl;
  
  exit(EXIT_SUCCESS);
 }
 
 string pocket_name;
 bool pocket_opt = false;
 
 string library_name;
 bool library_opt = false;
 
 string output_name;
 bool output_opt = false;
 
 for ( int i = 0; i < argc; i++ )
 {
  if ( !strcmp(argv[i],"-p") && i < argc ) { pocket_name  = string(argv[i+1]); pocket_opt  = true; }
  if ( !strcmp(argv[i],"-s") && i < argc ) { library_name = string(argv[i+1]); library_opt = true; }
  if ( !strcmp(argv[i],"-o") && i < argc ) { output_name  = string(argv[i+1]); output_opt  = true; }
  if ( !strcmp(argv[i],"-n") && i < argc ) { pktnum       = atoi(argv[i+1]);                       }
  if ( !strcmp(argv[i],"-m") && i < argc ) { smethod      = string(argv[i+1]);                     }
 }
 
 char * path1;
 
 path1 = getenv("EF_MOD"); if ( path1==NULL ) { cout << "EF_MOD is not set" << endl; exit(EXIT_FAILURE); }
 
 string model_path;
 model_path = getenv("EF_MOD");
 
 ifstream f01( (model_path+"/screen01SVM.model").c_str() );
 ifstream f02( (model_path+"/screen10SVM.model").c_str() );
 ifstream f03( (model_path+"/screenSVM.scale").c_str() );
 ifstream f04( (model_path+"/scoringSVM.model").c_str() );
 ifstream f05( (model_path+"/scoringSVM.scale").c_str() );
 
 if ( !f01 || !f02 || !f03 || !f04 || !f05 )
 {
  cout << "Could not find SVM models in " << model_path << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( smethod != "tst" && smethod != "tsa" && smethod != "tsc" && 
      smethod != "tmt" && smethod != "tma" && smethod != "tmc" && 
      smethod != "sum" && smethod != "max" && smethod != "min" && 
      smethod != "svm" )
 {
  cout << "Invalid scoring method: " << smethod << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !pocket_opt )
 {
  cout << "Provide eFindSite pockets" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !library_opt )
 {
  cout << "Provide screening library" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !output_opt )
 {
  cout << "Provide output filename" << endl;
  exit(EXIT_FAILURE);
 }
 else
 {
  ofstream outrank( (output_name).c_str() );
  outrank.close();
 }
 
 ModelSVM * model_svm;
 
 model_svm = new ModelSVM( false, false, false, false, false, false, false, false, false );
 
 cout << "Loading SVM models " << flush;
 
 time(&t_bench1);
 
 model_svm->loadModel( 7, model_path+"/screen01SVM.model" ); cout << '.' << flush;
 model_svm->loadModel( 8, model_path+"/screen10SVM.model" ); cout << '.' << flush;
 
 model_svm->loadScale( 7, model_path+"/screenSVM.scale" ); cout << '.' << flush;
 model_svm->loadScale( 8, model_path+"/screenSVM.scale" ); cout << '.' << flush;
 
 if ( smethod == "svm" )
 {
  model_svm->loadModel( 9, model_path+"/scoringSVM.model" ); cout << '.' << flush;
  
  model_svm->loadScale( 9, model_path+"/scoringSVM.scale" ); cout << '.' << flush;
 }
 
 time(&t_bench2);
 
 cout << " done (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 string line1;
 
 bool pkt1 = false;
 
 double profile1[MAXSMI];
 double profile2[MAXMAC];
 
 list<cmp1> cmps1;
 list<cmp2> cmps2;
 
 double mprop[5][2];
 
 for ( int ii1 = 0; ii1 < 5; ii1++ )
  for ( int ii2 = 0; ii2 < 2; ii2++ )
   mprop[ii1][ii2] = 0.0;
 
 ifstream pocket_file( pocket_name.c_str() );
 
 if ( !pocket_file.is_open() )  { cout << "Cannot open " << pocket_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(pocket_file,line1))
 {
  if ( line1.size() > 31 )
   if ( line1.substr(0,6) == "POCKET" )
    if ( atoi(line1.substr(7,4).c_str()) == pktnum )
     pkt1 = true;
  
  if ( line1 == "TER" )
   pkt1 = false;
  
  if ( pkt1 )
  {
   if ( line1.size() > 1051 )
    if ( line1.substr(0,6) == "SMILES" )
    {
     double lig1      = atof(line1.substr(19,8).c_str());
     std::string lig2 = line1.substr(28,MAXSMI).c_str();
     
     bitset<MAXSMI> lig3;
     
     string::iterator pkt2;
     int pkt3;
     
     for ( pkt2 = lig2.begin(), pkt3 = 0; pkt2 < lig2.end(); pkt2++, pkt3++ )
      switch(*pkt2)
      {
       case '0': lig3[MAXSMI-1-pkt3] = 0; break;
       case '1': lig3[MAXSMI-1-pkt3] = 1; break;
      }
     
     cmp1 lig4 = { lig1, lig3 };
     
     cmps1.push_back( lig4 );
    }
   
   if ( line1.size() > 195 )
    if ( line1.substr(0,5) == "MACCS" )
    {
     double lig1      = atof(line1.substr(19,8).c_str());
     std::string lig2 = line1.substr(28,MAXMAC).c_str();
     
     bitset<MAXMAC> lig3;
     
     string::iterator pkt2;
     int pkt3;
     
     for ( pkt2 = lig2.begin(), pkt3 = 0; pkt2 < lig2.end(); pkt2++, pkt3++ )
      switch(*pkt2)
      {
       case '0': lig3[MAXMAC-1-pkt3] = 0; break;
       case '1': lig3[MAXMAC-1-pkt3] = 1; break;
      }
     
     cmp2 lig4 = { lig1, lig3 };
     
     cmps2.push_back( lig4 );
    }
   
   if ( line1.size() > 5126 )
    if ( line1.substr(0,7) == "PROFSMI" )
     for ( int prf1 = 0; prf1 < MAXSMI; prf1++ )
      profile1[prf1] = atof(line1.substr(7+prf1*5,5).c_str());
   
   if ( line1.size() > 846 )
    if ( line1.substr(0,7) == "PROFMAC" )
     for ( int prf1 = 0; prf1 < MAXMAC; prf1++ )
      profile2[prf1] = atof(line1.substr(7+prf1*5,5).c_str());
   
   if ( line1.size() > 29 )
   {
         if ( line1.substr(0,9) == "PROP_MW  " )
    {
     mprop[0][0] = atof(line1.substr(10,10).c_str());
     mprop[0][1] = atof(line1.substr(20,10).c_str());
    }
    else if ( line1.substr(0,9) == "PROP_LOGP" )
    {
     mprop[1][0] = atof(line1.substr(10,10).c_str());
     mprop[1][1] = atof(line1.substr(20,10).c_str());
    }
    else if ( line1.substr(0,9) == "PROP_PSA " )
    {
     mprop[2][0] = atof(line1.substr(10,10).c_str());
     mprop[2][1] = atof(line1.substr(20,10).c_str());
    }
    else if ( line1.substr(0,9) == "PROP_HBD " )
    {
     mprop[3][0] = atof(line1.substr(10,10).c_str());
     mprop[3][1] = atof(line1.substr(20,10).c_str());
    }
    else if ( line1.substr(0,9) == "PROP_HBA " )
    {
     mprop[4][0] = atof(line1.substr(10,10).c_str());
     mprop[4][1] = atof(line1.substr(20,10).c_str());
    }
   }
  }
 }
 
 pocket_file.close();
 
 string line2;
 
 igzstream library_file( library_name.c_str() );
 
 if ( ! library_file.good() )
 {
  std::cerr << "ERROR: Opening file `" << library_name << "' failed\n";
  return EXIT_FAILURE;
 }
 
 map< std::string , cmp3 > library;
 
 std::string lib2[8];
 
 int lib8 = 0;
 
 while ( getline(library_file,line2) )
  if ( lib8++ < MAXSCR )
  {
   int lib3 = 0;
   
   istringstream lib4(line2);
   
   while (lib4)
    lib4 >> lib2[lib3++];
   
   double tmw   = atof(lib2[1].c_str());
   double tlogp = atof(lib2[2].c_str());
   double tpsa  = atof(lib2[3].c_str());
   
   int thbd = atoi(lib2[4].c_str());
   int thba = atoi(lib2[5].c_str());
   
   bitset<MAXSMI> tfpt1;
   bitset<MAXMAC> tfpt2;
   
   string::iterator lib6;
   
   int lib7;
   
   for ( lib6 = lib2[6].begin(), lib7 = 0; lib6 < lib2[6].end(); lib6++, lib7++ )
    switch(*lib6)
    {
     case '0': tfpt1[MAXSMI-1-lib7] = 0; break;
     case '1': tfpt1[MAXSMI-1-lib7] = 1; break;
    }
   
   for ( lib6 = lib2[7].begin(), lib7 = 0; lib6 < lib2[7].end(); lib6++, lib7++ )
    switch(*lib6)
    {
     case '0': tfpt2[MAXMAC-1-lib7] = 0; break;
     case '1': tfpt2[MAXMAC-1-lib7] = 1; break;
    }
   
   cmp3 tcmp = { tmw, tlogp, tpsa, thbd, thba, tfpt1, tfpt2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   
   library.insert( pair< std::string , cmp3 >(lib2[0],tcmp) );
  }
 
 library_file.close();
 
 if ( cmps1.size() < 1 || cmps2.size() < 1 )
 {
  cout << "Could not find ligands in eFindSite pocket file" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( library.size() < 1 )
 {
  cout << "Screening library seems empty" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << "Number of target ligand clusters (Daylight/MACCS): " << cmps1.size() << "/" << cmps2.size() << endl << endl
      << "Screening library size: " << library.size() << endl  << endl;
 
      if ( smethod == "tst" )
  cout << "Scoring function: classical Tanimoto coeff for Daylight" << endl;
 else if ( smethod == "tsa" )
  cout << "Scoring function: average Tanimoto coeff for Daylight" << endl;
 else if ( smethod == "tsc" )
  cout << "Scoring function: continuous Tanimoto coeff for Daylight" << endl;
 else if ( smethod == "tmt" )
  cout << "Scoring function: classical Tanimoto coeff for MACCS" << endl;
 else if ( smethod == "tma" )
  cout << "Scoring function: average Tanimoto coeff for MACCS" << endl;
 else if ( smethod == "tmc" )
  cout << "Scoring function: continuous Tanimoto coeff for MACCS" << endl;
 else if ( smethod == "sum" )
  cout << "Scoring function: data fusion using sum rule" << endl;
 else if ( smethod == "max" )
  cout << "Scoring function: data fusion using max rule" << endl;
 else if ( smethod == "min" )
  cout << "Scoring function: data fusion using min rule" << endl;
 else if ( smethod == "svm" )
  cout << "Scoring function: support vector machines" << endl;
 
 cout << endl << "Running ligand-based virtual screening ... " << flush;
 
 time(&t_bench1);
 
 map< std::string , cmp3 >::iterator i1;
 
 multimap< double , std::string > score1;
 
 multimap< double , std::string > score2;
 multimap< double , std::string > score3;
 multimap< double , std::string > score4;
 multimap< double , std::string > score5;
 multimap< double , std::string > score6;
 multimap< double , std::string > score7;
 
 double as  = 0.0;
 double as2 = 0.0;
 double asn = 0.0;
 double asd = 0.0;
 
 double bs[6];
 double bs2[6];
 double bsn[6];
 double bsd[6];
 
 for ( int ii3 = 0; ii3 < 6; ii3++ )
 {
  bs[ii3] = 0.0;
  bs2[ii3] = 0.0;
  bsn[ii3] = 0.0;
  bsd[ii3] = 0.0;
 }
 
 for ( i1 = library.begin() ; i1 != library.end(); i1++ )
 {
  if ( smethod == "tst" || smethod == "svm" || smethod == "sum" || smethod == "max" || smethod == "min" )
  {
   list<cmp1>::iterator i2;
   
   ((*i1).second).sco_tst = 0.0;
   
   for ( i2 = cmps1.begin() ; i2 != cmps1.end(); i2++ )
    ((*i1).second).sco_tst += (*i2).weight * getTanimoto1024(((*i1).second).fpt1, (*i2).fingerprint);
  }
  
  if ( smethod == "tsa" || smethod == "svm" || smethod == "sum" || smethod == "max" || smethod == "min" )
  {
   list<cmp1>::iterator i2;
   
   ((*i1).second).sco_tsa = 0.0;
   
   for ( i2 = cmps1.begin() ; i2 != cmps1.end(); i2++ )
    ((*i1).second).sco_tsa += (*i2).weight * getTanimotoAve1024(((*i1).second).fpt1, (*i2).fingerprint);
  }
  
  if ( smethod == "tsc" || smethod == "svm" || smethod == "sum" || smethod == "max" || smethod == "min" )
   ((*i1).second).sco_tsc = getTanimotoCnt1024( ((*i1).second).fpt1, profile1 );
  
  if ( smethod == "tmt" || smethod == "svm" || smethod == "sum" || smethod == "max" || smethod == "min" )
  {
   list<cmp2>::iterator i2;
   
   ((*i1).second).sco_tmt = 0.0;
   
   for ( i2 = cmps2.begin() ; i2 != cmps2.end(); i2++ )
    ((*i1).second).sco_tmt += (*i2).weight * getTanimoto166(((*i1).second).fpt2, (*i2).fingerprint);
  }
  
  if ( smethod == "tma" || smethod == "svm" || smethod == "sum" || smethod == "max" || smethod == "min" )
  {
   list<cmp2>::iterator i2;
   
   ((*i1).second).sco_tma = 0.0;
   
   for ( i2 = cmps2.begin() ; i2 != cmps2.end(); i2++ )
    ((*i1).second).sco_tma += (*i2).weight * getTanimotoAve166(((*i1).second).fpt2, (*i2).fingerprint);
  }
  
  if ( smethod == "tmc" || smethod == "svm" || smethod == "sum" || smethod == "max" || smethod == "min" )
   ((*i1).second).sco_tmc = getTanimotoCnt166( ((*i1).second).fpt2, profile2 );
  
  if ( smethod == "svm" )
  {
   ((*i1).second).sco_mw   = 0.5 * pow( ( mprop[0][0] - ((*i1).second).mw   ) / mprop[0][1], 2.0) - log( 1.0 / ( mprop[0][1] * sqrt( 2.0 * PI ) ) );
   ((*i1).second).sco_logp = 0.5 * pow( ( mprop[1][0] - ((*i1).second).logp ) / mprop[1][1], 2.0) - log( 1.0 / ( mprop[1][1] * sqrt( 2.0 * PI ) ) );
   ((*i1).second).sco_psa  = 0.5 * pow( ( mprop[2][0] - ((*i1).second).psa  ) / mprop[2][1], 2.0) - log( 1.0 / ( mprop[2][1] * sqrt( 2.0 * PI ) ) );
   
   ((*i1).second).sco_hbd  = 0.5 * pow( ( mprop[3][0] - ((*i1).second).hbd  ) / mprop[3][1], 2.0) - log( 1.0 / ( mprop[3][1] * sqrt( 2.0 * PI ) ) );
   ((*i1).second).sco_hba  = 0.5 * pow( ( mprop[4][0] - ((*i1).second).hba  ) / mprop[4][1], 2.0) - log( 1.0 / ( mprop[4][1] * sqrt( 2.0 * PI ) ) );
   
   double ssvm1[MAXSV8];
   
   ssvm1[0]  = ((*i1).second).sco_tst;
   ssvm1[1]  = ((*i1).second).sco_tsa;
   ssvm1[2]  = ((*i1).second).sco_tsc;
   ssvm1[3]  = ((*i1).second).sco_tmt;
   ssvm1[4]  = ((*i1).second).sco_tma;
   ssvm1[5]  = ((*i1).second).sco_tmc;
   ssvm1[6]  = ((*i1).second).sco_mw;
   ssvm1[7]  = ((*i1).second).sco_logp;
   ssvm1[8]  = ((*i1).second).sco_psa;
   ssvm1[9]  = ((*i1).second).sco_hbd;
   ssvm1[10] = ((*i1).second).sco_hba;
   
   ((*i1).second).sco_svm  = model_svm->SVMpredict( 9, ssvm1 );
  }
  
       if ( smethod == "tst" )
  {
   as  += ((*i1).second).sco_tst;
   as2 += ((*i1).second).sco_tst * ((*i1).second).sco_tst;
   asn += 1.0;
   
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tst,(*i1).first) );
  }
  else if ( smethod == "tsa" )
  {
   as  += ((*i1).second).sco_tsa;
   as2 += ((*i1).second).sco_tsa * ((*i1).second).sco_tsa;
   asn += 1.0;
   
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tsa,(*i1).first) );
  }
  else if ( smethod == "tsc" )
  {
   as  += ((*i1).second).sco_tsc;
   as2 += ((*i1).second).sco_tsc * ((*i1).second).sco_tsc;
   asn += 1.0;
   
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tsc,(*i1).first) );
  }
  else if ( smethod == "tmt" )
  {
   as  += ((*i1).second).sco_tmt;
   as2 += ((*i1).second).sco_tmt * ((*i1).second).sco_tmt;
   asn += 1.0;
   
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tmt,(*i1).first) );
  }
  else if ( smethod == "tma" )
  {
   as  += ((*i1).second).sco_tma;
   as2 += ((*i1).second).sco_tma * ((*i1).second).sco_tma;
   asn += 1.0;
   
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tma,(*i1).first) );
  }
  else if ( smethod == "tmc" )
  {
   as  += ((*i1).second).sco_tmc;
   as2 += ((*i1).second).sco_tmc * ((*i1).second).sco_tmc;
   asn += 1.0;
   
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tmc,(*i1).first) );
  }
  else if ( smethod == "svm" )
  {
   as  += ((*i1).second).sco_svm;
   as2 += ((*i1).second).sco_svm * ((*i1).second).sco_svm;
   asn += 1.0;
   
   score1.insert( pair< double , std::string >(-1.0*((*i1).second).sco_svm,(*i1).first) );
  }
  else if ( smethod == "sum" || smethod == "max" || smethod == "min" )
  {
   score2.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tst,(*i1).first) );
   score3.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tsa,(*i1).first) );
   score4.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tsc,(*i1).first) );
   score5.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tmt,(*i1).first) );
   score6.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tma,(*i1).first) );
   score7.insert( pair< double , std::string >(-1.0*((*i1).second).sco_tmc,(*i1).first) );
   
   if ( smethod == "sum" )
   {
    bs[0]  += ((*i1).second).sco_tst;
    bs2[0] += ((*i1).second).sco_tst * ((*i1).second).sco_tst;
    bsn[0] += 1.0;
    
    bs[1]  += ((*i1).second).sco_tsa;
    bs2[1] += ((*i1).second).sco_tsa * ((*i1).second).sco_tsa;
    bsn[1] += 1.0;
    
    bs[2]  += ((*i1).second).sco_tsc;
    bs2[2] += ((*i1).second).sco_tsc * ((*i1).second).sco_tsc;
    bsn[2] += 1.0;
    
    bs[3]  += ((*i1).second).sco_tmt;
    bs2[3] += ((*i1).second).sco_tmt * ((*i1).second).sco_tmt;
    bsn[3] += 1.0;
    
    bs[4]  += ((*i1).second).sco_tma;
    bs2[4] += ((*i1).second).sco_tma * ((*i1).second).sco_tma;
    bsn[4] += 1.0;
    
    bs[5]  += ((*i1).second).sco_tmc;
    bs2[5] += ((*i1).second).sco_tmc * ((*i1).second).sco_tmc;
    bsn[5] += 1.0;
   }
  }
 }
 
 cmps1.clear();
 cmps2.clear();
 
 if ( smethod == "sum" || smethod == "max" || smethod == "min" )
 {
  multimap< double , std::string >::iterator i3;
  
  int i4;
  
  for ( i3 = score2.begin(), i4 = 1; i3 != score2.end(); i3++, i4++ )
   (library.find((*i3).second)->second).ran_tst = i4 / (double) library.size();
  
  for ( i3 = score3.begin(), i4 = 1; i3 != score3.end(); i3++, i4++ )
   (library.find((*i3).second)->second).ran_tsa = i4 / (double) library.size();
  
  for ( i3 = score4.begin(), i4 = 1; i3 != score4.end(); i3++, i4++ )
   (library.find((*i3).second)->second).ran_tsc = i4 / (double) library.size();
  
  for ( i3 = score5.begin(), i4 = 1; i3 != score5.end(); i3++, i4++ )
   (library.find((*i3).second)->second).ran_tmt = i4 / (double) library.size();
  
  for ( i3 = score6.begin(), i4 = 1; i3 != score6.end(); i3++, i4++ )
   (library.find((*i3).second)->second).ran_tma = i4 / (double) library.size();
  
  for ( i3 = score7.begin(), i4 = 1; i3 != score7.end(); i3++, i4++ )
   (library.find((*i3).second)->second).ran_tmc = i4 / (double) library.size();
  
  map< std::string , cmp3 >::iterator i5;
  
  for ( i5 = library.begin() ; i5 != library.end(); i5++ )
  {
   double tsc1 = 0.0;
   
   if ( smethod == "sum" )
    tsc1 = ( ((*i5).second).ran_tst + ((*i5).second).ran_tsa + ((*i5).second).ran_tsc + ((*i5).second).ran_tmt + ((*i5).second).ran_tma + ((*i5).second).ran_tmc ) / 6.0;
   else if ( smethod == "max" )
   {
    tsc1 = ((*i5).second).ran_tst;
    
    if ( ((*i5).second).ran_tsa < tsc1 ) tsc1 = ((*i5).second).ran_tsa;
    if ( ((*i5).second).ran_tsc < tsc1 ) tsc1 = ((*i5).second).ran_tsc;
    if ( ((*i5).second).ran_tmt < tsc1 ) tsc1 = ((*i5).second).ran_tmt;
    if ( ((*i5).second).ran_tma < tsc1 ) tsc1 = ((*i5).second).ran_tma;
    if ( ((*i5).second).ran_tmc < tsc1 ) tsc1 = ((*i5).second).ran_tmc;
   }
   else if ( smethod == "min" )
   {
    tsc1 = ((*i5).second).ran_tst;
    
    if ( ((*i5).second).ran_tsa > tsc1 ) tsc1 = ((*i5).second).ran_tsa;
    if ( ((*i5).second).ran_tsc > tsc1 ) tsc1 = ((*i5).second).ran_tsc;
    if ( ((*i5).second).ran_tmt > tsc1 ) tsc1 = ((*i5).second).ran_tmt;
    if ( ((*i5).second).ran_tma > tsc1 ) tsc1 = ((*i5).second).ran_tma;
    if ( ((*i5).second).ran_tmc > tsc1 ) tsc1 = ((*i5).second).ran_tmc;
   }
   
   as  += tsc1;
   as2 += tsc1 * tsc1;
   asn += 1.0;
   
   score1.insert( pair< double , std::string >(tsc1,(*i5).first) );
  }
 }
 
 as  /= asn;
 as2 /= asn;
 
 asd = sqrt( as2 - as * as );
 
 time(&t_bench2);
 
 cout << "done (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 multimap< double , std::string >::iterator i6;
 
 int i7;
 
 double rtop[6];
 
 ofstream outrank( (output_name).c_str() );
 
 for ( i6 = score1.begin(), i7 = 1; i6 != score1.end(); i6++, i7++ )
 {
  double score1t = 0.0;
  double score2t = 0.0;
  
  if ( smethod == "sum" || smethod == "max" || smethod == "min" )
  {
   score1t = (*i6).first;
   score2t = ( as - (*i6).first ) / asd;
   
   if ( smethod == "sum" and i7 == 1 )
    for ( int i8 = 0; i8 < 6; i8++ )
    {
     bs[i8]  /= bsn[i8];
     bs2[i8] /= bsn[i8];
     
     bsd[i8] = sqrt( bs2[i8] - bs[i8] * bs[i8] );
     
          if ( i8 == 0 )
      rtop[i8] = ( bs[i8] - (library.find((*i6).second)->second).ran_tst ) / bsd[i8];
     else if ( i8 == 1 )
      rtop[i8] = ( bs[i8] - (library.find((*i6).second)->second).ran_tsa ) / bsd[i8];
     else if ( i8 == 2 )
      rtop[i8] = ( bs[i8] - (library.find((*i6).second)->second).ran_tsc ) / bsd[i8];
     else if ( i8 == 3 )
      rtop[i8] = ( bs[i8] - (library.find((*i6).second)->second).ran_tmt ) / bsd[i8];
     else if ( i8 == 4 )
      rtop[i8] = ( bs[i8] - (library.find((*i6).second)->second).ran_tma ) / bsd[i8];
     else if ( i8 == 5 )
      rtop[i8] = ( bs[i8] - (library.find((*i6).second)->second).ran_tmc ) / bsd[i8];
    }
  }
  else
  {
   score1t = -1.0 * (*i6).first;
   score2t = ( -1.0 * (*i6).first - as ) / asd;
  }
  
  outrank << i7 << " " << (*i6).second << " " << fixed << setprecision(6) << score1t << " " << score2t << endl;
 }
 
 outrank.close();
 
 library.clear();
 
 score1.clear();
 score2.clear();
 score3.clear();
 score4.clear();
 score5.clear();
 score6.clear();
 score7.clear();
 
 std::string conf1 = "low";
 std::string conf2 = "low";
 
 if ( model_svm->SVMpredict( 7, rtop ) >= ctres1 )
  conf1 = "high";
 
 if ( model_svm->SVMpredict( 8, rtop ) >= ctres2 )
  conf2 = "high";
 
 cout << "Confidence top 1%:  " << conf1 << endl
      << "Confidence top 10%: " << conf2 << endl << endl;
 
 time(&t_end);
 
 printTime( difftime(t_end, t_start) );
 
 return 0;
}
