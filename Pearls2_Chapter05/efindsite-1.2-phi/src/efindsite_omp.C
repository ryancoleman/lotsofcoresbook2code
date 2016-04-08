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


#include "efindsite_omp.h"

using namespace std;

int main(int argc, char *argv[])
{
 time_t t_start, t_end, t_bench1, t_bench2;
 time(&t_start);
 
 double       cut_tmscore   = 0.40;
 double       cut_seqid     = 1.00;
 double       cut_clustdis  = 8.00;
 unsigned int cut_templates = MAXTPL;
 double       cut_binrest   = 0.25;
 int          cut_binresn   = 1;
 double       cut_clustlig  = 0.70;
 std::string  met_clustlig  = "T";
 
#ifdef __linux
 
 std::string  met_clustdis  = "P";
 
#else
 
 std::string  met_clustdis  = "L";
 
#endif
 
 cout << "------------------------------------------------------------" << endl
      << "                         efindsite" << endl
      << "                        version 1.2" << endl
      << "              ligand binding pocket prediction" << endl << endl
      << "       report bugs and issues to michal@brylinski.org" << endl
      << "------------------------------------------------------------" << endl << endl;
 
 if ( argc < 7 )
 {
  cout << " efindsite -s <target structure in PDB format>" << endl
       << "           -t <templates detected by eThread>" << endl
       << "           -i <secondary structure profile by psipred>" << endl
       << "           -e <sequence profile>" << endl
       << "           -o <output filename>" << endl << endl
       << " additional options:" << endl << endl
       << "           -l <auxiliary ligands in SDF format>" << endl
       << "           -b <sequence identity threshold (default 1.0)>" << endl
       << "           -m <TMscore threshold (default 0.4)>" << endl
       << "           -x <max number of templates (default " << MAXTPL << ")>" << endl
       << "           -r <binding residue threshold (default 0.25)>" << endl
       << "           -n <min # of binding residues (default 1)>" << endl
       << "           -d <pocket clustering cutoff (default 8.0)>" << endl

#ifdef __linux

       << "           -g <pocket clustering method (Linux only, default P)>" << endl
       << "               P - affinity propagation" << endl
       << "               L - average linkage" << endl

#endif

       << "           -f <fingerprint clustering cutoff (default 0.7)>" << endl
       << "           -c <fingerprint clustering method (default T)>" << endl
       << "               T - classical Tanimoto coeff" << endl
       << "               A - average Tanimoto coeff" << endl << endl;
  
  exit(EXIT_SUCCESS);
 }
 
 string target_name;
 bool target_opt = false;
 
 string psipred_name;
 bool psipred_opt = false;
 
 string sequence_name;
 bool sequence_opt = false;
 
 string templates_name;
 bool templates_opt = false;
 
 string output_name;
 bool output_opt = false;
 
 string cmps_name;
 bool cmps_opt = false;
 
 for ( int i = 0; i < argc; i++ )
 {
  if ( !strcmp(argv[i],"-s") && i < argc ) { target_name    = string(argv[i+1]); target_opt    = true; }
  if ( !strcmp(argv[i],"-i") && i < argc ) { psipred_name   = string(argv[i+1]); psipred_opt   = true; }
  if ( !strcmp(argv[i],"-e") && i < argc ) { sequence_name  = string(argv[i+1]); sequence_opt  = true; }
  if ( !strcmp(argv[i],"-t") && i < argc ) { templates_name = string(argv[i+1]); templates_opt = true; }
  if ( !strcmp(argv[i],"-o") && i < argc ) { output_name    = string(argv[i+1]); output_opt    = true; }
  if ( !strcmp(argv[i],"-l") && i < argc ) { cmps_name      = string(argv[i+1]); cmps_opt      = true; }
  if ( !strcmp(argv[i],"-b") && i < argc ) { cut_seqid      = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-m") && i < argc ) { cut_tmscore    = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-x") && i < argc ) { cut_templates  = atoi(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-r") && i < argc ) { cut_binrest    = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-n") && i < argc ) { cut_binresn    = atoi(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-f") && i < argc ) { cut_clustlig   = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-c") && i < argc ) { met_clustlig   = string(argv[i+1]);                       }
  if ( !strcmp(argv[i],"-d") && i < argc ) { cut_clustdis   = atof(argv[i+1]);                         }
  
#ifdef __linux
  
  if ( !strcmp(argv[i],"-g") && i < argc ) { met_clustdis   = string(argv[i+1]);                       }
  
#endif
  
 }
 
 char * path1;
 
 path1 = getenv("EF_LIB"); if ( path1==NULL ) { cout << "EF_LIB is not set" << endl; exit(EXIT_FAILURE); }
 
 path1 = getenv("EF_MAP"); if ( path1==NULL ) { cout << "EF_MAP is not set" << endl; exit(EXIT_FAILURE); }
 
 path1 = getenv("EF_MOD"); if ( path1==NULL ) { cout << "EF_MOD is not set" << endl; exit(EXIT_FAILURE); }
 
#ifdef __linux
 
 path1 = getenv("EF_APL"); if ( path1==NULL ) { cout << "EF_APL is not set" << endl; exit(EXIT_FAILURE); }
 
#endif
 
 string lib_path;
 lib_path = getenv("EF_LIB");
 
 string cluster_def;
 cluster_def = getenv("EF_MAP");
 
 string model_path;
 model_path = getenv("EF_MOD");
 
#ifdef __linux
 
 string ap_lib;
 ap_lib = getenv("EF_APL");
 
#endif
 
 ifstream f01( (model_path+"/compositionSVM.model").c_str() );
 ifstream f02( (model_path+"/profileSVM.model").c_str() );
 ifstream f03( (model_path+"/residueSVM.model").c_str() );
 ifstream f04( (model_path+"/ligandSVM.model").c_str() );
 ifstream f05( (model_path+"/rank7SVM.model").c_str() );
 ifstream f06( (model_path+"/rank8SVM.model").c_str() );
 ifstream f07( (model_path+"/compositionSVM.scale").c_str() );
 ifstream f08( (model_path+"/profileSVM.scale").c_str() );
 ifstream f09( (model_path+"/residueSVM.scale").c_str() );
 ifstream f10( (model_path+"/ligandSVM.scale").c_str() );
 ifstream f11( (model_path+"/rank7SVM.scale").c_str() );
 ifstream f12( (model_path+"/rank8SVM.scale").c_str() );
 
 if ( !f01 || !f02 || !f03 || !f04 || !f05 || !f06 || !f07 || !f08 || !f09 || !f10 || !f11 || !f12 )
 {
  cout << "Could not find SVM models in " << model_path << endl;
  exit(EXIT_FAILURE);
 }
 
#ifdef __linux
 
 ifstream ap1( (ap_lib).c_str() );
 
 if ( !ap1 )
 {
  cout << "Could not find Affinity Propagation library: " << ap_lib << endl;
  exit(EXIT_FAILURE);
 }
 
#endif
 
 if ( !target_opt )
 {
  cout << "Provide target structure in PDB format" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !psipred_opt )
 {
  cout << "Provide secondary structure profile by psipred" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !sequence_opt )
 {
  cout << "Provide sequence profile" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !templates_opt )
 {
  cout << "Provide templates detected by eThread" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !output_opt )
 {
  cout << "Provide output filename" << endl;
  exit(EXIT_FAILURE);
 }
 else
 {
  ofstream outprot( (output_name+".templates.pdb").c_str() );
  outprot.close();
  
  ofstream outali( (output_name+".alignments.dat").c_str() );
  outali.close();
  
  ofstream outlig( (output_name+".ligands.sdf").c_str() );
  outlig.close();
  
  ofstream outpkt1( (output_name+".pockets.pdb").c_str() );
  outpkt1.close();
  
  ofstream outpkt2( (output_name+".pockets.dat").c_str() );
  outpkt2.close();
 }
 
 if ( cut_seqid < 1.0 )
  cout << "!!! Benchmarking mode activated with max sid of " << cut_seqid << " !!!" << endl << endl;
 
 if ( cut_tmscore < 0.4 )
  cout << "!!! TMscore of " << cut_tmscore << " is below the statistical significance threshold !!!" << endl << endl;
 
 if ( cut_binresn < 1 )
 {
  cout << "!!! Min # of binding residues must be >0, setting to 1 !!!" << endl << endl;
  
  cut_binresn = 1;
 }
 
 if ( cut_binrest < ( 1 / 1e6 ) )
 {
  cout << "!!! Threshold for binding residues must be >0, setting to 0.18 !!!" << endl << endl;
  
  cut_binrest = 0.18;
 }
 else if ( cut_binrest > 1 )
 {
  cout << "!!! Threshold for binding residues must be <=1, setting to 0.18 !!!" << endl << endl;
  
  cut_binrest = 0.18;
 }
 
 if ( cut_clustlig < ( 1 / 1e6 ) )
 {
  cout << "!!! Fingerprint clustering cutoff must be >0, setting to 0.1 !!!" << endl << endl;
  
  cut_clustlig = 0.1;
 }
 else if ( cut_clustlig > 1 )
 {
  cout << "!!! Fingerprint clustering cutoff must be <=1, setting to 1.0 !!!" << endl << endl;
  
  cut_clustlig = 1.0;
 }
 
 if ( met_clustlig != "T" && met_clustlig != "A" )
 {
  cout << "!!! Fingerprint clustering method must be either T or A, setting to T !!!" << endl << endl;
  
  met_clustlig = "T";
 }
 
#ifdef __linux
 
 if ( met_clustdis != "P" && met_clustdis != "L" )
 {
  cout << "!!! Pocket clustering method must be either P or L, setting to P !!!" << endl << endl;
  
  met_clustdis = "P";
 }
 
#endif
 
 if ( cut_templates > (int) MAXTPL )
 {
  cout << "!!! Max number of templates exceeded, setting to " << MAXTPL << " !!!" << endl << endl;
  
  cut_templates = MAXTPL;
 }
 
 if ( cut_templates < 1 )
 {
  cout << "!!! Max number of templates must be >0, setting to 1 !!!" << endl << endl;
  
  cut_templates = 1;
 }
 
#ifdef __linux
 
 cout << "Checking Affinity Propagation library ... " << flush;
 
 void *dlh = NULL;
 
 if ( !( dlh = dlopen( (ap_lib).c_str(), RTLD_LAZY ) ) )
 {
  cout << dlerror() << endl;
  exit(EXIT_FAILURE);
 }
 
 char *error;
 
 if ( ( error = dlerror() ) != NULL)
 {
  cout << error << endl;
  exit(EXIT_FAILURE);
 }
 
 dlclose(dlh);
 
 cout << "looks good" << endl << endl;
 
#endif
 
 Target * target;
 
 target = new Target( 0, 0 );
 
 if ( target->loadTarget(target_name) )
 {
  cout << "Cannot read target structure" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( target->loadPsipred(psipred_name) )
 {
  cout << "Cannot read psipred profile" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( target->loadSequence(sequence_name) )
 {
  cout << "Cannot read sequence profile" << endl;
  exit(EXIT_FAILURE);
 }
 
 Cmps * compounds;
 
 compounds = new Cmps( 0 );
 
 if ( cmps_opt )
 {
  if ( compounds->loadCompounds(cmps_name) )
  {
   cout << "Cannot read auxiliary compounds" << endl;
   exit(EXIT_FAILURE);
  }
 }
 
 list<string> template_list;
 
 map<string,double> template_prob1;
 map<string,double> template_prob2;
 
 getList( templates_name, cluster_def, template_list, template_prob1, template_prob2 );
 
 cout << "eFindSite library: " << lib_path <<  endl
      << "eFindSite mapping: " << cluster_def << endl << endl
 
      << "Number of ligand-bound templates: " << setw(5) << template_list.size() << endl << endl;
 
 cout << "Template filtering ... " << flush;
 
 time(&t_bench1);
 
 multimap< int, Template *, greater<int> > template_set;
 
 list<string>::iterator it1;
 
 for ( it1 = template_list.begin(); it1 != template_list.end(); it1++ )
  if ( template_set.size() < cut_templates )
  {
   Template * template_tmp = new Template( 0, 0, 0, template_prob1.find(*it1)->second, template_prob2.find(*it1)->second );
   
   bool load1 = template_tmp->loadTemplate( lib_path+"/data/"+(*it1).substr(1,2)+"/"+(*it1), template_list );
   
   if ( !load1 )
   {
    double sid1 = template_tmp->alignNW( target->getProteinSequence() );
    
    if ( sid1 <= cut_seqid )
     template_set.insert( std::pair< int,Template * >( template_tmp->getProteinResiduesTotal(), template_tmp ) );
   }
  }
 
 int t_offset[MAXTPL];
 int n_offset = 0;
 int s_offset = 0;
 
 std::multimap< int, Template *, greater<int> >::iterator tpl1;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); tpl1++ )
 {
  t_offset[n_offset++] = s_offset;
  
  s_offset += (*tpl1).first;
 }
 
 time(&t_bench2);
 
 cout << ". " << template_set.size() << " templates survived (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Estimating data size ... " << flush;
 
 time(&t_bench1);
 
 int n_tar = target->getProteinResiduesTotal();
 int n_off = s_offset;
 int n_set = template_set.size();
 
 int  t_len1;
 int *t_len2 = new int [n_set];
 
 char *t_seq1 = new char [n_tar];
 char *t_seq2 = new char [n_off];
 
 int *t_res1 = new int [n_tar];
 int *t_res2 = new int [n_off];
 
 double *t_xyz1 = new double [n_tar * 3];
 double *t_xyz2 = new double [n_off * 3];
 
 int    *t_sco1 = new int    [n_set];
 double *t_sco2 = new double [n_set];
 double *t_sco3 = new double [n_set];
 double *t_sco4 = new double [n_set];
 
 int *t_alig = new int [n_tar * n_set];
 
 double *t_rmat = new double [n_set * 12];
 
 t_len1 = target->getProteinResiduesTotal();
 
 char t_seqs[MAXPRO];
 
 strcpy(t_seqs, (target->getProteinSequence()).c_str());
 
 for ( int t_i = 0; t_i < target->getProteinResiduesTotal(); t_i++ )
  t_seq1[t_i] = t_seqs[t_i];
  
 for ( int t_i = 0; t_i < t_len1; t_i++ )
  t_res1[t_i] = t_i + 1;
 
 target->getProteinCoords1D(t_xyz1);
 
 int i_offset = 0;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); tpl1++ )
 {
  t_len2[i_offset] = ((*tpl1).second)->getProteinResiduesTotal();
  
  char t_seqt[MAXPRO];
  
  strcpy(t_seqt, (((*tpl1).second)->getProteinSequence()).c_str());
  
  for ( int t_i = 0; t_i < ((*tpl1).second)->getProteinResiduesTotal(); t_i++ )
   t_seq2[t_i+t_offset[i_offset]] = t_seqt[t_i];
  
  for ( int t_i = 0; t_i < ((*tpl1).second)->getProteinResiduesTotal(); t_i++ )
   t_res2[t_i+t_offset[i_offset]] = t_i + 1;
  
  double *t_xyzt = new double [((*tpl1).second)->getProteinResiduesTotal() * 3];
  
  ((*tpl1).second)->getProteinCoords1D(t_xyzt);
  
  for ( int t_i = 0; t_i < ((*tpl1).second)->getProteinResiduesTotal() * 3; t_i++ )
   t_xyz2[t_i+t_offset[i_offset]*3] = t_xyzt[t_i];
  
  delete [] t_xyzt;
  
  i_offset++;
 }
 
 int mem_len = sizeof(t_len1) + n_set * sizeof(t_len2);
 int mem_seq = n_tar * sizeof(t_seq1) + n_off * sizeof(t_seq2);
 int mem_res = n_tar * sizeof(t_res1) + n_off * sizeof(t_res2);
 int mem_xyz = n_tar * 3 * sizeof(t_xyz1) + n_off * 3 * sizeof(t_xyz2);
 int mem_sco = n_set * sizeof(t_sco1) + n_set * sizeof(t_sco2) + n_set * sizeof(t_sco3) + n_set * sizeof(t_sco4);
 int mem_ali = n_tar * n_set * sizeof(t_alig);
 int mem_mat = n_set * 12 * sizeof(t_rmat);
 
 time(&t_bench2);
 
 cout << setprecision(1) << ( (double) ( mem_len + mem_seq + mem_res + mem_xyz + mem_sco + mem_ali + mem_mat ) ) / 1e6 << "M (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Calculating alignments ... " << flush;
 
 time(&t_bench1);
 
 int tm_i;
 
 #pragma omp parallel for schedule(dynamic) private(tm_i)
 for ( tm_i = 0; tm_i < n_offset; tm_i++ )
 {
  int    o_sco1;
  double o_sco2, o_sco3, o_sco4;
  
  int o_len1 = t_len1;
  int o_len2 = t_len2[tm_i];
  
  char o_seq1[MAXPRO];
  char o_seq2[MAXPRO];
  
  int o_res1[MAXPRO];
  int o_res2[MAXPRO];
  
  double o_xyz1[MAXPRO][3];
  double o_xyz2[MAXPRO][3];
  
  for ( int t_i = 0; t_i < o_len1; t_i++ )
  {
   o_seq1[t_i] = t_seq1[t_i];
   
   o_res1[t_i] = t_res1[t_i];
   
   for ( int t_j = 0; t_j < 3; t_j++ )
    o_xyz1[t_i][t_j] = t_xyz1[t_i*3+t_j];
  }
  
  for ( int t_i = 0; t_i < o_len2; t_i++ )
  {
   o_seq2[t_i] = t_seq2[t_offset[tm_i]+t_i];
   
   o_res2[t_i] = t_res2[t_offset[tm_i]+t_i];
   
   for ( int t_j = 0; t_j < 3; t_j++ )
    o_xyz2[t_i][t_j] = t_xyz2[t_offset[tm_i]*3+t_i*3+t_j];
  }
  
  int o_alig[MAXPRO];
  double o_t[3];
  double o_u[3][3];
  
  frtmalign_( &o_sco1, &o_sco2, &o_sco3, &o_sco4, &o_len2, &o_len1, &o_seq2, &o_seq1, &o_alig, &o_res2, &o_res1, &o_xyz2, &o_xyz1, &o_t, &o_u, &o_len1 );
  
  t_sco1[tm_i] = o_sco1;
  t_sco2[tm_i] = o_sco2;
  t_sco3[tm_i] = o_sco3;
  t_sco4[tm_i] = 0.0;
  
  for ( int t_i = 0; t_i < o_len1; t_i++ )
   t_alig[tm_i*o_len1+t_i] = o_alig[t_i];
  
  for ( int t_i = 0; t_i < 3; t_i++ )
  {
   t_rmat[tm_i*12+t_i] = o_t[t_i];
   
   for ( int t_j = 0; t_j < 3; t_j++ )
    t_rmat[tm_i*12+3+3*t_i+t_j] = o_u[t_i][t_j];
  }
  
  for ( int t_i = 0; t_i < o_len1; t_i++ )
   if ( o_alig[t_i] > -1 )
    if ( o_seq1[t_i] == o_seq2[o_alig[t_i]] )
     t_sco4[tm_i]++;
  
  t_sco4[tm_i] = t_sco4[tm_i] / ( (double) t_sco1[tm_i] );
 }
 
 int j_offset = 0;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); tpl1++ )
 {
  ((*tpl1).second)->setProteinLengthTM( t_sco1[j_offset] );
  ((*tpl1).second)->setProteinRMSD( t_sco2[j_offset] );
  ((*tpl1).second)->setProteinTMscore( t_sco3[j_offset] );
  ((*tpl1).second)->setProteinSeqID2( t_sco4[j_offset] );
  
  int p_alig[MAXPRO];
  
  for ( int t_i = 0; t_i < t_len1; t_i++ )
   p_alig[t_i] = t_alig[j_offset*t_len1+t_i];
  
  ((*tpl1).second)->setTMalignment(p_alig, t_len1);
  
  double p_t[3];
  double p_u[3][3];
  
  for ( int t_i = 0; t_i < 3; t_i++ )
  {
   p_t[t_i] = t_rmat[j_offset*12+t_i];
   
   for ( int t_j = 0; t_j < 3; t_j++ )
    p_u[t_i][t_j] = t_rmat[j_offset*12+3+3*t_i+t_j];
  }
  
  ((*tpl1).second)->setMatrix(p_t, p_u);
  
  j_offset++;
 }
 
 delete [] t_len2;
 delete [] t_seq1;
 delete [] t_seq2;
 delete [] t_res1;
 delete [] t_res2;
 delete [] t_xyz1;
 delete [] t_xyz2;
 delete [] t_sco1;
 delete [] t_sco2;
 delete [] t_sco3;
 delete [] t_sco4;
 delete [] t_alig;
 delete [] t_rmat;
 
 list<string> template_list_filtered;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); )
 {
  std::multimap< int, Template *, greater<int> >::iterator tpl6 = tpl1++;
  
  if ( ((*tpl6).second)->getProteinTMscore() < cut_tmscore )
   template_set.erase(tpl6);
  else
   template_list_filtered.push_back( ((*tpl6).second)->getProteinID() );
 }
 
 int ltot = 0;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); tpl1++ )
 {
  ((*tpl1).second)->purgeAlignments(template_list_filtered);
  
  ((*tpl1).second)->calculateContacts();
  
  ltot += ((*tpl1).second)->getLigandsTotal();
 }
 
 if ( template_set.empty() )
 {
  cout << "no templates survived" << endl << endl;
  
  time(&t_end);
  
  printTime( difftime(t_end, t_start) );
  
  exit(EXIT_SUCCESS);
 }
 
 time(&t_bench2);
 
 cout << template_set.size() << "/" << ltot << " templates/ligands survived (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Detecting pockets ... " << flush;
 
 time(&t_bench1);
 
 double *sim1 = new double[ltot*ltot];
 
 std::multimap< int, Template *, greater<int> >::iterator tpl2;
 std::multimap< int, Template *, greater<int> >::iterator tpl3;
 
 int nl1 = 0;
 int nl2 = 0;
 
 for ( tpl2 = template_set.begin(); tpl2 != template_set.end(); tpl2++ )
  for ( int il1 = 0; il1 < ((*tpl2).second)->getLigandsTotal(); il1++ )
  {
   for ( tpl3 = template_set.begin(); tpl3 != template_set.end(); tpl3++ )
    for ( int il2 = 0; il2 < ((*tpl3).second)->getLigandsTotal(); il2++ )
    {
     sim1[nl1*ltot+nl2] = getDistance( 1, (*tpl2).second, il1, (*tpl3).second, il2 );
     
     nl2++;
    }
   
   nl2 = 0;
   
   nl1++;
  }
 
 int * clu1 = new int [nl1];
 
 int clu2;
 
#ifdef __linux
 
 if ( met_clustdis == "P" && template_set.size() > 4 )
  clu2 = cluster_ap( sim1, clu1, nl1, cut_clustdis, ap_lib);
 
 else
  clu2 = cluster_avelink( sim1, clu1, nl1, cut_clustdis, "min" );
 
#else
 
 clu2 = cluster_avelink( sim1, clu1, nl1, cut_clustdis, "min" );
 
#endif
 
 delete [] sim1;
 
 if ( clu2 < 1 )
 {
  cout << "no pockets found" << endl << endl;
  
  template_set.clear();
  
  time(&t_end);
  
  printTime( difftime(t_end, t_start) );
  
  exit(EXIT_SUCCESS);
 }
 
 clu2 = refine_pockets(template_set, clu2, clu1, nl1, cut_clustdis);
 
 time(&t_bench2);
 
 cout << clu2 << " pockets found (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 int nl3 = 0;
 
 std::multimap< int, Template *, greater<int> >::iterator tpl4;
 
 for ( tpl4 = template_set.begin(); tpl4 != template_set.end(); tpl4++ )
  for ( int il1 = 0; il1 < ((*tpl4).second)->getLigandsTotal(); il1++ )
   ((*tpl4).second)->setPocketNumber(il1, clu1[nl3++]);
 
 list< Pocket * > pocket_set;
 
 for ( int clu3 = 0; clu3 < clu2; clu3++ )
 {
  Pocket * pocket_tmp = new Pocket( clu3 );
  
  std::multimap< int, Template *, greater<int> >::iterator tpl5;
  
  for ( tpl5 = template_set.begin(); tpl5 != template_set.end(); tpl5++ )
   for ( int il1 = 0; il1 < ((*tpl5).second)->getLigandsTotal(); il1++ )
    if ( ((*tpl5).second)->getPocketNumber(il1) == clu3 )
     pocket_tmp->addTemplate((*tpl5).second);
  
  if ( pocket_tmp->getProteinsTotal() > 0 && pocket_tmp->getLigandsTotal() > 0 )
   pocket_set.push_back( pocket_tmp );
  else
   delete pocket_tmp;
 }
 
 delete [] clu1;
 
 cout << "Loading SVM models " << flush;
 
 ModelSVM * model_svm;
 
 model_svm = new ModelSVM( false, false, false, false, false, false, false, false, false );
 
 time(&t_bench1);
 
 model_svm->loadModel( 1, model_path+"/compositionSVM.model" ); cout << '.' << flush;
 model_svm->loadModel( 2, model_path+"/profileSVM.model" ); cout << '.' << flush;
 model_svm->loadModel( 3, model_path+"/residueSVM.model" ); cout << '.' << flush;
 model_svm->loadModel( 4, model_path+"/ligandSVM.model" ); cout << '.' << flush;
 model_svm->loadModel( 5, model_path+"/rank7SVM.model" ); cout << '.' << flush;
 model_svm->loadModel( 6, model_path+"/rank8SVM.model" ); cout << '.' << flush;
 
 model_svm->loadScale( 1, model_path+"/compositionSVM.scale" ); cout << '.' << flush;
 model_svm->loadScale( 2, model_path+"/profileSVM.scale" ); cout << '.' << flush;
 model_svm->loadScale( 3, model_path+"/residueSVM.scale" ); cout << '.' << flush;
 model_svm->loadScale( 4, model_path+"/ligandSVM.scale" ); cout << '.' << flush;
 model_svm->loadScale( 5, model_path+"/rank7SVM.scale" ); cout << '.' << flush;
 model_svm->loadScale( 6, model_path+"/rank8SVM.scale" ); cout << '.' << flush;
 
 time(&t_bench2);
 
 if ( target->compositionSVM(model_svm) )
 {
  cout << "SVM for aa composition failed" << endl;
  exit(EXIT_FAILURE);
 }
 
 cout << " done (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Predicting binding residues ... " << flush;
 
 time(&t_bench1);
 
 list< Pocket * > pocket_set_filtered;
 
 list< Pocket * >::iterator ipkt1;
 
 for ( ipkt1 = pocket_set.begin(); ipkt1 != pocket_set.end(); ipkt1++ )
 {
  (*ipkt1)->calculatePocketCenter();
  
  int bres1 = (*ipkt1)->calculateBindingResidues( target, model_svm, cut_binrest );
  
  if ( bres1 >= cut_binresn )
   pocket_set_filtered.push_back( *ipkt1 );
 }
 
 pocket_set.clear();
 
 if ( pocket_set_filtered.empty() )
 {
  cout << "no pockets survived" << endl << endl;
  
  template_set.clear();
  
  pocket_set_filtered.clear();
  
  time(&t_end);
  
  printTime( difftime(t_end, t_start) );
  
  exit(EXIT_SUCCESS);
 }
 
 time(&t_bench2);
 
 cout << pocket_set_filtered.size() << " pockets survived (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 double fra1 = 0.0;
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
  fra1 += (double) (*ipkt1)->getLigandsTotal();
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
 {
  double fra2 = 0.0;
  
  if ( fra1 > 0.0 )
   fra2 = ( (double) (*ipkt1)->getLigandsTotal() ) / fra1;
  
  (*ipkt1)->setPocketFraction( fra2 );
 }
 
 cout << "Constructing ligand fingerprints ... " << flush;
 
 time(&t_bench1);
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
 {
  (*ipkt1)->calculateFingerprintsSMILES( cut_clustlig, met_clustlig );
  (*ipkt1)->calculateFingerprintsMACCS( cut_clustlig, met_clustlig );
  
  if ( cmps_opt )
   (*ipkt1)->calculateCmpsScores( compounds, model_svm );
 }
 
 time(&t_bench2);
 
 cout << "done (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Ranking pockets ... " << flush;
 
 time(&t_bench1);
 
 double rank1 = 0.0;
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
 {
  double rank2 = (*ipkt1)->calculateConfidence( cmps_opt, model_svm );
  
  if ( rank2 > rank1 )
   rank1 = rank2;
 }
 
 time(&t_bench2);
 
 cout << "top-ranked pocket has a confidence index of " << fixed << setprecision(1) << rank1 * 100 << "% (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 multimap<double,Pocket *> pocket_map_sorted;
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
  pocket_map_sorted.insert( pair<double,Pocket *>(-1.0*(*ipkt1)->getConfidence(),*ipkt1) );
 
 list< Pocket * > pocket_set_sorted;
 
 multimap<double,Pocket *>::iterator ipkt3;
 
 for ( ipkt3 = pocket_map_sorted.begin() ; ipkt3 != pocket_map_sorted.end(); ipkt3++ )
  pocket_set_sorted.push_back( (*ipkt3).second );
 
 pocket_set_filtered.clear();
 
 map<string,bool> chk1;
 map<string,bool> chk2;
 
 int ipkt2 = 1;
 
 for ( ipkt1 = pocket_set_sorted.begin(); ipkt1 != pocket_set_sorted.end(); ipkt1++ )
 {
  (*ipkt1)->setCenter( cut_binrest, cut_clustdis );
  
  (*ipkt1)->dumpProteinsAlignments( output_name, chk1, target );
  
  (*ipkt1)->dumpPocket( output_name, target, cut_binrest, ipkt2 );
  
  (*ipkt1)->dumpLigands( output_name, chk2, ipkt2 );
  
  ipkt2++;
 }
 
 template_set.clear();
 
 pocket_set_sorted.clear();
 
 time(&t_end);
 
 printTime( difftime(t_end, t_start) );
 
 return 0;
}
