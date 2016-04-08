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


#include "list.h"

using namespace std;

void getList( string templates_name, string cluster_def, list<string> &template_list, map<string,double> &template_prob1, map<string,double> &template_prob2 )
{
 multimap<string,string> cluster_map;
 
 multimap<string,double> template_tmp1;
 multimap<string,double> template_tmp2;
 
 string line1;
 
 ifstream clusters_file( cluster_def.c_str() );
 
 if ( !clusters_file.is_open() )  { cout << "Cannot open " << cluster_def << endl; exit(EXIT_FAILURE); }
 
 while (getline(clusters_file,line1))
 {
  string tpl1;
  int    tpl3;
  
  tpl1 = line1.substr(0,15);
  
  string::size_type tpl6 = 0;
  
  bool tpl7 = true;
  
  while( tpl7 )
  {
   tpl6 = tpl1.find(" ");
   if( tpl6 != string::npos )
    tpl1.erase( tpl6, 1 );
   else
    tpl7 = false;
  }
  
  tpl3 = atoi(line1.substr(15,5).c_str());
  
  for ( int tpl4 = 0; tpl4 < tpl3; tpl4++ )
  {
   string tpl5;
   
   tpl5 = line1.substr(21+6*tpl4,5);
   
   cluster_map.insert( pair<string,string>(tpl1,tpl5) );
  }
 }
 
 clusters_file.close();
 
 multimap<int,string,std::greater<int> > template_fun;
 
 string line2;
 
 ifstream template_file( templates_name.c_str() );
 
 if ( !template_file.is_open() )  { cout << "Cannot open " << templates_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(template_file,line2))
  if ( line2.length() > 71 )
  {
   if ( line2.compare(0, 1, "#") != 0 )
   {
    string tpl1;
    int    tpl2;
    double tpl3;
    double tpl4;
    
    tpl1 = line2.substr(0,15);
    
    string::size_type tpl6 = 0;
    
    bool tpl7 = true;
    
    while( tpl7 )
    {
     tpl6 = tpl1.find(" ");
     if( tpl6 != string::npos )
      tpl1.erase( tpl6, 1 );
     else
      tpl7 = false;
    }
    
    tpl2 = atoi(line2.substr(57,6).c_str());
    
    tpl3 = atof(line2.substr(55,8).c_str());
    tpl4 = atof(line2.substr(64,8).c_str());
    
    multimap<string,string>::iterator it1;
    pair<multimap<string,string>::iterator,multimap<string,string>::iterator> it2;
    
    it2 = cluster_map.equal_range(tpl1);
    
    for ( it1 = it2.first; it1 != it2.second; it1++ )
    {
     template_fun.insert( pair<int,string>(tpl2,(*it1).second) );
     
     template_tmp1.insert( pair<string,double>((*it1).second, tpl3) );
     template_tmp2.insert( pair<string,double>((*it1).second, tpl4) );
    }
   }
  }
  else if ( line2.length() )
  {
   if ( line2.compare(0, 1, "#") != 0 )
   {
    string tpl1;
    int    tpl2;
    double tpl3;
    double tpl4;
    
    tpl1 = line2;
    
    string::size_type tpl6 = 0;
    
    bool tpl7 = true;
    
    while( tpl7 )
    {
     tpl6 = tpl1.find(" ");
     if( tpl6 != string::npos )
      tpl1.erase( tpl6, 1 );
     else
      tpl7 = false;
    }
    
    tpl2 = rand();
    
    tpl3 = 1.0;
    tpl4 = 1.0;
    
    multimap<string,string>::iterator it1;
    pair<multimap<string,string>::iterator,multimap<string,string>::iterator> it2;
    
    it2 = cluster_map.equal_range(tpl1);
    
    for ( it1 = it2.first; it1 != it2.second; it1++ )
    {
     template_fun.insert( pair<int,string>(tpl2,(*it1).second) );
     
     template_tmp1.insert( pair<string,double>((*it1).second, tpl3) );
     template_tmp2.insert( pair<string,double>((*it1).second, tpl4) );
    }
   }
  }
 
 template_file.close();
 
 multimap<int,string>::iterator it3;
 
 for ( it3 = template_fun.begin() ; it3 != template_fun.end(); it3++ )
 {
  bool w1 = false;
  
  if ( template_list.size() )
  {
   list<string>::iterator it4;
   
   it4 = find(template_list.begin(), template_list.end(), (*it3).second);
   
   if ( it4 == template_list.end() )
    w1 = true;
  }
  else
   w1 = true;
  
  if ( w1 )
  {
   double prob1 = 0.0;
   double prob2 = 0.0;
   
   multimap<string,double>::iterator it5;
   pair<multimap<string,double>::iterator,multimap<string,double>::iterator> it6;
   
   it6 = template_tmp1.equal_range( (*it3).second );
   
   for ( it5 = it6.first; it5 != it6.second; it5++ )
    if ( (*it5).second > prob1 )
     prob1 = (*it5).second;
   
   it6 = template_tmp2.equal_range( (*it3).second );
   
   for ( it5 = it6.first; it5 != it6.second; it5++ )
    if ( (*it5).second > prob2 )
     prob2 = (*it5).second;
   
   template_prob1.insert( pair<string,double>((*it3).second, prob1) );
   template_prob2.insert( pair<string,double>((*it3).second, prob2) );
   
   template_list.push_back( (*it3).second );
  }
 }
}
