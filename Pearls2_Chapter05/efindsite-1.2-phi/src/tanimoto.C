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


#include "tanimoto.h"

using namespace std;

double getTanimoto1024( bitset<MAXSMI> &tc_fpt1, bitset<MAXSMI> &tc_fpt2 )
{
 bitset<MAXSMI> xftp1 = tc_fpt1;
 bitset<MAXSMI> xftp2 = tc_fpt2;
 
 int v1 = (xftp1^=xftp2).count();
 int v2 = (xftp1&=xftp2).count();
 int v3 = (xftp1|=xftp2).count();
 
 double tc1 = 0.0;
 
 if ( ( v1 - v2 + v3 ) > 0 )
  tc1 = ( ( (double) ( v3 - v2 ) ) / ( (double) ( v1 - v2 + v3 ) ) );
 
 return tc1;
}

double getTanimotoAve1024( bitset<MAXSMI> &tc_fpt1, bitset<MAXSMI> &tc_fpt2 )
{
 bitset<MAXSMI> xftp1 = tc_fpt1;
 bitset<MAXSMI> xftp2 = tc_fpt2;
 
 bitset<MAXSMI> xftp3 = tc_fpt1; xftp3 = xftp3.flip();
 bitset<MAXSMI> xftp4 = tc_fpt2; xftp4 = xftp4.flip();
 
 int v1r = (xftp1^=xftp2).count();
 int v2r = (xftp1&=xftp2).count();
 int v3r = (xftp1|=xftp2).count();
 
 int v1f = (xftp3^=xftp4).count();
 int v2f = (xftp3&=xftp4).count();
 int v3f = (xftp3|=xftp4).count();
 
 double tc1 = 0.0;
 
 if ( ( v1r - v2r + v3r ) > 0 )
  tc1 = ( 0.5 * ( ( ( (double) ( v3r - v2r ) ) / ( (double) ( v1r - v2r + v3r ) ) ) + ( ( (double) ( v3f - v2f ) ) / ( (double) ( v1f - v2f + v3f ) ) ) ) );
 
 return tc1;
}

double getTanimotoCnt1024( bitset<MAXSMI> &tc_fpt1, double tc_fpt2[] )
{
 bitset<MAXSMI> xftp1 = tc_fpt1;
 
 std::string xftp2 = tc_fpt1.to_string<char,char_traits<char>,allocator<char> >();
 
 string::iterator xftp3;
 
 int xftp4;
 
 double v1 = 0;
 double v2 = 0;
 double v3 = 0;
 
 for ( xftp3 = xftp2.begin(), xftp4 = 0; xftp3 < xftp2.end(); xftp3++, xftp4++ )
 {
  if ( (*xftp3) == '1' )
  {
   v3 += tc_fpt2[xftp4];
   
   v2++;
  }
  
  v1 += ( tc_fpt2[xftp4] * tc_fpt2[xftp4] );
 }
 
 double tc1 = 0.0;
 
 if ( ( v1 + v2 - v3 ) > 0 )
  tc1 = ( v3 / ( v1 + v2 - v3 ) );
 
 return tc1;
}

double getTanimoto166( bitset<MAXMAC> &tc_fpt1, bitset<MAXMAC> &tc_fpt2 )
{
 bitset<MAXMAC> xftp1 = tc_fpt1;
 bitset<MAXMAC> xftp2 = tc_fpt2;
 
 int v1 = (xftp1^=xftp2).count();
 int v2 = (xftp1&=xftp2).count();
 int v3 = (xftp1|=xftp2).count();
 
 double tc1 = 0.0;
 
 if ( ( v1 - v2 + v3 ) > 0 )
  tc1 = ( ( (double) ( v3 - v2 ) ) / ( (double) ( v1 - v2 + v3 ) ) );
 
 return tc1;
}

double getTanimotoAve166( bitset<MAXMAC> &tc_fpt1, bitset<MAXMAC> &tc_fpt2 )
{
 bitset<MAXMAC> xftp1 = tc_fpt1;
 bitset<MAXMAC> xftp2 = tc_fpt2;
 
 bitset<MAXMAC> xftp3 = tc_fpt1; xftp3 = xftp3.flip();
 bitset<MAXMAC> xftp4 = tc_fpt2; xftp4 = xftp4.flip();
 
 int v1r = (xftp1^=xftp2).count();
 int v2r = (xftp1&=xftp2).count();
 int v3r = (xftp1|=xftp2).count();
 
 int v1f = (xftp3^=xftp4).count();
 int v2f = (xftp3&=xftp4).count();
 int v3f = (xftp3|=xftp4).count();
 
 double tc1 = 0.0;
 
 if ( ( v1r - v2r + v3r ) > 0 )
  tc1 = ( 0.5 * ( ( ( (double) ( v3r - v2r ) ) / ( (double) ( v1r - v2r + v3r ) ) ) + ( ( (double) ( v3f - v2f ) ) / ( (double) ( v1f - v2f + v3f ) ) ) ) );
 
 return tc1;
}

double getTanimotoCnt166( bitset<MAXMAC> &tc_fpt1, double tc_fpt2[] )
{
 bitset<MAXMAC> xftp1 = tc_fpt1;
 
 std::string xftp2 = tc_fpt1.to_string<char,char_traits<char>,allocator<char> >();
 
 string::iterator xftp3;
 
 int xftp4;
 
 double v1 = 0;
 double v2 = 0;
 double v3 = 0;
 
 for ( xftp3 = xftp2.begin(), xftp4 = 0; xftp3 < xftp2.end(); xftp3++, xftp4++ )
 {
  if ( (*xftp3) == '1' )
  {
   v3 += tc_fpt2[xftp4];
   
   v2++;
  }
  
  v1 += ( tc_fpt2[xftp4] * tc_fpt2[xftp4] );
 }
 
 double tc1 = 0.0;
 
 if ( ( v1 + v2 - v3 ) > 0 )
  tc1 = ( v3 / ( v1 + v2 - v3 ) );
 
 return tc1;
}

