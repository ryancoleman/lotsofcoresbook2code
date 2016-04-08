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


#include "data.h"

using namespace std;


double getBkgFreq1( string resnam1 )
{
 // Release notes for UniProtKB/Swiss-Prot release 2011_09 - September 2011
 
      if ( resnam1 == "ALA" || resnam1 == "A" ) { return 0.0826; }
 else if ( resnam1 == "CYS" || resnam1 == "C" ) { return 0.0136; }
 else if ( resnam1 == "ASP" || resnam1 == "D" ) { return 0.0546; }
 else if ( resnam1 == "GLU" || resnam1 == "E" ) { return 0.0675; }
 else if ( resnam1 == "PHE" || resnam1 == "F" ) { return 0.0386; }
 else if ( resnam1 == "GLY" || resnam1 == "G" ) { return 0.0708; }
 else if ( resnam1 == "HIS" || resnam1 == "H" ) { return 0.0227; }
 else if ( resnam1 == "ILE" || resnam1 == "I" ) { return 0.0597; }
 else if ( resnam1 == "LYS" || resnam1 == "K" ) { return 0.0585; }
 else if ( resnam1 == "LEU" || resnam1 == "L" ) { return 0.0966; }
 else if ( resnam1 == "MET" || resnam1 == "M" ) { return 0.0242; }
 else if ( resnam1 == "ASN" || resnam1 == "N" ) { return 0.0406; }
 else if ( resnam1 == "PRO" || resnam1 == "P" ) { return 0.0469; }
 else if ( resnam1 == "GLN" || resnam1 == "Q" ) { return 0.0393; }
 else if ( resnam1 == "ARG" || resnam1 == "R" ) { return 0.0553; }
 else if ( resnam1 == "SER" || resnam1 == "S" ) { return 0.0655; }
 else if ( resnam1 == "THR" || resnam1 == "T" ) { return 0.0534; }
 else if ( resnam1 == "VAL" || resnam1 == "V" ) { return 0.0687; }
 else if ( resnam1 == "TRP" || resnam1 == "W" ) { return 0.0108; }
 else if ( resnam1 == "TYR" || resnam1 == "Y" ) { return 0.0292; }
 
 else
 {
  cout << "Unknown residue passed to getBackFreq1: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

std::string three2oneS( std::string resnam1 )
{
      if ( resnam1 == "ALA" ) { return "A"; }
 else if ( resnam1 == "CYS" ) { return "C"; }
 else if ( resnam1 == "ASP" ) { return "D"; }
 else if ( resnam1 == "GLU" ) { return "E"; }
 else if ( resnam1 == "PHE" ) { return "F"; }
 else if ( resnam1 == "GLY" ) { return "G"; }
 else if ( resnam1 == "HIS" ) { return "H"; }
 else if ( resnam1 == "ILE" ) { return "I"; }
 else if ( resnam1 == "LYS" ) { return "K"; }
 else if ( resnam1 == "LEU" ) { return "L"; }
 else if ( resnam1 == "MET" ) { return "M"; }
 else if ( resnam1 == "ASN" ) { return "N"; }
 else if ( resnam1 == "PRO" ) { return "P"; }
 else if ( resnam1 == "GLN" ) { return "Q"; }
 else if ( resnam1 == "ARG" ) { return "R"; }
 else if ( resnam1 == "SER" ) { return "S"; }
 else if ( resnam1 == "THR" ) { return "T"; }
 else if ( resnam1 == "VAL" ) { return "V"; }
 else if ( resnam1 == "TRP" ) { return "W"; }
 else if ( resnam1 == "TYR" ) { return "Y"; }
 
 else
 {
  cout << "Unknown residue passed to three2one: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

char three2oneC( std::string resnam1 )
{
      if ( resnam1 == "ALA" ) { return 'A'; }
 else if ( resnam1 == "CYS" ) { return 'C'; }
 else if ( resnam1 == "ASP" ) { return 'D'; }
 else if ( resnam1 == "GLU" ) { return 'E'; }
 else if ( resnam1 == "PHE" ) { return 'F'; }
 else if ( resnam1 == "GLY" ) { return 'G'; }
 else if ( resnam1 == "HIS" ) { return 'H'; }
 else if ( resnam1 == "ILE" ) { return 'I'; }
 else if ( resnam1 == "LYS" ) { return 'K'; }
 else if ( resnam1 == "LEU" ) { return 'L'; }
 else if ( resnam1 == "MET" ) { return 'M'; }
 else if ( resnam1 == "ASN" ) { return 'N'; }
 else if ( resnam1 == "PRO" ) { return 'P'; }
 else if ( resnam1 == "GLN" ) { return 'Q'; }
 else if ( resnam1 == "ARG" ) { return 'R'; }
 else if ( resnam1 == "SER" ) { return 'S'; }
 else if ( resnam1 == "THR" ) { return 'T'; }
 else if ( resnam1 == "VAL" ) { return 'V'; }
 else if ( resnam1 == "TRP" ) { return 'W'; }
 else if ( resnam1 == "TYR" ) { return 'Y'; }
 
 else
 {
  cout << "Unknown residue passed to three2one: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

std::string one2three( std::string resnam1 )
{
      if ( resnam1 == "A" ) { return "ALA"; }
 else if ( resnam1 == "C" ) { return "CYS"; }
 else if ( resnam1 == "D" ) { return "ASP"; }
 else if ( resnam1 == "E" ) { return "GLU"; }
 else if ( resnam1 == "F" ) { return "PHE"; }
 else if ( resnam1 == "G" ) { return "GLY"; }
 else if ( resnam1 == "H" ) { return "HIS"; }
 else if ( resnam1 == "I" ) { return "ILE"; }
 else if ( resnam1 == "K" ) { return "LYS"; }
 else if ( resnam1 == "L" ) { return "LEU"; }
 else if ( resnam1 == "M" ) { return "MET"; }
 else if ( resnam1 == "N" ) { return "ASN"; }
 else if ( resnam1 == "P" ) { return "PRO"; }
 else if ( resnam1 == "Q" ) { return "GLN"; }
 else if ( resnam1 == "R" ) { return "ARG"; }
 else if ( resnam1 == "S" ) { return "SER"; }
 else if ( resnam1 == "T" ) { return "THR"; }
 else if ( resnam1 == "V" ) { return "VAL"; }
 else if ( resnam1 == "W" ) { return "TRP"; }
 else if ( resnam1 == "Y" ) { return "TYR"; }
 
 else
 {
  cout << "Unknown residue passed to one2three: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

int one2num( char resnam1 )
{
      if ( resnam1 == 'A' ) { return  0; }
 else if ( resnam1 == 'C' ) { return  1; }
 else if ( resnam1 == 'D' ) { return  2; }
 else if ( resnam1 == 'E' ) { return  3; }
 else if ( resnam1 == 'F' ) { return  4; }
 else if ( resnam1 == 'G' ) { return  5; }
 else if ( resnam1 == 'H' ) { return  6; }
 else if ( resnam1 == 'I' ) { return  7; }
 else if ( resnam1 == 'K' ) { return  8; }
 else if ( resnam1 == 'L' ) { return  9; }
 else if ( resnam1 == 'M' ) { return 10; }
 else if ( resnam1 == 'N' ) { return 11; }
 else if ( resnam1 == 'P' ) { return 12; }
 else if ( resnam1 == 'Q' ) { return 13; }
 else if ( resnam1 == 'R' ) { return 14; }
 else if ( resnam1 == 'S' ) { return 15; }
 else if ( resnam1 == 'T' ) { return 16; }
 else if ( resnam1 == 'V' ) { return 17; }
 else if ( resnam1 == 'W' ) { return 18; }
 else if ( resnam1 == 'Y' ) { return 19; }
 
 else
 {
  cout << "Unknown residue passed to one2num: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

std::string num2one( int resnum1 )
{
      if ( resnum1 ==  0 ) { return "A"; }
 else if ( resnum1 ==  1 ) { return "C"; }
 else if ( resnum1 ==  2 ) { return "D"; }
 else if ( resnum1 ==  3 ) { return "E"; }
 else if ( resnum1 ==  4 ) { return "F"; }
 else if ( resnum1 ==  5 ) { return "G"; }
 else if ( resnum1 ==  6 ) { return "H"; }
 else if ( resnum1 ==  7 ) { return "I"; }
 else if ( resnum1 ==  8 ) { return "K"; }
 else if ( resnum1 ==  9 ) { return "L"; }
 else if ( resnum1 == 10 ) { return "M"; }
 else if ( resnum1 == 11 ) { return "N"; }
 else if ( resnum1 == 12 ) { return "P"; }
 else if ( resnum1 == 13 ) { return "Q"; }
 else if ( resnum1 == 14 ) { return "R"; }
 else if ( resnum1 == 15 ) { return "S"; }
 else if ( resnum1 == 16 ) { return "T"; }
 else if ( resnum1 == 17 ) { return "V"; }
 else if ( resnum1 == 18 ) { return "W"; }
 else if ( resnum1 == 19 ) { return "Y"; }
 
 else
 {
  cout << "Unknown residue passed to num2one: " << resnum1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

double one2plb( std::string resnam1 )
{
      if ( resnam1 == "A" ) { return 0.701; }
 else if ( resnam1 == "C" ) { return 1.650; }
 else if ( resnam1 == "D" ) { return 1.015; }
 else if ( resnam1 == "E" ) { return 0.956; }
 else if ( resnam1 == "F" ) { return 1.952; }
 else if ( resnam1 == "G" ) { return 0.788; }
 else if ( resnam1 == "H" ) { return 2.286; }
 else if ( resnam1 == "I" ) { return 1.006; }
 else if ( resnam1 == "K" ) { return 0.468; }
 else if ( resnam1 == "L" ) { return 1.045; }
 else if ( resnam1 == "M" ) { return 1.894; }
 else if ( resnam1 == "N" ) { return 0.811; }
 else if ( resnam1 == "P" ) { return 0.212; }
 else if ( resnam1 == "Q" ) { return 0.669; }
 else if ( resnam1 == "R" ) { return 0.916; }
 else if ( resnam1 == "S" ) { return 0.883; }
 else if ( resnam1 == "T" ) { return 0.730; }
 else if ( resnam1 == "V" ) { return 0.884; }
 else if ( resnam1 == "W" ) { return 3.084; }
 else if ( resnam1 == "Y" ) { return 1.672; }
 
 else
 {
  cout << "Unknown residue passed to one2plb: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}
