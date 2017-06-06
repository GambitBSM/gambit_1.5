//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
//   
//   \file Slha files for specific model points 
//   
//   \author Rose Kudzman-Blais
//   
//   *********************************************

#include <iostream>
#include <sstream>
#include <string>
#include <slhaea.h>

using namespace std;
using namespace SLHAea;

      int main() {

	double M0 = 70.;
	double M12 = 350.;
	double TanBeta = 10.;
	double SignMu = 1.;
	double A0 = 0.;
  	
	Coll out1;
  	string block = "MINPAR";
  	out1[block][""] << "BLOCK" << block;
  	out1[block][""] <<      1  << M0      << "# M0";
  	out1[block][""] <<      2  << M12     << "# M12";
  	out1[block][""] <<      3  << TanBeta << "# TanBeta";
  	out1[block][""] <<      4  << SignMu  << "# SignMu";
  	out1[block][""] <<      5  << A0      << "# A0";
  
	cout << out1;

      }
