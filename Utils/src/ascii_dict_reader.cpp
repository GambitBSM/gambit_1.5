//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Simple reader for ASCII tables
///
///  *********************************************
///
///  Authors (add name and date if you modify):
//
///  \author Janina Renk
///          <janina.renk@fysik.su.se>
///  \date Oct 2018
///
///  *********************************************

#include "gambit/Utils/ascii_dict_reader.hpp"

namespace Gambit
{


  int ASCIIdictReader::read(std::string filename)
  { 
    std::ifstream in(filename.c_str(), std::ios::binary);
    if (in.fail())
    {
      std::cout << "ERROR. Failed loading: " << filename << std::endl;
      // TODO: Throw proper IO error
      exit(-1);
    }
    std::string line;
    while(std::getline(in, line))
    {
      if (line[0] == '#') continue;  // Ignore comments lines, starting with "#"
      std::stringstream ss(line);

      double tmp;
      std::string key;
      std::vector<double> data_tmp;

      ss >> key;
      while(ss >> tmp)
      {
        data_tmp.push_back(tmp);
      }
  

      if (std::find(keys.begin(), keys.end(), key) != keys.end())
      {
        duplicate = true;
      }

      //else {
      dict[key] = data_tmp;
      keys.push_back(key); 
      //}
      
    }
    in.close();
    return 0;
  }



}
