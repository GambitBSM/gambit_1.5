///
///  \author Rose Kudzman-Blais
///  \date 2017 May
///
///  *********************************************

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <hdf5.h>

using namespace std;

namespace Gambit {
  namespace ColliderBit {


    class Perf_Plot {
    private:

      const char* _outfilename;
      size_t _numvariables;   
      vector<const char*> _variables;
      vector<vector<double>> _values;

    public:

      ~Perf_Plot() {

      }


      Perf_Plot(string outFileName, vector<const char*>* varNames) {

        string temp = outFileName+".h5";
	_outfilename = temp.c_str();
	_variables = *varNames;
	_numvariables = _variables.size(); 
      }

      void fill(vector<double>* varValues) {
	_values.push_back(*varValues);
      }

      void createFile() {
	hid_t file;

        file = H5Fcreate(_outfilename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

        int nvalues = _values.size();

        for (size_t iVal=0;iVal<_numvariables;iVal++) { 

          hid_t dataset, dataspace;
          hsize_t dims[1];
          herr_t status;

          dims[0] = nvalues;
          dataspace = H5Screate_simple(1, dims, NULL);

          dataset = H5Dcreate2(file, _variables.at(iVal), H5T_NATIVE_DOUBLE, dataspace, H5S_ALL, H5S_ALL, H5P_DEFAULT); 

          double data[nvalues];
          for (size_t iVal2=0;iVal2<nvalues;iVal2++)data[iVal2]=_values.at(iVal2).at(iVal);

          status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

          H5Sclose(dataspace);
          H5Dclose(dataset);

        }  

        H5Fclose(file);	
	
      }

    };
  }
}
