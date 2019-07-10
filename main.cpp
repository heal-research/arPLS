#include "arPLS.hpp"
#include "linalg.h"
#include "solvers.h"
#include "ap.h"
#include "statistics.h"
#include <math.h>
#include <iostream>
#include <fstream>

int main(int argc, char** argv) {
  try {
    double ratio = 0.001;
    double lambda = 1;

    if(argc < 2) {
      std::cout << "Usage: arPLS <spectrum> [<lambda> [<ratio>]]" << std::endl;
      std::cout << "  lambda: smoothness parameter (larger than zero, default=1000, range: 10^2 .. 10^8)" << std::endl;
      std::cout << "  ratio: stopping parameter. Algorithm stops when the change of relative weights becomes small (|w - wt|/|w| < ratio). (default=0.001)" << std::endl;
      return 1;
    }

    if(argc >= 3) lambda = atof(argv[2]);
    if(argc == 4) ratio = atof(argv[3]);
    
    alglib::real_2d_array spectrum;
    alglib::read_csv(argv[1], ' ', alglib::CSV_DEFAULT, spectrum);
    int n = spectrum.rows();
    
    alglib::real_1d_array y; y.setlength(n);
    alglib::vmove(&y[0], 1, &spectrum[0][1], spectrum.getstride(), n);

    // arPLS smoothing
    alglib::real_1d_array ySmoothed = smooth(y, lambda, ratio);

    // write
    std::ofstream outFile;
    std::string prefix = std::string("filtered_");
    std::string filename =  std::string(argv[1]);
    outFile.open(prefix+filename);
    for(int i=0;i<n;i++) {
      outFile << spectrum[i][0] << " " << ySmoothed[i] << std::endl;
    }
    outFile.close();
    
    return 0;
  } catch(alglib::ap_error & e) {
    std::cout << e.msg << std::endl;
    return 2;
  }
}
