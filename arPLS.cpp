#include "linalg.h"
#include "solvers.h"
#include "ap.h"
#include "statistics.h"
#include <math.h>
#include <iostream>
#include <fstream>

alglib::real_1d_array smooth(alglib::real_1d_array& y, double lambda, double ratio) {
  int n = y.length();

  // calculate H
  double c[] = {0, 0, 1, -2, 1, 0, 0};
  alglib::sparsematrix H;
  alglib::sparsecreatesksband(n, n, 5, H);
  // for each row
  for(int i=0;i<n;i++) {
    int rowOffset; int colOffset;
    if(i == 0) { rowOffset = 4; }
    else if(i == 1) { rowOffset = 3; }
    else if(i == n - 2) { rowOffset = 1; }
    else if(i == n - 1) { rowOffset = 0; }
    else { rowOffset = 2; }
    // for three bands
    for(int b = 0; b<3; b++) {
      colOffset = rowOffset - b;
      double s = alglib::vdotproduct(&c[rowOffset], &c[colOffset], 3);
      if(i+b < n) {
        alglib::sparseset(H, i, i+b, lambda * s);
        // std::cout << s << " ";
      }
    }
    //std::cout << std::endl;
  }

  alglib::real_1d_array w; w.setlength(n);
  for(int i = 0;i<n;i++) w[i] = 1.0;

  alglib::sparsematrix C;
  alglib::real_1d_array t; t.setlength(n);
  alglib::real_1d_array z; z.setlength(n);
  alglib::real_1d_array d; d.setlength(n);
  alglib::real_1d_array dn; dn.setlength(n);
  alglib::real_1d_array wt; wt.setlength(n);
  alglib::real_1d_array delta_w; delta_w.setlength(n);
  double w_norm, delta_w_norm;
  int nIters = 0;
  do {
    alglib::sparsecopybuf(H, C);
    for(int i = 0;i<n;i++) {
      alglib::sparseset(C, i, i, alglib::sparseget(C, i,i) + w[i]);
      t[i] = w[i] * y[i];
    }
    alglib::sparsesolverreport rep;
    alglib::sparsecholeskyskyline(C, n, true);
    alglib::sparsecholeskysolvesks(C, n, true, t, rep, z);
    alglib::vmove(&d[0], &z[0], n, -1.0);
    alglib::vadd(&d[0], &y[0], n);
    int num_neg = 0;
    for(int i=0;i<n;i++) {
      if(d[i] < 0) dn[num_neg++] = d[i];
    }
    double m = alglib::samplemean(dn, num_neg);
    double s = std::sqrt(alglib::samplevariance(dn, num_neg));

    for(int i=0;i<n;i++) {
      wt[i] = 1.0 / (1 + std::exp(2 * (d[i] - (2*s-m))/s));
    }

    w_norm = alglib::vdotproduct(&w[0], &w[0], n);
    alglib::vmove(&delta_w[0], &w[0], n);
    alglib::vadd(&delta_w[0], &wt[0], n, -1.0);
    delta_w_norm = alglib::vdotproduct(&delta_w[0], &delta_w[0], n);
    alglib::vmove(&w[0], &wt[0], n);

    std::cout << "|w - wt| / |w|=" << delta_w_norm / w_norm << std::endl;
  } while(delta_w_norm / w_norm >= ratio && nIters++ < 50);

  std::cout << std::endl;

  return d;
}

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

