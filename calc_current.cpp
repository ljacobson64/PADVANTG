#include <chrono>

#include "boost/multi_array.hpp"
#include "cnpy.h"

#define PI 3.14159265358979323846
#define FOURPI (4.0 * PI)

using boost::extents;
using boost::multi_array;
using boost::multi_array_ref;
using cnpy::npy_load;
using cnpy::npy_save;
using cnpy::NpyArray;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

// Read pickle from file
NpyArray read_pickle(string fname) {
  cout << "Reading pickles/" << fname << ".npy from file...";
  cout.flush();
  high_resolution_clock::time_point begin = high_resolution_clock::now();
  NpyArray pickle = npy_load("pickles/" + fname + ".npy");
  high_resolution_clock::time_point end = high_resolution_clock::now();
  double elapsed = duration_cast<duration<double>>(end - begin).count();
  cout << " took " << elapsed << " seconds" << endl;
  cout.flush();
  return pickle;
}

// Write pickle to file
template <typename T, size_t NDIMS>
void write_pickle(string fname, multi_array<T, NDIMS> pickle) {
  cout << "Writing pickles/" << fname << ".npy to file...";
  cout.flush();
  high_resolution_clock::time_point begin = high_resolution_clock::now();
  const long unsigned int *shape_arr = pickle.shape();
  vector<long unsigned int> shape(shape_arr, shape_arr + NDIMS);
  npy_save("pickles/" + fname + ".npy", pickle.origin(), shape, "w");
  high_resolution_clock::time_point end = high_resolution_clock::now();
  double elapsed = duration_cast<duration<double>>(end - begin).count();
  cout << " took " << elapsed << " seconds" << endl;
  cout.flush();
}

int main() {
  // Load pickle data
  NpyArray angular_flux_fwd_np = read_pickle("angular_flux_fwd");
  NpyArray angular_flux_adj_np = read_pickle("angular_flux_adj");
  NpyArray quadrature_weights_np = read_pickle("quadrature_weights");
  NpyArray angles_np = read_pickle("angles");

  // Dimensions
  int nz = angular_flux_fwd_np.shape[0];  // Number of Z intervals
  int ny = angular_flux_fwd_np.shape[1];  // Number of Y intervals
  int nx = angular_flux_fwd_np.shape[2];  // Number of X intervals
  int na = angular_flux_fwd_np.shape[3];  // Number of angles
  int ngf = angular_flux_fwd_np.shape[4]; // Number of energy groups in flux

  // Convert pickle data to Boost MultiArrays
  multi_array_ref<double, 5> angular_flux_fwd{
      angular_flux_fwd_np.data<double>(), extents[nz][ny][nx][na][ngf]};
  multi_array_ref<double, 5> angular_flux_adj{
      angular_flux_adj_np.data<double>(), extents[nz][ny][nx][na][ngf]};
  multi_array_ref<double, 1> quadrature_weights{
      quadrature_weights_np.data<double>(), extents[na]};
  multi_array_ref<double, 2> angles{angles_np.data<double>(), extents[na][3]};

  // Allocate arrays for current
  multi_array<double, 5> current_fwd{extents[nz][ny][nx][ngf][3]};
  multi_array<double, 5> current_adj{extents[nz][ny][nx][ngf][3]};
  multi_array<double, 5> current_con{extents[nz][ny][nx][ngf][3]};

  // Calculate current
  high_resolution_clock::time_point begin = high_resolution_clock::now();
  for (int iz = 0; iz < nz; iz++) {                    // Z mesh index
    for (int iy = 0; iy < ny; iy++) {                  // Y mesh index
      for (int ix = 0; ix < nx; ix++) {                // X mesh index
        for (int ia = 0; ia < na; ia++) {              // Angle index
          double qw = quadrature_weights[ia] / FOURPI; // Quadrature weight
          for (int igf = 0; igf < ngf; igf++) {        // Energy index
            double ffwd = angular_flux_fwd[iz][iy][ix][ia][igf];
            double fadj = angular_flux_adj[iz][iy][ix][ia][igf];
            for (int id = 0; id < 3; id++) { // Dimension index
              current_fwd[iz][iy][ix][igf][id] += ffwd * qw * angles[ia][id];
              current_adj[iz][iy][ix][igf][id] += fadj * qw * angles[ia][id];
            }
          }
        }
      }
    }
  }
  for (int iz = 0; iz < nz; iz++)
    for (int iy = 0; iy < ny; iy++)
      for (int ix = 0; ix < nx; ix++)
        for (int igf = 0; igf < ngf; igf++)
          for (int id = 0; id < 3; id++)
            current_con[iz][iy][ix][igf][id] =
                current_fwd[iz][iy][ix][igf][id] *
                current_adj[iz][iy][ix][igf][id];
  high_resolution_clock::time_point end = high_resolution_clock::now();
  double elapsed = duration_cast<duration<double>>(end - begin).count();
  cout << "Calculated current in " << elapsed << " seconds" << endl;

  // Write pickles to file
  write_pickle("current_fwd", current_fwd);
  write_pickle("current_adj", current_adj);
  write_pickle("current_con", current_con);

  return 0;
}
