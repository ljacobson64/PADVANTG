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
  NpyArray scalar_flux_fwd_np = read_pickle("scalar_flux_fwd");
  NpyArray scalar_flux_adj_np = read_pickle("scalar_flux_adj");
  NpyArray sigma_t_pert_np = read_pickle("sigma_t_pert");
  NpyArray sigma_s_pert_np = read_pickle("sigma_s_pert");
  NpyArray material_map_np = read_pickle("material_map");
  NpyArray mesh_g_np = read_pickle("mesh_g");

  // Dimensions
  int n_mix = sigma_t_pert_np.shape[0];  // Number of mixed materials
  int nm = sigma_t_pert_np.shape[1];     // Number of pure materials
  int nz = scalar_flux_fwd_np.shape[0];  // Number of Z intervals
  int ny = scalar_flux_fwd_np.shape[1];  // Number of Y intervals
  int nx = scalar_flux_fwd_np.shape[2];  // Number of X intervals
  int ngf = scalar_flux_fwd_np.shape[3]; // Number of energy groups in flux
  int ngx = sigma_t_pert_np.shape[2];    // Number of energy groups in XS

  // Convert pickle data to Boost MultiArrays
  multi_array_ref<double, 4> scalar_flux_fwd{scalar_flux_fwd_np.data<double>(),
                                             extents[nz][ny][nx][ngf]};
  multi_array_ref<double, 4> scalar_flux_adj{scalar_flux_adj_np.data<double>(),
                                             extents[nz][ny][nx][ngf]};
  multi_array_ref<double, 3> sigma_t_pert{sigma_t_pert_np.data<double>(),
                                          extents[n_mix][nm][ngx]};
  multi_array_ref<double, 4> sigma_s_pert{sigma_s_pert_np.data<double>(),
                                          extents[n_mix][nm][ngx][ngx]};
  multi_array_ref<int, 3> material_map{material_map_np.data<int>(),
                                       extents[nz][ny][nx]};
  multi_array_ref<int, 1> mesh_g{mesh_g_np.data<int>(), extents[ngf]};

  // First energy group
  int g0 = mesh_g[0];

  // Allocate array for dR
  multi_array<double, 4> dR_scalar{extents[nm][nz][ny][nx]};

  // Calculate dR
  high_resolution_clock::time_point begin = high_resolution_clock::now();
  for (int im = 0; im < nm; im++) {             // Material index
    for (int iz = 0; iz < nz; iz++) {           // Z mesh index
      for (int iy = 0; iy < ny; iy++) {         // Y mesh index
        for (int ix = 0; ix < nx; ix++) {       // X mesh index
          int i_mix = material_map[iz][iy][ix]; // Mixed material index
          for (int igf = 0; igf < ngf; igf++) { // Energy index (in flux data)
            int igx = igf + g0; // Energy index (in cross section data)
            // Total component of dHphi
            double dHphi_t =
                sigma_t_pert[i_mix][im][igx] * scalar_flux_fwd[iz][iy][ix][igf];
            // Scattering component of dHphi
            double dHphi_s = 0.0;
            for (int jgf = 0; jgf < ngf; jgf++) {
              int jgx = jgf + g0;
              dHphi_s += -sigma_s_pert[i_mix][im][igx][jgx] *
                         scalar_flux_fwd[iz][iy][ix][jgf];
            }
            double dHphi = dHphi_t + dHphi_s;
            double dR_comp = -scalar_flux_adj[iz][iy][ix][igf] * dHphi;
            dR_scalar[im][iz][iy][ix] += dR_comp;
          }
        }
      }
    }
  }
  high_resolution_clock::time_point end = high_resolution_clock::now();
  double elapsed = duration_cast<duration<double>>(end - begin).count();
  cout << "Calculated dR for all materials in " << elapsed << " seconds"
       << endl;

  // Write pickles to file
  write_pickle("dR_scalar", dR_scalar);

  return 0;
}
