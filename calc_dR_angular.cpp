#include "boost/multi_array.hpp"
#include "cnpy.h"

#define PI 3.14159265358979323846

using boost::extents;
using boost::multi_array;
using boost::multi_array_ref;
using cnpy::npy_load;
using cnpy::npy_save;
using cnpy::NpyArray;
using std::cout;
using std::endl;
using std::string;

NpyArray read_pickle(std::string fname) {
  cout << "Reading pickles/" << fname << ".npy from file" << endl;
  NpyArray pickle = npy_load("pickles/" + fname + ".npy");
  return pickle;
}

int main() {
  // Load pickle data
  NpyArray flux_fwd_np = read_pickle("angular_flux_fwd");
  NpyArray flux_adj_np = read_pickle("angular_flux_adj");
  NpyArray quad_weights_np = read_pickle("quadrature_weights");
  NpyArray sigma_t_pert_np = read_pickle("sigma_t_pert");
  NpyArray sigma_s_pert_np = read_pickle("sigma_s_pert");
  NpyArray material_map_np = read_pickle("material_map");
  NpyArray mesh_g_np = read_pickle("mesh_g");

  // Dimensions
  int n_mix = sigma_t_pert_np.shape[0]; // Number of mixed materials
  int nm = sigma_t_pert_np.shape[1];    // Number of pure materials
  int nz = flux_fwd_np.shape[0];        // Number of Z intervals
  int ny = flux_fwd_np.shape[1];        // Number of Y intervals
  int nx = flux_fwd_np.shape[2];        // Number of X intervals
  int na = flux_fwd_np.shape[3];        // Number of angles
  int ngf = flux_fwd_np.shape[4];       // Number of energy groups in flux
  int ngx = sigma_t_pert_np.shape[2];   // Number of energy groups in XS

  // Convert pickle data to Boost MultiArrays
  multi_array_ref<double, 5> flux_fwd{flux_fwd_np.data<double>(),
                                      extents[nz][ny][nx][na][ngf]};
  multi_array_ref<double, 5> flux_adj{flux_adj_np.data<double>(),
                                      extents[nz][ny][nx][na][ngf]};
  multi_array_ref<double, 1> quad_weights{quad_weights_np.data<double>(),
                                          extents[na]};
  multi_array_ref<double, 3> sigma_t_pert{sigma_t_pert_np.data<double>(),
                                          extents[n_mix][nm][ngx]};
  multi_array_ref<double, 4> sigma_s_pert{sigma_s_pert_np.data<double>(),
                                          extents[n_mix][nm][ngx][ngx]};
  multi_array_ref<long long int, 3> material_map{
      material_map_np.data<long long int>(), extents[nz][ny][nx]};
  multi_array_ref<long long int, 1> mesh_g{mesh_g_np.data<long long int>(),
                                           extents[ngf]};

  // First energy group
  int g0 = mesh_g[0];

  // Allocate array for dR
  multi_array<double, 4> dR{boost::extents[nm][nz][ny][nx]};

  // Calculate dR
  for (int im = 0; im < nm; im++) { // Material index
    cout << "Calculating dR for material " << im + 1 << " of " << nm << endl;
    for (int iz = 0; iz < nz; iz++) {             // Z mesh index
      for (int iy = 0; iy < ny; iy++) {           // Y mesh index
        for (int ix = 0; ix < nx; ix++) {         // X mesh index
          int i_mix = material_map[iz][iy][ix];   // Mixed material index
          for (int ia = 0; ia < na; ia++) {       // Angle index
            double qw = quad_weights[ia];         // Quadrature weight
            for (int igf = 0; igf < ngf; igf++) { // Energy index (in flux data)
              int igx = igf + g0; // Energy index (in cross section data)
              // Total component of dHphi
              double dHphi_t =
                  sigma_t_pert[i_mix][im][igx] * flux_fwd[iz][iy][ix][ia][igf];
              // Scattering component of dHphi
              double dHphi_s = 0.0;
              for (int igf2 = 0; igf2 < ngf; igf2++) {
                int igx2 = igf2 + g0;
                dHphi_s += -sigma_s_pert[i_mix][im][igx][igx2] *
                           flux_fwd[iz][iy][ix][ia][igf2];
              }
              // dR
              double dHphi = dHphi_t + dHphi_s;
              double dR_comp =
                  -flux_adj[iz][iy][ix][ia][igf] * dHphi * qw / (4.0 * PI);
              dR[im][iz][iy][ix] += dR_comp;
            }
          }
        }
      }
    }
  }

  // Write dR to file
  string fname = "pickles/dR_angular.npy";
  cout << "Writing " << fname << " to file" << endl;
  npy_save(fname, &dR[0][0][0][0], {nm, nz, ny, nx}, "w");

  return 0;
}
