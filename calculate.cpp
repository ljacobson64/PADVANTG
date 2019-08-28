#include <chrono>
#include <cstdlib>

#include "boost/multi_array.hpp"
#include "cnpy.h"

#define PI 3.14159265358979323846
#define FOURPI 4.0 * PI

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
  // Variables for timing
  high_resolution_clock::time_point begin, end;
  double elapsed;

  // Load pickle data
  NpyArray mesh_x_np = read_pickle("mesh_x");
  NpyArray mesh_y_np = read_pickle("mesh_y");
  NpyArray mesh_z_np = read_pickle("mesh_z");
  NpyArray mesh_g_np = read_pickle("mesh_g");
  NpyArray angles_np = read_pickle("angles");
  NpyArray quadrature_weights_np = read_pickle("quadrature_weights");
  NpyArray source_indices_np = read_pickle("source_indices");
  NpyArray source_strengths_np = read_pickle("source_strengths");
  NpyArray source_spectrum_np = read_pickle("source_spectrum");
  NpyArray response_indices_np = read_pickle("response_indices");
  NpyArray response_strengths_np = read_pickle("response_strengths");
  NpyArray response_spectrum_np = read_pickle("response_spectrum");
  NpyArray mix_table_np = read_pickle("mix_table");
  NpyArray material_map_np = read_pickle("material_map");
  NpyArray sigma_t_np = read_pickle("sigma_t");
  NpyArray sigma_s_np = read_pickle("sigma_s");
  NpyArray angular_flux_fwd_np = read_pickle("angular_flux_fwd");
  NpyArray angular_flux_adj_np = read_pickle("angular_flux_adj");

  // Dimensions
  int n_mix = mix_table_np.shape[0];        // Number of mixed materials
  int nm = mix_table_np.shape[1];           // Number of pure materials
  int nz = angular_flux_fwd_np.shape[0];    // Number of Z intervals
  int ny = angular_flux_fwd_np.shape[1];    // Number of Y intervals
  int nx = angular_flux_fwd_np.shape[2];    // Number of X intervals
  int na = angular_flux_fwd_np.shape[3];    // Number of angles
  int ngf = angular_flux_fwd_np.shape[4];   // Number of energy groups in flux
  int ngx = sigma_t_np.shape[1];            // Number of energy groups in XS
  int n_src = source_indices_np.shape[0];   // Number of source indices
  int n_res = response_indices_np.shape[0]; // Number of response indices

  // Convert pickle data to Boost MultiArrays
  multi_array_ref<double, 1> mesh_x{mesh_x_np.data<double>(), extents[nx]};
  multi_array_ref<double, 1> mesh_y{mesh_y_np.data<double>(), extents[ny]};
  multi_array_ref<double, 1> mesh_z{mesh_z_np.data<double>(), extents[nz]};
  multi_array_ref<int, 1> mesh_g{mesh_g_np.data<int>(), extents[ngf]};
  multi_array_ref<double, 2> angles{angles_np.data<double>(), extents[na][3]};
  multi_array_ref<double, 1> quadrature_weights{
      quadrature_weights_np.data<double>(), extents[na]};
  multi_array_ref<int, 1> source_indices{source_indices_np.data<int>(),
                                         extents[n_src]};
  multi_array_ref<double, 1> source_strengths{
      source_strengths_np.data<double>(), extents[n_src]};
  multi_array_ref<double, 1> source_spectrum{source_spectrum_np.data<double>(),
                                             extents[ngx]};
  multi_array_ref<int, 1> response_indices{response_indices_np.data<int>(),
                                           extents[n_res]};
  multi_array_ref<double, 1> response_strengths{
      response_strengths_np.data<double>(), extents[n_res]};
  multi_array_ref<double, 1> response_spectrum{
      response_spectrum_np.data<double>(), extents[ngx]};
  multi_array_ref<double, 2> mix_table{mix_table_np.data<double>(),
                                       extents[n_mix][nm]};
  multi_array_ref<int, 3> material_map{material_map_np.data<int>(),
                                       extents[nz][ny][nx]};
  multi_array_ref<double, 2> sigma_t{sigma_t_np.data<double>(),
                                     extents[nm][ngx]};
  multi_array_ref<double, 3> sigma_s{sigma_s_np.data<double>(),
                                     extents[nm][ngx][ngx]};
  multi_array_ref<double, 5> angular_flux_fwd{
      angular_flux_fwd_np.data<double>(), extents[nz][ny][nx][na][ngf]};
  multi_array_ref<double, 5> angular_flux_adj{
      angular_flux_adj_np.data<double>(), extents[nz][ny][nx][na][ngf]};

  // First energy group in flux
  int g0 = mesh_g[0];

  // Calculate reverse angle map
  cout << "Calculating reverse angle map...";
  cout.flush();
  begin = high_resolution_clock::now();
  multi_array<int, 1> reverse_angle_map{extents[na]};
  for (int ia = 0; ia < na; ia++) {
    double ix = angles[ia][0];
    double iy = angles[ia][1];
    double iz = angles[ia][2];
    for (int ja = 0; ja < na; ja++) {
      double jx = angles[ja][0];
      double jy = angles[ja][1];
      double jz = angles[ja][2];
      if (ix == -jx and iy == -jy and iz == -jz) {
        reverse_angle_map[ia] = ja;
        cout << ia << " " << ja << endl;
        break;
      }
    }
  }
  end = high_resolution_clock::now();
  elapsed = duration_cast<duration<double>>(end - begin).count();
  cout << " took " << elapsed << " seconds" << endl;
  cout.flush();

  // Calculate source and response
  cout << "Calculating source and response...";
  cout.flush();
  begin = high_resolution_clock::now();
  multi_array<double, 4> source{extents[nz][ny][nx][ngf]};
  multi_array<double, 4> response{extents[nz][ny][nx][ngf]};
  for (int i = 0; i < n_src; i++) {
    int ind = source_indices[i];
    int iz = ind / (ny * nx);
    int iy = (ind / nx) % ny;
    int ix = ind % nx;
    double strength = source_strengths[i];
    for (int igf = 0; igf < ngf; igf++)
      source[iz][iy][ix][igf] = strength * source_spectrum[g0 + igf];
  }
  for (int i = 0; i < n_res; i++) {
    int ind = response_indices[i];
    int iz = ind / (ny * nx);
    int iy = (ind / nx) % ny;
    int ix = ind % nx;
    double strength = response_strengths[i];
    for (int igf = 0; igf < ngf; igf++)
      response[iz][iy][ix][igf] = strength * response_spectrum[g0 + igf];
  }
  end = high_resolution_clock::now();
  elapsed = duration_cast<duration<double>>(end - begin).count();
  cout << " took " << elapsed << " seconds" << endl;
  cout.flush();

  // Calculate scalar flux
  cout << "Calculating scalar flux...";
  cout.flush();
  begin = high_resolution_clock::now();
  multi_array<double, 4> scalar_flux_fwd{extents[nz][ny][nx][ngf]};
  multi_array<double, 4> scalar_flux_adj{extents[nz][ny][nx][ngf]};
  multi_array<double, 4> scalar_flux_con{extents[nz][ny][nx][ngf]};
  for (int iz = 0; iz < nz; iz++) {                    // Z mesh index
    for (int iy = 0; iy < ny; iy++) {                  // Y mesh index
      for (int ix = 0; ix < nx; ix++) {                // X mesh index
        for (int ia = 0; ia < na; ia++) {              // Angle index
          double qw = quadrature_weights[ia] / FOURPI; // Quadrature weight
          int ja = reverse_angle_map[ia];              // Reverse angle index
          for (int igf = 0; igf < ngf; igf++) {        // Energy index
            double ffwd = angular_flux_fwd[iz][iy][ix][ia][igf] * qw;
            double fadj = angular_flux_adj[iz][iy][ix][ja][igf] * qw;
            scalar_flux_fwd[iz][iy][ix][igf] += ffwd;
            scalar_flux_adj[iz][iy][ix][igf] += fadj;
            scalar_flux_con[iz][iy][ix][igf] += ffwd * fadj;
          }
        }
      }
    }
  }
  end = high_resolution_clock::now();
  elapsed = duration_cast<duration<double>>(end - begin).count();
  cout << " took " << elapsed << " seconds" << endl;
  cout.flush();

  // Calculate cross sections for mixed materials
  cout << "Calculating cross sections for mixed materials...";
  cout.flush();
  begin = high_resolution_clock::now();
  multi_array<double, 2> sigma_t_mixed{extents[n_mix][ngx]};
  multi_array<double, 3> sigma_s_mixed{extents[n_mix][ngx][ngx]};
  for (int i_mix = 0; i_mix < n_mix; i_mix++) { // Mixed material index
    for (int im = 0; im < nm; im++) {           // Material index
      double vol_frac = mix_table[i_mix][im];   // Volume fraction
      if (vol_frac == 0.0)
        continue;
      for (int igx = 0; igx < ngx; igx++) { // Energy index
        sigma_t_mixed[i_mix][igx] += sigma_t[im][igx] * vol_frac;
        for (int jgx = 0; jgx < ngx; jgx++) { // Scattering energy index
          sigma_s_mixed[i_mix][igx][jgx] += sigma_s[im][igx][jgx] * vol_frac;
        }
      }
    }
  }
  end = high_resolution_clock::now();
  elapsed = duration_cast<duration<double>>(end - begin).count();
  cout << " took " << elapsed << " seconds" << endl;
  cout.flush();

  // Calculate perturbed cross sections
  cout << "Calculating perturbed cross sections...";
  cout.flush();
  begin = high_resolution_clock::now();
  multi_array<double, 3> sigma_t_pert{extents[n_mix][nm][ngx]};
  multi_array<double, 4> sigma_s_pert{extents[n_mix][nm][ngx][ngx]};
  for (int i_mix = 0; i_mix < n_mix; i_mix++) { // Mixed material index
    for (int im = 0; im < nm; im++) {           // Material index
      for (int igx = 0; igx < ngx; igx++) {     // Energy index
        sigma_t_pert[i_mix][im][igx] =
            sigma_t[im][igx] - sigma_t_mixed[i_mix][igx];
        for (int jgx = 0; jgx < ngx; jgx++) { // Scattering energy index
          sigma_s_pert[i_mix][im][igx][jgx] =
              sigma_s[im][igx][jgx] - sigma_s_mixed[i_mix][igx][jgx];
        }
      }
    }
  }
  end = high_resolution_clock::now();
  elapsed = duration_cast<duration<double>>(end - begin).count();
  cout << " took " << elapsed << " seconds" << endl;
  cout.flush();

  // Calculate current
  cout << "Calculating current...";
  cout.flush();
  begin = high_resolution_clock::now();
  multi_array<double, 5> current_fwd{extents[nz][ny][nx][ngf][3]};
  multi_array<double, 5> current_adj{extents[nz][ny][nx][ngf][3]};
  multi_array<double, 5> current_con{extents[nz][ny][nx][ngf][3]};
  for (int iz = 0; iz < nz; iz++) {                    // Z mesh index
    for (int iy = 0; iy < ny; iy++) {                  // Y mesh index
      for (int ix = 0; ix < nx; ix++) {                // X mesh index
        for (int ia = 0; ia < na; ia++) {              // Angle index
          double qw = quadrature_weights[ia] / FOURPI; // Quadrature weight
          for (int igf = 0; igf < ngf; igf++) {        // Energy index
            double ffwd = angular_flux_fwd[iz][iy][ix][ia][igf] * qw;
            double fadj = angular_flux_adj[iz][iy][ix][ia][igf] * qw;
            for (int id = 0; id < 3; id++) { // Dimension index
              current_fwd[iz][iy][ix][igf][id] += ffwd * angles[ia][id];
              current_adj[iz][iy][ix][igf][id] += fadj * angles[ia][id];
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
  end = high_resolution_clock::now();
  elapsed = duration_cast<duration<double>>(end - begin).count();
  cout << " took " << elapsed << " seconds" << endl;
  cout.flush();

  // Calculate dR using scalar flux
  cout << "Calculating dR using scalar flux...";
  cout.flush();
  begin = high_resolution_clock::now();
  multi_array<double, 4> dR_scalar{extents[nm][nz][ny][nx]};
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
  end = high_resolution_clock::now();
  elapsed = duration_cast<duration<double>>(end - begin).count();
  cout << " took " << elapsed << " seconds" << endl;
  cout.flush();

  // Calculate dR using angular flux
  cout << "Calculating dR using angular flux... 0";
  cout.flush();
  begin = high_resolution_clock::now();
  multi_array<double, 4> dR_angular{extents[nm][nz][ny][nx]};
  for (int im = 0; im < nm; im++) {                      // Material index
    for (int iz = 0; iz < nz; iz++) {                    // Z mesh index
      for (int iy = 0; iy < ny; iy++) {                  // Y mesh index
        for (int ix = 0; ix < nx; ix++) {                // X mesh index
          int i_mix = material_map[iz][iy][ix];          // Mixed material index
          for (int ia = 0; ia < na; ia++) {              // Angle index
            double qw = quadrature_weights[ia] / FOURPI; // Quadrature weight
            int ja = reverse_angle_map[ia];              // Reverse angle index
            for (int igf = 0; igf < ngf; igf++) { // Energy index (in flux data)
              int igx = igf + g0; // Energy index (in cross section data)
              // Total component of dHphi
              double dHphi_t = sigma_t_pert[i_mix][im][igx] *
                               angular_flux_fwd[iz][iy][ix][ia][igf];
              // Scattering component of dHphi
              double dHphi_s = 0.0;
              for (int jgf = 0; jgf < ngf; jgf++) {
                int jgx = jgf + g0;
                dHphi_s += -sigma_s_pert[i_mix][im][igx][jgx] *
                           angular_flux_fwd[iz][iy][ix][ia][jgf];
              }
              double dHphi = dHphi_t + dHphi_s;
              double dR_comp =
                  -angular_flux_adj[iz][iy][ix][ja][igf] * dHphi * qw;
              dR_angular[im][iz][iy][ix] += dR_comp;
            }
          }
        }
      }
    }
    cout << " " << (100 * (im + 1)) / nm;
    cout.flush();
  }
  end = high_resolution_clock::now();
  elapsed = duration_cast<duration<double>>(end - begin).count();
  cout << " took " << elapsed << " seconds" << endl;
  cout.flush();

  // Write pickles
  write_pickle("reverse_angle_map", reverse_angle_map);
  write_pickle("source", source);
  write_pickle("response", response);
  write_pickle("scalar_flux_fwd", scalar_flux_fwd);
  write_pickle("scalar_flux_adj", scalar_flux_adj);
  write_pickle("scalar_flux_con", scalar_flux_con);
  write_pickle("sigma_t_mixed", sigma_t_mixed);
  write_pickle("sigma_s_mixed", sigma_s_mixed);
  write_pickle("sigma_t_pert", sigma_t_pert);
  write_pickle("sigma_s_pert", sigma_s_pert);
  write_pickle("current_fwd", current_fwd);
  write_pickle("current_adj", current_adj);
  write_pickle("current_con", current_con);
  write_pickle("dR_scalar", dR_scalar);
  write_pickle("dR_angular", dR_angular);

  return 0;
}
