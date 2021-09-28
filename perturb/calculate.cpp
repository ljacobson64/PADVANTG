#include "H5Cpp.h"
#include "boost/multi_array.hpp"
#include <chrono>
#include <iostream>
#include <map>
#include <typeindex>

#define PI 3.14159265358979323846
#define FOURPI (4.0 * PI)

#define TIME_START(s)                                                          \
  {                                                                            \
    printf("%-49s", s);                                                        \
    fflush(stdout);                                                            \
    begin = high_resolution_clock::now();                                      \
  }

#define TIME_END()                                                             \
  {                                                                            \
    end = high_resolution_clock::now();                                        \
    elapsed = duration_cast<duration<double>>(end - begin).count();            \
    printf(" took %7.3f seconds\n", elapsed);                                  \
    fflush(stdout);                                                            \
  }

using boost::array;
using boost::extents;
using boost::multi_array;
using boost::detail::multi_array::multi_array_base;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

struct MixData {
  int row, col;
  double val;
};

// Setup mapping between C++ type and HDF5 data type
H5::DataType get_mixtable_type() {
  H5::CompType mixtable_type(sizeof(MixData));
  mixtable_type.insertMember("row", HOFFSET(MixData, row),
                             H5::PredType::NATIVE_INT32);
  mixtable_type.insertMember("col", HOFFSET(MixData, col),
                             H5::PredType::NATIVE_INT32);
  mixtable_type.insertMember("val", HOFFSET(MixData, val),
                             H5::PredType::NATIVE_DOUBLE);
  return mixtable_type;
}
template <typename T> struct get_hdf5_data_type {
  static const H5::DataType type;
};
template <>
const H5::DataType get_hdf5_data_type<int>::type = H5::PredType::NATIVE_INT32;
template <>
const H5::DataType get_hdf5_data_type<double>::type =
    H5::PredType::NATIVE_DOUBLE;
template <>
const H5::DataType get_hdf5_data_type<MixData>::type = get_mixtable_type();

template <typename T, std::size_t NDIMS>
void read_hdf5_array(H5::H5File hf, const std::string &dname,
                     multi_array<T, NDIMS> &A) {
  /* Read an HDF5 array into a boost::multi_array.
   */
  H5::DataType dtype = get_hdf5_data_type<T>::type;
  H5::DataSet dset = hf.openDataSet(dname.c_str());
  H5::DataSpace dspace = dset.getSpace();
  hsize_t dims[NDIMS];
  dspace.getSimpleExtentDims(dims);
  array<multi_array_base::index, NDIMS> A_extents;
  for (int i = 0; i < NDIMS; i++) {
    A_extents[i] = dims[i];
  }
  A.resize(A_extents);
  dset.read(A.origin(), dtype, dspace, dspace);
}

template <typename T, std::size_t NDIMS>
void write_hdf5_array(H5::H5File hf, const std::string &dname,
                      const multi_array<T, NDIMS> &A) {
  /* Write a boost::multi_array to HDF5.
   */
  H5::DataType dtype = get_hdf5_data_type<T>::type;
  std::vector<hsize_t> dims(A.shape(), A.shape() + NDIMS);
  H5::DataSpace dspace(NDIMS, dims.data());
  H5::DataSet dset = hf.createDataSet(dname, dtype, dspace);
  dset.write(A.data(), dtype);
}

int main() {
  // Variables for timing
  high_resolution_clock::time_point begin, end;
  double elapsed;

  // Create HDF5 file objects
  std::string fname_fwd_inp = "fwd_solution/denovo-forward.inp.h5";
  std::string fname_adj_inp = "adj_solution/denovo-adjoint.inp.h5";
  std::string fname_fwd_out = "fwd_solution/denovo-forward.out.h5";
  std::string fname_adj_out = "adj_solution/denovo-adjoint.out.h5";
  std::string fname_xs = "custom_output/xs.h5";
  H5::H5File hf_fwd_inp(fname_fwd_inp.c_str(), H5F_ACC_RDONLY);
  H5::H5File hf_adj_inp(fname_adj_inp.c_str(), H5F_ACC_RDONLY);
  H5::H5File hf_fwd_out(fname_fwd_out.c_str(), H5F_ACC_RDONLY);
  H5::H5File hf_adj_out(fname_adj_out.c_str(), H5F_ACC_RDONLY);
  H5::H5File hf_xs(fname_xs.c_str(), H5F_ACC_RDONLY);

  // Create arrays for data to be read from HDF5
  multi_array<int, 1> mesh_g;
  multi_array<int, 3> matids, source_ids, response_ids;
  multi_array<double, 1> group_bounds_n, group_bounds_p, mesh_x, mesh_y, mesh_z,
      quadrature_weights;
  multi_array<double, 2> source_spectra, response_spectra, quadrature_angles,
      sigma_t;
  multi_array<double, 3> source_strength, response_strength, sigma_s;
  multi_array<double, 5> angular_flux_fwd, angular_flux_adj;
  multi_array<MixData, 1> mixtable;

  // Read data from HDF5
  TIME_START("Reading data from HDF5...");
  read_hdf5_array(hf_fwd_inp, "group_bounds_n", group_bounds_n);
  read_hdf5_array(hf_fwd_inp, "group_bounds_p", group_bounds_p);
  read_hdf5_array(hf_fwd_inp, "mixtable", mixtable);
  read_hdf5_array(hf_fwd_inp, "matids", matids);
  read_hdf5_array(hf_fwd_inp, "volsrc/spectra", source_spectra);
  read_hdf5_array(hf_fwd_inp, "volsrc/strength", source_strength);
  read_hdf5_array(hf_fwd_inp, "volsrc/ids", source_ids);
  read_hdf5_array(hf_adj_inp, "volsrc/spectra", response_spectra);
  read_hdf5_array(hf_adj_inp, "volsrc/strength", response_strength);
  read_hdf5_array(hf_adj_inp, "volsrc/ids", response_ids);
  read_hdf5_array(hf_fwd_out, "denovo/mesh_x", mesh_x);
  read_hdf5_array(hf_fwd_out, "denovo/mesh_y", mesh_y);
  read_hdf5_array(hf_fwd_out, "denovo/mesh_z", mesh_z);
  read_hdf5_array(hf_fwd_out, "denovo/mesh_g", mesh_g);
  read_hdf5_array(hf_fwd_out, "denovo/quadrature_angles", quadrature_angles);
  read_hdf5_array(hf_fwd_out, "denovo/quadrature_weights", quadrature_weights);
  read_hdf5_array(hf_fwd_out, "denovo/angular_flux", angular_flux_fwd);
  read_hdf5_array(hf_adj_out, "denovo/angular_flux", angular_flux_adj);
  read_hdf5_array(hf_xs, "sigma_t", sigma_t);
  read_hdf5_array(hf_xs, "sigma_s", sigma_s);
  TIME_END();

  // +--------------------------------------------------------------------+
  // |                    Array extents for reference                     |
  // +--------------------+-----------------------+-----------------------+
  // | group_bounds_n     | (28,)                 |                       |
  // | group_bounds_p     | (20,)                 |                       |
  // | mixtable           | (17,)                 |                       |
  // | matids             | (44, 44, 43)          | [nz][ny][nx]          |
  // | source_spectra     | (1, 46)               | [0][ngx]              |
  // | source_strength    | (44, 44, 43)          | [nz][ny][nx]          |
  // | source_ids         | (44, 44, 43)          | [nz][ny][nx]          |
  // | response_spectra   | (1, 46)               | [0][ngx]              |
  // | response_strength  | (44, 44, 43)          | [nz][ny][nx]          |
  // | response_ids       | (44, 44, 43)          | [nz][ny][nx]          |
  // | mesh_x             | (44,)                 | [nx + 1]              |
  // | mesh_y             | (45,)                 | [ny + 1]              |
  // | mesh_z             | (45,)                 | [nz + 1]              |
  // | mesh_g             | (26,)                 | [ngf]                 |
  // | quadrature_angles  | (128, 3)              | [na][3]               |
  // | quadrature_weights | (128,)                | [na]                  |
  // | angular_flux_fwd   | (26, 44, 44, 43, 128) | [ngf][nz][ny][nx][na] |
  // | angular_flux_adj   | (26, 44, 44, 43, 128) | [ngf][nz][ny][nx][na] |
  // | sigma_t            | (17, 46)              | [nm][ngx]             |
  // | sigma_s            | (17, 46, 46)          | [nm][ngx][ngx]        |
  // +--------------------+-----------------------+-----------------------+

  // Relevant dimensions
  int nm = sigma_t.shape()[0];           // Number of pure materials
  int ngx = sigma_t.shape()[1];          // Number of energy groups in XS
  int ngf = angular_flux_fwd.shape()[0]; // Number of energy groups in flux
  int nz = angular_flux_fwd.shape()[1];  // Number of Z intervals
  int ny = angular_flux_fwd.shape()[2];  // Number of Y intervals
  int nx = angular_flux_fwd.shape()[3];  // Number of X intervals
  int na = angular_flux_fwd.shape()[4];  // Number of angles

  // First energy group in flux
  int g0 = mesh_g[0];

  // Number of mixed materials
  int nmix = mixtable[mixtable.shape()[0] - 1].row + 1;

  // Calculate reverse angle map
  TIME_START("Calculating reverse angle map...");
  multi_array<int, 1> reverse_angle_map{extents[na]};
  for (int ia = 0; ia < na; ia++) {
    double ix = quadrature_angles[ia][0];
    double iy = quadrature_angles[ia][1];
    double iz = quadrature_angles[ia][2];
    for (int ja = 0; ja < na; ja++) {
      double jx = quadrature_angles[ja][0];
      double jy = quadrature_angles[ja][1];
      double jz = quadrature_angles[ja][2];
      if (ix == -jx and iy == -jy and iz == -jz) {
        reverse_angle_map[ia] = ja;
        break;
      }
    }
  }
  TIME_END();

  // Calculate source and response
  TIME_START("Calculating source and response...");
  multi_array<double, 4> source{extents[ngf][nz][ny][nx]};
  multi_array<double, 4> response{extents[ngf][nz][ny][nx]};
  for (int igf = 0; igf < ngf; igf++) { // Energy index
    double spec_src = source_spectra[0][g0 + igf];
    double spec_res = response_spectra[0][g0 + igf];
    for (int iz = 0; iz < nz; iz++) {     // Z mesh index
      for (int iy = 0; iy < ny; iy++) {   // Y mesh index
        for (int ix = 0; ix < nx; ix++) { // X mesh index
          source[igf][iz][iy][ix] = source_strength[iz][iy][ix] * spec_src;
          response[igf][iz][iy][ix] = response_strength[iz][iy][ix] * spec_res;
        }
      }
    }
  }
  TIME_END();

  // Calculate cross sections for mixed materials
  TIME_START("Calculating cross sections for mixed materials...");
  multi_array<double, 2> sigma_t_mixed{extents[nmix][ngx]};
  multi_array<double, 3> sigma_s_mixed{extents[nmix][ngx][ngx]};
  for (int i = 0; i < mixtable.shape()[0]; i++) { // Mix table entry index
    int imix = mixtable[i].row;                   // Mixed material index
    int im = mixtable[i].col;                     // Pure material index
    double vol_frac = mixtable[i].val;            // Volume fraction
    for (int igx = 0; igx < ngx; igx++) {         // Energy index
      sigma_t_mixed[imix][igx] += sigma_t[im][igx] * vol_frac;
      for (int jgx = 0; jgx < ngx; jgx++) { // Scattering energy index
        sigma_s_mixed[imix][igx][jgx] += sigma_s[im][igx][jgx] * vol_frac;
      }
    }
  }
  TIME_END();

  // Calculate perturbations in cross sections
  TIME_START("Calculating perturbations in cross sections...");
  multi_array<double, 3> sigma_t_pert{extents[nmix][nm][ngx]};
  multi_array<double, 4> sigma_s_pert{extents[nmix][nm][ngx][ngx]};
  for (int imix = 0; imix < nmix; imix++) { // Mixed material index
    for (int im = 0; im < nm; im++) {       // Material index
      for (int igx = 0; igx < ngx; igx++) { // Energy index
        sigma_t_pert[imix][im][igx] =
            sigma_t[im][igx] - sigma_t_mixed[imix][igx];
        for (int jgx = 0; jgx < ngx; jgx++) { // Scattering energy index
          sigma_s_pert[imix][im][igx][jgx] =
              sigma_s[im][igx][jgx] - sigma_s_mixed[imix][igx][jgx];
        }
      }
    }
  }
  TIME_END();

  // Calculate scalar flux
  TIME_START("Calculating scalar flux...");
  multi_array<double, 4> scalar_flux_fwd{extents[ngf][nz][ny][nx]};
  multi_array<double, 4> scalar_flux_adj{extents[ngf][nz][ny][nx]};
  multi_array<double, 4> scalar_flux_con{extents[ngf][nz][ny][nx]};
  for (int igf = 0; igf < ngf; igf++) {                  // Energy index
    for (int iz = 0; iz < nz; iz++) {                    // Z mesh index
      for (int iy = 0; iy < ny; iy++) {                  // Y mesh index
        for (int ix = 0; ix < nx; ix++) {                // X mesh index
          for (int ia = 0; ia < na; ia++) {              // Angle index
            double qw = quadrature_weights[ia] / FOURPI; // Quadrature weight
            int ja = reverse_angle_map[ia];              // Reverse angle index
            double ffwd = angular_flux_fwd[igf][iz][iy][ix][ia] * qw;
            double fadj = angular_flux_adj[igf][iz][iy][ix][ja] * qw;
            scalar_flux_fwd[igf][iz][iy][ix] += ffwd;
            scalar_flux_adj[igf][iz][iy][ix] += fadj;
            scalar_flux_con[igf][iz][iy][ix] += ffwd * fadj;
          }
        }
      }
    }
  }
  TIME_END();

  // Calculate current
  TIME_START("Calculating current...");
  multi_array<double, 5> current_fwd{extents[ngf][nz][ny][nx][3]};
  multi_array<double, 5> current_adj{extents[ngf][nz][ny][nx][3]};
  multi_array<double, 5> current_con{extents[ngf][nz][ny][nx][3]};
  for (int igf = 0; igf < ngf; igf++) {                  // Energy index
    for (int iz = 0; iz < nz; iz++) {                    // Z mesh index
      for (int iy = 0; iy < ny; iy++) {                  // Y mesh index
        for (int ix = 0; ix < nx; ix++) {                // X mesh index
          for (int ia = 0; ia < na; ia++) {              // Angle index
            double qw = quadrature_weights[ia] / FOURPI; // Quadrature weight
            double affwd = angular_flux_fwd[igf][iz][iy][ix][ia] * qw;
            double afadj = angular_flux_adj[igf][iz][iy][ix][ia] * qw;
            for (int id = 0; id < 3; id++) { // Dimension index
              double qa = quadrature_angles[ia][id];
              current_fwd[igf][iz][iy][ix][id] += affwd * qa;
              current_adj[igf][iz][iy][ix][id] += afadj * qa;
            }
          }
          for (int id = 0; id < 3; id++) {
            current_con[igf][iz][iy][ix][id] =
                current_fwd[igf][iz][iy][ix][id] *
                current_adj[igf][iz][iy][ix][id];
          }
        }
      }
    }
  }
  TIME_END();

  // Calculate dR using angular flux
  TIME_START("Calculating dR...");
  multi_array<double, 4> dR{extents[nm][nz][ny][nx]};
  for (int im = 0; im < nm; im++) {       // Pure material index
    for (int iz = 0; iz < nz; iz++) {     // Z mesh index
      for (int iy = 0; iy < ny; iy++) {   // Y mesh index
        for (int ix = 0; ix < nx; ix++) { // X mesh index
          int imix = matids[iz][iy][ix];  // Mixed material index
          // Total component of dR
          double dR_t = 0.0;
          for (int igf = 0; igf < ngf; igf++) { // Energy index (in flux data)
            int igx = igf + g0; // Energy index (in cross section data)
            // Scalar contributon flux
            double sfc = scalar_flux_con[igf][iz][iy][ix];
            // Perturbation in total cross section
            double dst = sigma_t_pert[imix][im][igx];
            dR_t += sfc * dst;
          }
          // Scattering component of dR
          double dR_s = 0.0;
          for (int igf = 0; igf < ngf; igf++) { // Energy index (in flux data)
            int igx = igf + g0; // Energy index (in cross section data)
            // Scalar adjoint flux (at E)
            double sfa = scalar_flux_adj[igf][iz][iy][ix];
            for (int jgf = 0; jgf < ngf; jgf++) { // E' index (in flux data)
              int jgx = jgf + g0; // E' index (in cross section data)
              // Perturbation in scattering cross section (E' -> E)
              double dss = sigma_s_pert[imix][im][igx][jgx];
              // Scalar forward flux (at E')
              double sff = scalar_flux_fwd[jgf][iz][iy][ix];
              dR_s += dss * sff * sfa;
            }
          }
          // Calculate dR
          dR[im][iz][iy][ix] = dR_s - dR_t;
        }
      }
    }
  }
  TIME_END();

  // Write data to HDF5
  TIME_START("Writing data to HDF5...");
  std::string fnameo = "custom_output/data.h5";
  H5::H5File hfo(fnameo, H5F_ACC_TRUNC);
  write_hdf5_array(hfo, "reverse_angle_map", reverse_angle_map);
  write_hdf5_array(hfo, "source", source);
  write_hdf5_array(hfo, "response", response);
  write_hdf5_array(hfo, "sigma_t_mixed", sigma_t_mixed);
  write_hdf5_array(hfo, "sigma_s_mixed", sigma_s_mixed);
  write_hdf5_array(hfo, "sigma_t_pert", sigma_t_pert);
  write_hdf5_array(hfo, "sigma_s_pert", sigma_s_pert);
  write_hdf5_array(hfo, "scalar_flux_fwd", scalar_flux_fwd);
  write_hdf5_array(hfo, "scalar_flux_adj", scalar_flux_adj);
  write_hdf5_array(hfo, "scalar_flux_con", scalar_flux_con);
  write_hdf5_array(hfo, "current_fwd", current_fwd);
  write_hdf5_array(hfo, "current_adj", current_adj);
  write_hdf5_array(hfo, "current_con", current_con);
  write_hdf5_array(hfo, "dR", dR);
  TIME_END();

  return 0;
}
