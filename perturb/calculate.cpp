#include "H5Cpp.h"
#include "boost/multi_array.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <typeindex>
#include <vector>

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
using std::vector;
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

  // Get tally IDs
  std::string fname = "tally_list.txt";
  std::ifstream infile(fname);
  std::string buffer;
  vector<int> tallies;
  while (std::getline(infile, buffer)) {
    int tally = std::stof(buffer);
    tallies.push_back(tally);
  }
  int nas = tallies.size();

  // Create HDF5 file objects
  H5::H5File hf_xs("hdf5/advantg_xs.h5", H5F_ACC_RDONLY);
  H5::H5File hf_fi("hdf5/advantg_fwd_inp.h5", H5F_ACC_RDONLY);
  H5::H5File hf_fo("hdf5/advantg_fwd_out.h5", H5F_ACC_RDONLY);
  vector<H5::H5File> hf_ai(nas), hf_ao(nas);
  for (int ias = 0; ias < nas; ias++) {
    std::string fname_ai =
        "hdf5/advantg_adj_" + std::to_string(tallies[ias]) + "_inp.h5";
    std::string fname_ao =
        "hdf5/advantg_adj_" + std::to_string(tallies[ias]) + "_out.h5";
    hf_ai[ias].openFile(fname_ai.c_str(), H5F_ACC_RDONLY);
    hf_ao[ias].openFile(fname_ao.c_str(), H5F_ACC_RDONLY);
  }

  // Create arrays for miscellaneous data
  multi_array<int, 3> matids;
  multi_array<double, 1> quadrature_weights;
  multi_array<double, 2> sigma_t, quadrature_angles;
  multi_array<double, 3> sigma_s;
  multi_array<MixData, 1> mixtable;

  // Create arrays for forward data
  multi_array<int, 1> mesh_g_fwd;
  multi_array<double, 2> source_spectra_fwd;
  multi_array<double, 3> source_strength_fwd;
  multi_array<double, 5> angular_flux_fwd;

  // Create vectors of arrays for adjoint data
  vector<multi_array<int, 1>> mesh_g_adj(nas);
  vector<multi_array<double, 2>> source_spectra_adj(nas);
  vector<multi_array<double, 3>> source_strength_adj(nas);
  vector<multi_array<double, 5>> angular_flux_adj(nas);

  // Read data from HDF5
  TIME_START("Reading data from HDF5...");
  read_hdf5_array(hf_xs, "sigma_t", sigma_t);
  read_hdf5_array(hf_xs, "sigma_s", sigma_s);
  read_hdf5_array(hf_fi, "mixtable", mixtable);
  read_hdf5_array(hf_fi, "matids", matids);
  read_hdf5_array(hf_fo, "denovo/quadrature_angles", quadrature_angles);
  read_hdf5_array(hf_fo, "denovo/quadrature_weights", quadrature_weights);
  read_hdf5_array(hf_fi, "volsrc/spectra", source_spectra_fwd);
  read_hdf5_array(hf_fi, "volsrc/strength", source_strength_fwd);
  read_hdf5_array(hf_fo, "denovo/mesh_g", mesh_g_fwd);
  read_hdf5_array(hf_fo, "denovo/angular_flux", angular_flux_fwd);
  for (int ias = 0; ias < nas; ias++) {
    read_hdf5_array(hf_ai[ias], "volsrc/spectra", source_spectra_adj[ias]);
    read_hdf5_array(hf_ai[ias], "volsrc/strength", source_strength_adj[ias]);
    read_hdf5_array(hf_ao[ias], "denovo/mesh_g", mesh_g_adj[ias]);
    read_hdf5_array(hf_ao[ias], "denovo/angular_flux", angular_flux_adj[ias]);
  }
  TIME_END();

  // +----------------------------------------------------------------------+
  // |                      Array extents for reference                     |
  // +------------------------+-----------------------+---------------------+
  // | sigma_t                | (39, 46)              | nm,ngx              |
  // | sigma_s                | (39, 46, 46)          | nm,ngx,ngx          |
  // | mixtable               | (1281,)               |                     |
  // | matids                 | (45, 45, 45)          | nz,ny,nx            |
  // | quadrature_angles      | (128, 3)              | na,3                |
  // | quadrature_weights     | (128,)                | na                  |
  // | source_spectra_fwd     | (1, 46)               | 0,ngx               |
  // | source_strength_fwd    | (45, 45, 45)          | nz,ny,nx            |
  // | mesh_g_fwd             | (46,)                 | ngff                |
  // | angular_flux_fwd       | (46, 45, 45, 45, 128) | ngff,nz,ny,nx,na    |
  // | source_spectra_adj [0] | (1, 46)               | 0,ngx               |
  // | source_strength_adj[0] | (45, 45, 45)          | nz,ny,nx            |
  // | mesh_g_adj         [0] | (20,)                 | ngfa[0]             |
  // | angular_flux_adj   [0] | (20, 45, 45, 45, 128) | ngfa[0],nz,ny,nx,na |
  // +------------------------+-----------------------+---------------------+

  // Relevant dimensions
  int nm = sigma_t.shape()[0];          // Number of pure materials
  int ngx = sigma_t.shape()[1];         // Number of energy groups in XS
  int nz = angular_flux_fwd.shape()[1]; // Number of Z intervals
  int ny = angular_flux_fwd.shape()[2]; // Number of Y intervals
  int nx = angular_flux_fwd.shape()[3]; // Number of X intervals
  int na = angular_flux_fwd.shape()[4]; // Number of angles

  // Number of energy groups
  int ngff = angular_flux_fwd.shape()[0]; // Forward flux
  vector<int> ngfa(nas);
  for (int ias = 0; ias < nas; ias++) {
    ngfa[ias] = angular_flux_adj[ias].shape()[0]; // Adjoint flux
  }

  // Number of mixed materials
  int nmix = mixtable[mixtable.shape()[0] - 1].row + 1;

  // First energy group in flux
  int g0f = mesh_g_fwd[0]; // Forward flux
  vector<int> g0a(nas);
  for (int ias = 0; ias < nas; ias++) {
    g0a[ias] = mesh_g_adj[ias][0];
  }

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

  // Calculate 4D sources
  TIME_START("Calculating 4D sources...");
  multi_array<double, 4> source_fwd{extents[ngx][nz][ny][nx]};
  multi_array<double, 5> source_adj{extents[nas][ngx][nz][ny][nx]};
  for (int igx = 0; igx < ngx; igx++) { // Energy index
    double spec_fwd = source_spectra_fwd[0][igx];
    for (int iz = 0; iz < nz; iz++) {     // Z mesh index
      for (int iy = 0; iy < ny; iy++) {   // Y mesh index
        for (int ix = 0; ix < nx; ix++) { // X mesh index
          source_fwd[igx][iz][iy][ix] =
              source_strength_fwd[iz][iy][ix] * spec_fwd;
        }
      }
    }
  }
  for (int ias = 0; ias < nas; ias++) {   // Adjoint source index
    for (int igx = 0; igx < ngx; igx++) { // Energy index
      double spec_adj = source_spectra_adj[ias][0][igx];
      for (int iz = 0; iz < nz; iz++) {     // Z mesh index
        for (int iy = 0; iy < ny; iy++) {   // Y mesh index
          for (int ix = 0; ix < nx; ix++) { // X mesh index
            source_adj[ias][igx][iz][iy][ix] =
                source_strength_adj[ias][iz][iy][ix] * spec_adj;
          }
        }
      }
    }
  }
  TIME_END();

  // Calculate scalar flux
  TIME_START("Calculating scalar flux...");
  multi_array<double, 4> scalar_flux_fwd{extents[ngx][nz][ny][nx]};
  multi_array<double, 5> scalar_flux_adj{extents[nas][ngx][nz][ny][nx]};
  multi_array<double, 5> scalar_flux_con{extents[nas][ngx][nz][ny][nx]};
  for (int igx = 0; igx < ngx; igx++) {     // Energy index
    for (int iz = 0; iz < nz; iz++) {       // Z mesh index
      for (int iy = 0; iy < ny; iy++) {     // Y mesh index
        for (int ix = 0; ix < nx; ix++) {   // X mesh index
          for (int ia = 0; ia < na; ia++) { // Angle index
            double qw = quadrature_weights[ia] / FOURPI;
            double ffwd = 0.0;
            if (igx - g0f >= 0 && igx - g0f < ngff) {
              ffwd = angular_flux_fwd[igx - g0f][iz][iy][ix][ia] * qw;
            }
            scalar_flux_fwd[igx][iz][iy][ix] += ffwd;
          }
        }
      }
    }
  }
  for (int ias = 0; ias < nas; ias++) {       // Adjoint source index
    for (int igx = 0; igx < ngx; igx++) {     // Energy index
      for (int iz = 0; iz < nz; iz++) {       // Z mesh index
        for (int iy = 0; iy < ny; iy++) {     // Y mesh index
          for (int ix = 0; ix < nx; ix++) {   // X mesh index
            for (int ia = 0; ia < na; ia++) { // Angle index
              double qw = quadrature_weights[ia] / FOURPI;
              int ja = reverse_angle_map[ia];
              double ffwd = 0.0;
              double fadj = 0.0;
              if (igx - g0f >= 0 && igx - g0f < ngff) {
                ffwd = angular_flux_fwd[igx - g0f][iz][iy][ix][ia] * qw;
              }
              if (igx - g0a[ias] >= 0 && igx - g0a[ias] < ngfa[ias]) {
                fadj =
                    angular_flux_adj[ias][igx - g0a[ias]][iz][iy][ix][ja] * qw;
              }
              scalar_flux_adj[ias][igx][iz][iy][ix] += fadj;
              scalar_flux_con[ias][igx][iz][iy][ix] += ffwd * fadj;
            }
          }
        }
      }
    }
  }
  TIME_END();

  // Calculate current
  TIME_START("Calculating current...");
  multi_array<double, 5> current_fwd{extents[ngx][nz][ny][nx][3]};
  multi_array<double, 6> current_adj{extents[nas][ngx][nz][ny][nx][3]};
  multi_array<double, 6> current_con{extents[nas][ngx][nz][ny][nx][3]};
  for (int igx = 0; igx < ngx; igx++) {     // Energy index
    for (int iz = 0; iz < nz; iz++) {       // Z mesh index
      for (int iy = 0; iy < ny; iy++) {     // Y mesh index
        for (int ix = 0; ix < nx; ix++) {   // X mesh index
          for (int ia = 0; ia < na; ia++) { // Angle index
            double qw = quadrature_weights[ia] / FOURPI;
            double affwd = 0.0;
            if (igx - g0f >= 0 && igx - g0f < ngff) {
              affwd = angular_flux_fwd[igx - g0f][iz][iy][ix][ia] * qw;
            }
            for (int id = 0; id < 3; id++) { // Dimension index
              double qa = quadrature_angles[ia][id];
              current_fwd[igx][iz][iy][ix][id] += affwd * qa;
            }
          }
        }
      }
    }
  }
  for (int ias = 0; ias < nas; ias++) {       // Adjoint source index
    for (int igx = 0; igx < ngx; igx++) {     // Energy index
      for (int iz = 0; iz < nz; iz++) {       // Z mesh index
        for (int iy = 0; iy < ny; iy++) {     // Y mesh index
          for (int ix = 0; ix < nx; ix++) {   // X mesh index
            for (int ia = 0; ia < na; ia++) { // Angle index
              double qw = quadrature_weights[ia] / FOURPI;
              double afadj = 0.0;
              if (igx - g0a[ias] >= 0 && igx - g0a[ias] < ngfa[ias]) {
                afadj =
                    angular_flux_adj[ias][igx - g0a[ias]][iz][iy][ix][ia] * qw;
              }
              for (int id = 0; id < 3; id++) { // Dimension index
                double qa = quadrature_angles[ia][id];
                current_adj[ias][igx][iz][iy][ix][id] += afadj * qa;
              }
            }
            for (int id = 0; id < 3; id++) { // Dimension index
              current_con[ias][igx][iz][iy][ix][id] =
                  current_fwd[igx][iz][iy][ix][id] *
                  current_adj[ias][igx][iz][iy][ix][id];
            }
          }
        }
      }
    }
  }
  TIME_END();

  // Calculate response
  multi_array<double, 1> response{extents[nas]};
  TIME_START("Calculating response...");
  for (int igx = 0; igx < ngx; igx++) {   // Energy index
    for (int iz = 0; iz < nz; iz++) {     // Z mesh index
      for (int iy = 0; iy < ny; iy++) {   // Y mesh index
        for (int ix = 0; ix < nx; ix++) { // X mesh index
          double sff = scalar_flux_fwd[igx][iz][iy][ix];
          for (int ias = 0; ias < nas; ias++) {
            double sa = source_adj[ias][igx][iz][iy][ix];
            response[ias] += sff * sa;
          }
        }
      }
    }
  }
  TIME_END();

  // Calculate dR using angular flux
  TIME_START("Calculating dR...");
  multi_array<double, 5> dR{extents[nas][nm][nz][ny][nx]};
  for (int ias = 0; ias < nas; ias++) {     // Adjoint source index
    for (int im = 0; im < nm; im++) {       // Pure material index
      for (int iz = 0; iz < nz; iz++) {     // Z mesh index
        for (int iy = 0; iy < ny; iy++) {   // Y mesh index
          for (int ix = 0; ix < nx; ix++) { // X mesh index
            int imix = matids[iz][iy][ix];  // Mixed material index
            // Total component of dR
            double dR_t = 0.0;
            for (int igx = 0; igx < ngx; igx++) { // Energy index
              // Scalar contributon flux
              double sfc = scalar_flux_con[ias][igx][iz][iy][ix];
              // Perturbation in total cross section
              double dst = sigma_t_pert[imix][im][igx];
              dR_t += sfc * dst;
            }
            // Scattering component of dR
            double dR_s = 0.0;
            for (int igx = 0; igx < ngx; igx++) { // Energy index
              // Scalar adjoint flux (at E)
              double sfa = scalar_flux_adj[ias][igx][iz][iy][ix];
              for (int jgx = 0; jgx < ngx; jgx++) { // E' index
                // Perturbation in scattering cross section (E' -> E)
                double dss = sigma_s_pert[imix][im][igx][jgx];
                // Scalar forward flux (at E')
                double sff = scalar_flux_fwd[jgx][iz][iy][ix];
                dR_s += dss * sff * sfa;
              }
            }
            // Calculate dR
            dR[ias][im][iz][iy][ix] = dR_s - dR_t;
          }
        }
      }
    }
  }
  TIME_END();

  // Write data to HDF5
  TIME_START("Writing data to HDF5...");
  H5::H5File hf_o("hdf5/data.h5", H5F_ACC_TRUNC);
  write_hdf5_array(hf_o, "sigma_t_mixed", sigma_t_mixed);
  write_hdf5_array(hf_o, "sigma_s_mixed", sigma_s_mixed);
  write_hdf5_array(hf_o, "sigma_t_pert", sigma_t_pert);
  write_hdf5_array(hf_o, "sigma_s_pert", sigma_s_pert);
  write_hdf5_array(hf_o, "source_fwd", source_fwd);
  write_hdf5_array(hf_o, "source_adj", source_adj);
  write_hdf5_array(hf_o, "scalar_flux_fwd", scalar_flux_fwd);
  write_hdf5_array(hf_o, "scalar_flux_adj", scalar_flux_adj);
  write_hdf5_array(hf_o, "scalar_flux_con", scalar_flux_con);
  write_hdf5_array(hf_o, "current_fwd", current_fwd);
  write_hdf5_array(hf_o, "current_adj", current_adj);
  write_hdf5_array(hf_o, "current_con", current_con);
  write_hdf5_array(hf_o, "response", response);
  write_hdf5_array(hf_o, "dR", dR);
  TIME_END();

  return 0;
}
