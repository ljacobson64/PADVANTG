#include "H5Cpp.h"
#include "boost/multi_array.hpp"
#include <chrono>
#include <iostream>

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

using boost::extents;
using boost::multi_array;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

typedef multi_array<int, 1> ArrayInt1D;
typedef multi_array<int, 3> ArrayInt3D;
typedef multi_array<double, 1> ArrayDouble1D;
typedef multi_array<double, 2> ArrayDouble2D;
typedef multi_array<double, 3> ArrayDouble3D;
typedef multi_array<double, 4> ArrayDouble4D;
typedef multi_array<double, 5> ArrayDouble5D;

// template <typename T>
// void load_int_array(H5::H5File hf, std::string dname, T array) {
//  H5::DataSet dset = hf.openDataSet(dname.c_str());
//  H5::DataSpace dspace = dset.getSpace();
//  hsize_t rank = dspace.getSimpleExtentNdims();
//  hsize_t dims[rank];
//  dspace.getSimpleExtentDims(dims);
//  if (rank == 1)
//    array.resize(extents[dims[0]]);
//  else if (rank == 2)
//    array.resize(extents[dims[0]][dims[1]]);
//  else if (rank == 3)
//    array.resize(extents[dims[0]][dims[1]][dims[2]]);
//  else if (rank == 4)
//    array.resize(extents[dims[0]][dims[1]][dims[2]]);
//  else if (rank == 5)
//    array.resize(extents[dims[0]][dims[1]][dims[2]][dims[3]][dims[4]]);
//  dset.read(array.origin(), H5::PredType::NATIVE_INT, dspace, dspace);
//}

ArrayInt1D load_int_1D(H5::H5File hf, std::string dname) {
  H5::DataSet dset = hf.openDataSet(dname.c_str());
  H5::DataSpace dspace = dset.getSpace();
  hsize_t rank = dspace.getSimpleExtentNdims();
  hsize_t dims[rank];
  dspace.getSimpleExtentDims(dims);
  ArrayInt1D array;
  array.resize(extents[dims[0]]);
  dset.read(array.origin(), H5::PredType::NATIVE_INT, dspace, dspace);
  return array;
}

ArrayInt3D load_int_3D(H5::H5File hf, std::string dname) {
  H5::DataSet dset = hf.openDataSet(dname.c_str());
  H5::DataSpace dspace = dset.getSpace();
  hsize_t rank = dspace.getSimpleExtentNdims();
  hsize_t dims[rank];
  dspace.getSimpleExtentDims(dims);
  ArrayInt3D array;
  array.resize(extents[dims[0]][dims[1]][dims[2]]);
  dset.read(array.origin(), H5::PredType::NATIVE_INT, dspace, dspace);
  return array;
}

ArrayDouble1D load_double_1D(H5::H5File hf, std::string dname) {
  H5::DataSet dset = hf.openDataSet(dname.c_str());
  H5::DataSpace dspace = dset.getSpace();
  hsize_t rank = dspace.getSimpleExtentNdims();
  hsize_t dims[rank];
  dspace.getSimpleExtentDims(dims);
  ArrayDouble1D array;
  array.resize(extents[dims[0]]);
  dset.read(array.origin(), H5::PredType::NATIVE_DOUBLE, dspace, dspace);
  return array;
}

ArrayDouble2D load_double_2D(H5::H5File hf, std::string dname) {
  H5::DataSet dset = hf.openDataSet(dname.c_str());
  H5::DataSpace dspace = dset.getSpace();
  hsize_t rank = dspace.getSimpleExtentNdims();
  hsize_t dims[rank];
  dspace.getSimpleExtentDims(dims);
  ArrayDouble2D array;
  array.resize(extents[dims[0]][dims[1]]);
  dset.read(array.origin(), H5::PredType::NATIVE_DOUBLE, dspace, dspace);
  return array;
}

ArrayDouble3D load_double_3D(H5::H5File hf, std::string dname) {
  H5::DataSet dset = hf.openDataSet(dname.c_str());
  H5::DataSpace dspace = dset.getSpace();
  hsize_t rank = dspace.getSimpleExtentNdims();
  hsize_t dims[rank];
  dspace.getSimpleExtentDims(dims);
  ArrayDouble3D array;
  array.resize(extents[dims[0]][dims[1]][dims[2]]);
  dset.read(array.origin(), H5::PredType::NATIVE_DOUBLE, dspace, dspace);
  return array;
}

ArrayDouble5D load_double_5D(H5::H5File hf, std::string dname) {
  H5::DataSet dset = hf.openDataSet(dname.c_str());
  H5::DataSpace dspace = dset.getSpace();
  hsize_t rank = dspace.getSimpleExtentNdims();
  hsize_t dims[rank];
  dspace.getSimpleExtentDims(dims);
  ArrayDouble5D array;
  array.resize(extents[dims[0]][dims[1]][dims[2]][dims[3]][dims[4]]);
  dset.read(array.origin(), H5::PredType::NATIVE_DOUBLE, dspace, dspace);
  return array;
}

template <typename T, std::size_t DIMENSIONS, typename hdf5_data_type>
void write_hdf5_array(H5::H5File hf, const std::string &dname,
                      const boost::multi_array<T, DIMENSIONS> &array,
                      hdf5_data_type &dtype) {
  /* Write a boost::multi_array to HDF5.
     Source: https://stackoverflow.com/a/15221213
   */
  dtype.setOrder(H5T_ORDER_LE);
  std::vector<hsize_t> dimensions(array.shape(), array.shape() + DIMENSIONS);
  H5::DataSpace dataspace(DIMENSIONS, dimensions.data());
  H5::DataSet dataset = hf.createDataSet(dname, dtype, dataspace);
  dataset.write(array.data(), dtype);
}

int main() {
  // Variables for timing
  high_resolution_clock::time_point begin, end;
  double elapsed;

  std::string fname1 = "custom_output/data1.h5";
  std::string fname2 = "custom_output/data2.h5";
  H5::H5File hf1(fname1.c_str(), H5F_ACC_RDONLY);
  H5::H5File hf2(fname2.c_str(), H5F_ACC_RDONLY);

  // Read data from HDF5
  TIME_START("Reading data from HDF5...");
  ArrayDouble1D mesh_x = load_double_1D(hf2, "mesh_x");
  ArrayDouble1D mesh_y = load_double_1D(hf2, "mesh_y");
  ArrayDouble1D mesh_z = load_double_1D(hf2, "mesh_z");
  ArrayInt1D mesh_g = load_int_1D(hf2, "mesh_g");
  ArrayDouble2D angles = load_double_2D(hf2, "angles");
  ArrayDouble1D quadrature_weights = load_double_1D(hf2, "quadrature_weights");
  ArrayInt1D source_indices = load_int_1D(hf1, "source_indices");
  ArrayDouble1D source_strengths = load_double_1D(hf1, "source_strengths");
  ArrayDouble1D source_spectrum = load_double_1D(hf1, "source_spectrum");
  ArrayInt1D response_indices = load_int_1D(hf1, "response_indices");
  ArrayDouble1D response_strengths = load_double_1D(hf1, "response_strengths");
  ArrayDouble1D response_spectrum = load_double_1D(hf1, "response_spectrum");
  ArrayDouble2D mix_table = load_double_2D(hf1, "mix_table");
  ArrayInt3D material_map = load_int_3D(hf1, "material_map");
  ArrayDouble2D sigma_t = load_double_2D(hf1, "sigma_t");
  ArrayDouble3D sigma_s = load_double_3D(hf1, "sigma_s");
  ArrayDouble5D angular_flux_fwd = load_double_5D(hf2, "angular_flux_fwd");
  ArrayDouble5D angular_flux_adj = load_double_5D(hf2, "angular_flux_adj");
  TIME_END();

  // Dimensions
  int n_mix = mix_table.shape()[0];        // Number of mixed materials
  int nm = mix_table.shape()[1];           // Number of pure materials
  int nz = angular_flux_fwd.shape()[0];    // Number of Z intervals
  int ny = angular_flux_fwd.shape()[1];    // Number of Y intervals
  int nx = angular_flux_fwd.shape()[2];    // Number of X intervals
  int na = angular_flux_fwd.shape()[3];    // Number of angles
  int ngf = angular_flux_fwd.shape()[4];   // Number of energy groups in flux
  int ngx = sigma_t.shape()[1];            // Number of energy groups in XS
  int n_src = source_indices.shape()[0];   // Number of source indices
  int n_res = response_indices.shape()[0]; // Number of response indices

  // First energy group in flux
  int g0 = mesh_g[0];

  // Calculate reverse angle map
  TIME_START("Calculating reverse angle map...");
  ArrayInt1D reverse_angle_map{extents[na]};
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
        break;
      }
    }
  }
  TIME_END();

  // Calculate source and response
  TIME_START("Calculating source and response...");
  ArrayDouble4D source{extents[nz][ny][nx][ngf]};
  ArrayDouble4D response{extents[nz][ny][nx][ngf]};
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
  TIME_END();

  // Calculate scalar flux
  TIME_START("Calculating scalar flux...");
  ArrayDouble4D scalar_flux_fwd{extents[nz][ny][nx][ngf]};
  ArrayDouble4D scalar_flux_adj{extents[nz][ny][nx][ngf]};
  ArrayDouble4D scalar_flux_con{extents[nz][ny][nx][ngf]};
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
  TIME_END();

  // Calculate cross sections for mixed materials
  TIME_START("Calculating cross sections for mixed materials...");
  ArrayDouble2D sigma_t_mixed{extents[n_mix][ngx]};
  ArrayDouble3D sigma_s_mixed{extents[n_mix][ngx][ngx]};
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
  TIME_END();

  // Calculate perturbations in cross sections
  TIME_START("Calculating perturbations in cross sections...");
  ArrayDouble3D sigma_t_pert{extents[n_mix][nm][ngx]};
  ArrayDouble4D sigma_s_pert{extents[n_mix][nm][ngx][ngx]};
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
  TIME_END();

  // Calculate current
  TIME_START("Calculating current...");
  ArrayDouble5D current_fwd{extents[nz][ny][nx][ngf][3]};
  ArrayDouble5D current_adj{extents[nz][ny][nx][ngf][3]};
  ArrayDouble5D current_con{extents[nz][ny][nx][ngf][3]};
  for (int iz = 0; iz < nz; iz++) {                    // Z mesh index
    for (int iy = 0; iy < ny; iy++) {                  // Y mesh index
      for (int ix = 0; ix < nx; ix++) {                // X mesh index
        for (int ia = 0; ia < na; ia++) {              // Angle index
          double qw = quadrature_weights[ia] / FOURPI; // Quadrature weight
          for (int igf = 0; igf < ngf; igf++) {        // Energy index
            double affwd = angular_flux_fwd[iz][iy][ix][ia][igf] * qw;
            double afadj = angular_flux_adj[iz][iy][ix][ia][igf] * qw;
            for (int id = 0; id < 3; id++) { // Dimension index
              current_fwd[iz][iy][ix][igf][id] += affwd * angles[ia][id];
              current_adj[iz][iy][ix][igf][id] += afadj * angles[ia][id];
            }
          }
        }
      }
    }
  }
  for (int iz = 0; iz < nz; iz++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        for (int igf = 0; igf < ngf; igf++) {
          for (int id = 0; id < 3; id++) {
            current_con[iz][iy][ix][igf][id] =
                current_fwd[iz][iy][ix][igf][id] *
                current_adj[iz][iy][ix][igf][id];
          }
        }
      }
    }
  }
  TIME_END();

  // Calculate dR using angular flux
  TIME_START("Calculating dR...");
  ArrayDouble4D dR{extents[nm][nz][ny][nx]};
  for (int im = 0; im < nm; im++) {             // Material index
    for (int iz = 0; iz < nz; iz++) {           // Z mesh index
      for (int iy = 0; iy < ny; iy++) {         // Y mesh index
        for (int ix = 0; ix < nx; ix++) {       // X mesh index
          int i_mix = material_map[iz][iy][ix]; // Mixed material index
          // Total component of dR
          double dR_t = 0.0;
          for (int igf = 0; igf < ngf; igf++) { // Energy index (in flux data)
            int igx = igf + g0; // Energy index (in cross section data)
            // Scalar contributon flux
            double sfc = scalar_flux_con[iz][iy][ix][igf];
            // Perturbation in total cross section
            double dst = sigma_t_pert[i_mix][im][igx];
            dR_t += sfc * dst;
          }
          // Scattering component of dR
          double dR_s = 0.0;
          for (int igf = 0; igf < ngf; igf++) { // Energy index (in flux data)
            int igx = igf + g0; // Energy index (in cross section data)
            // Scalar adjoint flux (at E)
            double sfa = scalar_flux_adj[iz][iy][ix][igf];
            for (int jgf = 0; jgf < ngf; jgf++) { // E' index (in flux data)
              int jgx = jgf + g0; // E' index (in cross section data)
              // Perturbation in scattering cross section (E' -> E)
              double dss = sigma_s_pert[i_mix][im][igx][jgx];
              // Scalar forward flux (at E')
              double sff = scalar_flux_fwd[iz][iy][ix][jgf];
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
  auto int_type = H5::PredType::NATIVE_INT;
  auto double_type = H5::PredType::NATIVE_DOUBLE;
  std::string fnameo = "custom_output/data3.h5";
  H5::H5File hfo(fnameo, H5F_ACC_TRUNC);
  write_hdf5_array(hfo, "reverse_angle_map", reverse_angle_map, int_type);
  write_hdf5_array(hfo, "source", source, double_type);
  write_hdf5_array(hfo, "response", response, double_type);
  write_hdf5_array(hfo, "scalar_flux_fwd", scalar_flux_fwd, double_type);
  write_hdf5_array(hfo, "scalar_flux_adj", scalar_flux_adj, double_type);
  write_hdf5_array(hfo, "scalar_flux_con", scalar_flux_con, double_type);
  write_hdf5_array(hfo, "sigma_t_mixed", sigma_t_mixed, double_type);
  write_hdf5_array(hfo, "sigma_s_mixed", sigma_s_mixed, double_type);
  write_hdf5_array(hfo, "sigma_t_pert", sigma_t_pert, double_type);
  write_hdf5_array(hfo, "sigma_s_pert", sigma_s_pert, double_type);
  write_hdf5_array(hfo, "current_fwd", current_fwd, double_type);
  write_hdf5_array(hfo, "current_adj", current_adj, double_type);
  write_hdf5_array(hfo, "current_con", current_con, double_type);
  write_hdf5_array(hfo, "dR", dR, double_type);
  TIME_END();

  return 0;
}
