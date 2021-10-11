#include "H5Cpp.h"
#include "boost/multi_array.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#define THREAD_NUM omp_get_thread_num()
#else
#define THREAD_NUM 0
#endif

#define PI 3.14159265358979323846
#define FOURPI (4.0 * PI)

#define TIME_START(s)                                                          \
  {                                                                            \
    printf("%-52s", s);                                                        \
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
using std::size_t;
using std::string;
using std::to_string;
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

template <typename T, size_t NDIMS>
void read_hdf5_array(H5::H5File hf, const string &dname,
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

template <typename T, size_t NDIMS>
void write_hdf5_array(H5::H5File hf, const string &dname,
                      const multi_array<T, NDIMS> &A) {
  /* Write a boost::multi_array to HDF5.
   */
  H5::DataType dtype = get_hdf5_data_type<T>::type;
  vector<hsize_t> dims(A.shape(), A.shape() + NDIMS);
  H5::DataSpace dspace(NDIMS, dims.data());
  H5::DataSet dset = hf.createDataSet(dname, dtype, dspace);
  dset.write(A.data(), dtype);
}

class PADVANTG {
  // Flag to calculate current
  bool calculate_current;
  // Flag to calculate and write source, flux, and current data to file
  bool write_more;

  // Variables for timing
  high_resolution_clock::time_point begin, end;
  double elapsed;

  // Dimensions
  int nas;  // Number of adjoint sources
  int nm;   // Number of pure materials
  int nmix; // Number of mixed materials
  int ngx;  // Number of energy groups in XS
  int ngff; // Number of energy groups in forward flux
  int ngfa; // Number of energy groups in adjoint flux
  int nz;   // Number of Z intervals
  int ny;   // Number of Y intervals
  int nx;   // Number of X intervals
  int na;   // Number of angles
  int g0f;  // First energy group in forward flux
  int g0a;  // First energy group in adjoint flux

  // Tally IDs
  vector<int> tallies;

  // Data arrays (input)
  multi_array<double, 2> sigma_t;             // Total cross sections
  multi_array<double, 3> sigma_s;             // Scattering cross sections
  multi_array<MixData, 1> mixtable;           // Mix table
  multi_array<int, 3> matids;                 // Material map
  multi_array<double, 2> quadrature_angles;   // Quadrature angles
  multi_array<double, 1> quadrature_weights;  // Quadrature weights
  multi_array<double, 2> source_spectra_fwd;  // Forward source spectra
  multi_array<double, 3> source_strength_fwd; // Forward source strength
  multi_array<int, 1> mesh_g_fwd;             // Forward energy bin indices
  multi_array<double, 5> angular_flux_fwd;    // Forward angular flux
  multi_array<int, 1> mesh_g_adj;             // Adjoint source spectra
  multi_array<double, 2> source_spectra_adj;  // Adjoint source strength
  multi_array<double, 3> source_strength_adj; // Adjoint energy bin indices
  multi_array<double, 5> angular_flux_adj;    // Adjoint angular flux

  // Data arrays (output)
  multi_array<int, 1> reverse_angle_map; // Reverse angle map
  multi_array<double, 2> sigma_t_mixed;  // Total cross sections (mixed)
  multi_array<double, 3> sigma_s_mixed;  // Scattering cross sections (mixed)
  multi_array<double, 3> sigma_t_pert;   // Total cross sections (perturbed)
  multi_array<double, 4> sigma_s_pert; // Scattering cross sections (perturbed)
  multi_array<double, 4> source_fwd;   // 4D forward source
  multi_array<double, 5> source_adj;   // 4D adjoint source
  multi_array<double, 4> scalar_flux_fwd; // Forward scalar flux
  multi_array<double, 5> scalar_flux_adj; // Adjoint scalar flux
  multi_array<double, 5> scalar_flux_con; // Contributon scalar flux
  multi_array<double, 5> current_fwd;     // Forward current
  multi_array<double, 6> current_adj;     // Adjoint current
  multi_array<double, 6> current_con;     // Contributon current
  multi_array<double, 1> response;        // Response
  multi_array<double, 5> dR;              // Perturbation in response

  // Functions
  void read_tally_ids();
  void read_xs_data();
  void read_forward_data();
  void read_adjoint_data(int);
  void write_all_data();
  void allocate_xs_arrays();
  void allocate_forward_arrays();
  void allocate_adjoint_arrays();
  void calculate_reverse_angle_map();
  void calculate_sigma_mixed();
  void calculate_sigma_pert();
  void calculate_forward_source();
  void calculate_adjoint_source(int);
  void calculate_forward_scalar_flux();
  void calculate_adjoint_scalar_flux(int);
  void calculate_forward_current();
  void calculate_adjoint_current(int);
  void calculate_response(int);
  void calculate_dR(int);

public:
  PADVANTG(bool);
  void run_all();
};
