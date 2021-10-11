#include <boost/program_options.hpp>

#include "calculate.hpp"

namespace po = boost::program_options;

PADVANTG::PADVANTG(bool _write_more) {
  calculate_current = _write_more;
  write_more = _write_more;
}

void PADVANTG::read_tally_ids() {
  TIME_START("Reading tally IDs from file...");
  string fname = "tally_list.txt";
  std::ifstream infile(fname);
  string buffer;
  while (std::getline(infile, buffer)) {
    int tally = std::stof(buffer);
    tallies.push_back(tally);
  }
  nas = tallies.size();
  TIME_END();
}

void PADVANTG::read_xs_data() {
  TIME_START("Reading cross section data from HDF5...");
  H5::H5File hf_xs("hdf5/advantg_xs.h5", H5F_ACC_RDONLY);
  read_hdf5_array(hf_xs, "sigma_t", sigma_t);
  read_hdf5_array(hf_xs, "sigma_s", sigma_s);
  nm = sigma_t.shape()[0];
  ngx = sigma_t.shape()[1];
  TIME_END();
}
void PADVANTG::read_forward_data() {
  TIME_START("Reading forward flux data from HDF5...");
  H5::H5File hf_fi("hdf5/advantg_fwd_inp.h5", H5F_ACC_RDONLY);
  H5::H5File hf_fo("hdf5/advantg_fwd_out.h5", H5F_ACC_RDONLY);
  read_hdf5_array(hf_fi, "mixtable", mixtable);
  read_hdf5_array(hf_fi, "matids", matids);
  read_hdf5_array(hf_fo, "denovo/quadrature_angles", quadrature_angles);
  read_hdf5_array(hf_fo, "denovo/quadrature_weights", quadrature_weights);
  read_hdf5_array(hf_fi, "volsrc/spectra", source_spectra_fwd);
  read_hdf5_array(hf_fi, "volsrc/strength", source_strength_fwd);
  read_hdf5_array(hf_fo, "denovo/mesh_g", mesh_g_fwd);
  read_hdf5_array(hf_fo, "denovo/angular_flux", angular_flux_fwd);
  nmix = mixtable[mixtable.shape()[0] - 1].row + 1;
  ngff = angular_flux_fwd.shape()[0];
  nz = angular_flux_fwd.shape()[1];
  ny = angular_flux_fwd.shape()[2];
  nx = angular_flux_fwd.shape()[3];
  na = angular_flux_fwd.shape()[4];
  g0f = mesh_g_fwd[0];
  TIME_END();
}

void PADVANTG::read_adjoint_data(int ias) {
  TIME_START(("Reading adjoint flux data from HDF5 for tally " +
              to_string(tallies[ias]) + "...")
                 .c_str());
  int tally = tallies[ias];
  H5::H5File hf_ai(("hdf5/advantg_adj_" + to_string(tally) + "_inp.h5").c_str(),
                   H5F_ACC_RDONLY);
  H5::H5File hf_ao(("hdf5/advantg_adj_" + to_string(tally) + "_out.h5").c_str(),
                   H5F_ACC_RDONLY);
  read_hdf5_array(hf_ai, "volsrc/spectra", source_spectra_adj);
  read_hdf5_array(hf_ai, "volsrc/strength", source_strength_adj);
  read_hdf5_array(hf_ao, "denovo/mesh_g", mesh_g_adj);
  read_hdf5_array(hf_ao, "denovo/angular_flux", angular_flux_adj);
  ngfa = angular_flux_adj.shape()[0];
  g0a = mesh_g_adj[0];
  TIME_END();
}

void PADVANTG::write_all_data() {
  TIME_START("Writing all data to HDF5...");
  // H5::H5File hf_xs("hdf5/xs.h5", H5F_ACC_TRUNC);
  // write_hdf5_array(hf_xs, "sigma_t_mixed", sigma_t_mixed);
  // write_hdf5_array(hf_xs, "sigma_s_mixed", sigma_s_mixed);
  // write_hdf5_array(hf_xs, "sigma_t_pert", sigma_t_pert);
  // write_hdf5_array(hf_xs, "sigma_s_pert", sigma_s_pert);
  H5::H5File hf_main("hdf5/main.h5", H5F_ACC_TRUNC);
  write_hdf5_array(hf_main, "response", response);
  write_hdf5_array(hf_main, "dR", dR);
  if (write_more) {
    H5::H5File hf_plot("hdf5/more_data.h5", H5F_ACC_TRUNC);
    write_hdf5_array(hf_plot, "source_fwd", source_fwd);
    write_hdf5_array(hf_plot, "source_adj", source_adj);
    write_hdf5_array(hf_plot, "scalar_flux_fwd", scalar_flux_fwd);
    write_hdf5_array(hf_plot, "scalar_flux_adj", scalar_flux_adj);
    write_hdf5_array(hf_plot, "scalar_flux_con", scalar_flux_con);
    if (calculate_current) {
      write_hdf5_array(hf_plot, "current_fwd", current_fwd);
      write_hdf5_array(hf_plot, "current_adj", current_adj);
      write_hdf5_array(hf_plot, "current_con", current_con);
    }
  }
  TIME_END();
}

void PADVANTG::allocate_xs_arrays() {
  TIME_START("Allocating cross section arrays...");
  sigma_t_mixed.resize(extents[nmix][ngx]);
  sigma_s_mixed.resize(extents[nmix][ngx][ngx]);
  sigma_t_pert.resize(extents[nmix][nm][ngx]);
  sigma_s_pert.resize(extents[nmix][nm][ngx][ngx]);
  TIME_END();
}

void PADVANTG::allocate_forward_arrays() {
  TIME_START("Allocating forward arrays...");
  source_fwd.resize(extents[ngx][nz][ny][nx]);
  scalar_flux_fwd.resize(extents[ngx][nz][ny][nx]);
  if (calculate_current) {
    current_fwd.resize(extents[ngx][nz][ny][nx][3]);
  }
  TIME_END();
}

void PADVANTG::allocate_adjoint_arrays() {
  TIME_START("Allocating adjoint arrays...");
  source_adj.resize(extents[nas][ngx][nz][ny][nx]);
  scalar_flux_adj.resize(extents[nas][ngx][nz][ny][nx]);
  scalar_flux_con.resize(extents[nas][ngx][nz][ny][nx]);
  if (calculate_current) {
    current_adj.resize(extents[nas][ngx][nz][ny][nx][3]);
    current_con.resize(extents[nas][ngx][nz][ny][nx][3]);
  }
  response.resize(extents[nas]);
  dR.resize(extents[nas][nm][nz][ny][nx]);
  TIME_END();
}

void PADVANTG::calculate_reverse_angle_map() {
  TIME_START("Calculating reverse angle map...");
  reverse_angle_map.resize(extents[na]);
#pragma omp parallel for schedule(guided)
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
}

void PADVANTG::calculate_sigma_mixed() {
  TIME_START("Calculating cross sections for mixed materials...");
#pragma omp parallel for schedule(guided)
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
}

void PADVANTG::calculate_sigma_pert() {
  TIME_START("Calculating perturbations in cross sections...");
#pragma omp parallel for schedule(guided)
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
}

void PADVANTG::calculate_forward_source() {
  TIME_START("Calculating forward 4D source...");
#pragma omp parallel for schedule(guided)
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
  TIME_END();
}

void PADVANTG::calculate_adjoint_source(int ias) {
  TIME_START(("Calculating adjoint 4D source for tally " +
              to_string(tallies[ias]) + "...")
                 .c_str());
#pragma omp parallel for schedule(guided)
  for (int igx = 0; igx < ngx; igx++) { // Energy index
    double spec_adj = source_spectra_adj[0][igx];
    for (int iz = 0; iz < nz; iz++) {     // Z mesh index
      for (int iy = 0; iy < ny; iy++) {   // Y mesh index
        for (int ix = 0; ix < nx; ix++) { // X mesh index
          source_adj[ias][igx][iz][iy][ix] =
              source_strength_adj[iz][iy][ix] * spec_adj;
        }
      }
    }
  }
  TIME_END();
}

void PADVANTG::calculate_forward_scalar_flux() {
  TIME_START("Calculating forward scalar flux...");
#pragma omp parallel for schedule(guided)
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
  TIME_END();
}

void PADVANTG::calculate_adjoint_scalar_flux(int ias) {
  TIME_START(("Calculating adjoint scalar flux for tally " +
              to_string(tallies[ias]) + "...")
                 .c_str());
#pragma omp parallel for schedule(guided)
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
            if (igx - g0a >= 0 && igx - g0a < ngfa) {
              fadj = angular_flux_adj[igx - g0a][iz][iy][ix][ja] * qw;
            }
            scalar_flux_adj[ias][igx][iz][iy][ix] += fadj;
            scalar_flux_con[ias][igx][iz][iy][ix] += ffwd * fadj;
          }
        }
      }
    }
  }
  TIME_END();
}

void PADVANTG::calculate_forward_current() {
  TIME_START("Calculating forward current...");
#pragma omp parallel for schedule(guided)
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
  TIME_END();
}

void PADVANTG::calculate_adjoint_current(int ias) {
  TIME_START(("Calculating adjoint current for tally " +
              to_string(tallies[ias]) + "...")
                 .c_str());
#pragma omp parallel for schedule(guided)
  for (int igx = 0; igx < ngx; igx++) {     // Energy index
    for (int iz = 0; iz < nz; iz++) {       // Z mesh index
      for (int iy = 0; iy < ny; iy++) {     // Y mesh index
        for (int ix = 0; ix < nx; ix++) {   // X mesh index
          for (int ia = 0; ia < na; ia++) { // Angle index
            double qw = quadrature_weights[ia] / FOURPI;
            double afadj = 0.0;
            if (igx - g0a >= 0 && igx - g0a < ngfa) {
              afadj = angular_flux_adj[igx - g0a][iz][iy][ix][ia] * qw;
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
  TIME_END();
}

void PADVANTG::calculate_response(int ias) {
  TIME_START(
      ("Calculating response for tally " + to_string(tallies[ias]) + "...")
          .c_str());
#pragma omp parallel for schedule(guided)
  for (int igx = 0; igx < ngx; igx++) {   // Energy index
    for (int iz = 0; iz < nz; iz++) {     // Z mesh index
      for (int iy = 0; iy < ny; iy++) {   // Y mesh index
        for (int ix = 0; ix < nx; ix++) { // X mesh index
          double sff = scalar_flux_fwd[igx][iz][iy][ix];
          double sa = source_adj[ias][igx][iz][iy][ix];
          response[ias] += sff * sa;
        }
      }
    }
  }
  TIME_END();
}

void PADVANTG::calculate_dR(int ias) {
  TIME_START(
      ("Calculating dR for tally " + to_string(tallies[ias]) + "...").c_str());
#pragma omp parallel for schedule(guided)
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
  TIME_END();
}

void PADVANTG::run_all() {
  // Load, calculate, and write all data
  read_tally_ids();
  read_xs_data();
  read_forward_data();
  allocate_xs_arrays();
  allocate_forward_arrays();
  calculate_reverse_angle_map();
  calculate_sigma_mixed();
  calculate_sigma_pert();
  calculate_forward_source();
  calculate_forward_scalar_flux();
  if (calculate_current) {
    calculate_forward_current();
  }
  allocate_adjoint_arrays();
  for (int i = 0; i < nas; i++) {
    read_adjoint_data(i);
    calculate_adjoint_source(i);
    calculate_adjoint_scalar_flux(i);
    if (calculate_current) {
      calculate_adjoint_current(i);
    }
    calculate_response(i);
    calculate_dR(i);
  }
  write_all_data();
}

po::variables_map parse_args(int argc, char **argv) {
  // Define command line options
  po::variables_map vm;
  po::options_description desc("Allowed Options");
  desc.add_options()("help,h", po::bool_switch()->default_value(false),
                     "Display this information.")(
      "write_more,w", po::bool_switch()->default_value(false),
      "Write source, flux, and current data to HDF5.")(
      "num_threads,j", po::value<int>()->default_value(1),
      "Number of threads to use.");

  // Parse arguments and look for invalid options
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (std::exception &e) {
    std::cout << "Error: " << e.what() << std::endl;
    std::cout << desc << std::endl;
    exit(1);
  }

  // Display help message if requested
  if (vm["help"].as<bool>()) {
    std::cout << desc << std::endl;
    exit(0);
  }

  // Set number of OMP threads and error-check
  int num_threads = vm["num_threads"].as<int>();
#ifdef _OPENMP
  if (num_threads < 1) {
    std::cout << "Error: number of threads must be positive." << std::endl;
    exit(1);
  }
  omp_set_num_threads(num_threads);
#else
  if (num_threads != 1) {
    std::cout << "Warning: the program was not compiled with OMP support but "
                 "multiple threads were requested in the command line. "
                 "Calculation will proceed with a single thread."
              << std::endl;
  }
#endif

  // Return options
  return vm;
}

int main(int argc, char **argv) {
  po::variables_map vm = parse_args(argc, argv);
  PADVANTG p = PADVANTG(vm["write_more"].as<bool>());
  p.run_all();
  return 0;
}
