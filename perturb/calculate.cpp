#include "calculate.hpp"

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
}

void PADVANTG::calculate_sigma_mixed() {
  TIME_START("Calculating cross sections for mixed materials...");
  sigma_t_mixed.resize(extents[nmix][ngx]);
  sigma_s_mixed.resize(extents[nmix][ngx][ngx]);
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
  sigma_t_pert.resize(extents[nmix][nm][ngx]);
  sigma_s_pert.resize(extents[nmix][nm][ngx][ngx]);
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

void PADVANTG::calculate_reverse_angle_map() {
  TIME_START("Calculating reverse angle map...");
  reverse_angle_map.resize(extents[na]);
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

void PADVANTG::allocate_forward_arrays() {
  TIME_START("Allocating forward arrays...");
  source_fwd.resize(extents[ngx][nz][ny][nx]);
  scalar_flux_fwd.resize(extents[ngx][nz][ny][nx]);
  current_fwd.resize(extents[ngx][nz][ny][nx][3]);
  TIME_END();
}

void PADVANTG::allocate_adjoint_arrays() {
  TIME_START("Allocating adjoint arrays...");
  source_adj.resize(extents[nas][ngx][nz][ny][nx]);
  scalar_flux_adj.resize(extents[nas][ngx][nz][ny][nx]);
  scalar_flux_con.resize(extents[nas][ngx][nz][ny][nx]);
  current_adj.resize(extents[nas][ngx][nz][ny][nx][3]);
  current_con.resize(extents[nas][ngx][nz][ny][nx][3]);
  response.resize(extents[nas]);
  dR.resize(extents[nas][nm][nz][ny][nx]);
  TIME_END();
}

void PADVANTG::calculate_forward_source() {
  TIME_START("Calculating forward 4D source...");
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
  calculate_sigma_mixed();
  calculate_sigma_pert();
  allocate_forward_arrays();
  calculate_forward_source();
  calculate_forward_scalar_flux();
  calculate_forward_current();
  calculate_reverse_angle_map();
  allocate_adjoint_arrays();
  for (int i = 0; i < nas; i++) {
    read_adjoint_data(i); // ngfa, g0a
    calculate_adjoint_source(i);
    calculate_adjoint_scalar_flux(i);
    calculate_adjoint_current(i);
    calculate_response(i);
    calculate_dR(i);
  }
  write_all_data();
}

int main() {
  PADVANTG p = PADVANTG();
  p.run_all();
  return 0;
}
