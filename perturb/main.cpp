#include <boost/program_options.hpp>

#include "CalculateDR.hpp"

namespace po = boost::program_options;

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
