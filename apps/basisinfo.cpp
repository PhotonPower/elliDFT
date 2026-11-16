// ellidft-basisinfo: prints the size of the basis of a molecule as JSON.
//
// Usage: ellidft-basisinfo <geometry.xyz> <basis name> [--spherical|--cartesian] [--bohr]

#include <cstdlib>
#include <exception>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "basis/basis_set.hpp"
#include "core/elements.hpp"
#include "core/xyz.hpp"

int main(int argc, char** argv) {
  using namespace ellidft;

  std::vector<std::string> positional;
  std::optional<bool> spherical;
  LengthUnit unit = LengthUnit::Angstrom;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--spherical")
      spherical = true;
    else if (arg == "--cartesian")
      spherical = false;
    else if (arg == "--bohr")
      unit = LengthUnit::Bohr;
    else if (!arg.empty() && arg[0] == '-') {
      std::cerr << "unknown option " << arg << "\n";
      return EXIT_FAILURE;
    } else
      positional.push_back(arg);
  }
  if (positional.size() != 2) {
    std::cerr << "usage: ellidft-basisinfo <geometry.xyz> <basis name> "
                 "[--spherical|--cartesian] [--bohr]\n";
    return EXIT_FAILURE;
  }

  try {
    const Molecule molecule = read_xyz_file(positional[0], unit);
    const BasisLibrary library = BasisLibrary::load(positional[1]);
    const BasisSet basis = BasisSet::build(molecule, library, spherical);

    std::vector<int> shells_per_l(static_cast<std::size_t>(basis.max_l() + 1), 0);
    for (const Shell& s : basis.shells()) ++shells_per_l[static_cast<std::size_t>(s.l)];

    nlohmann::ordered_json symbols = nlohmann::ordered_json::array();
    nlohmann::ordered_json geometry = nlohmann::ordered_json::array();
    for (const Atom& a : molecule.atoms) {
      symbols.push_back(element_symbol(a.Z));
      geometry.push_back({a.r[0], a.r[1], a.r[2]});
    }

    nlohmann::ordered_json out;
    out["molecule"] = {{"name", molecule.name},
                       {"charge", molecule.charge},
                       {"n_atoms", molecule.n_atoms()},
                       {"n_electrons", molecule.n_electrons()},
                       {"symbols", symbols},
                       {"geometry_bohr", geometry}};
    out["basis"] = {{"name", basis.name()},
                    {"spherical", basis.spherical()},
                    {"n_shells", basis.n_shells()},
                    {"n_shells_per_l", shells_per_l},
                    {"n_primitives", basis.n_primitives()},
                    {"n_functions", basis.n_functions()}};
    std::cout << out.dump(2) << "\n";
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
