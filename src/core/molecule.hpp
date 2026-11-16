#pragma once

#include <array>
#include <string>
#include <vector>

namespace ellidft {

using Vec3 = std::array<double, 3>;

struct Atom {
  int Z = 0;   // atomic number
  Vec3 r{};    // position in Bohr
};

// A molecule: nuclei with positions in Bohr and a total charge.
// Only closed-shell systems are in scope, so there is no spin state.
struct Molecule {
  std::string name;
  int charge = 0;
  std::vector<Atom> atoms;

  std::size_t n_atoms() const { return atoms.size(); }

  // Sum of nuclear charges.
  int nuclear_charge() const;

  // Number of electrons = nuclear charge - total charge.
  int n_electrons() const;

  // Coulomb repulsion of the point nuclei, sum_{A<B} Z_A Z_B / |R_A - R_B|,
  // in Hartree.
  double nuclear_repulsion_energy() const;
};

}  // namespace ellidft
