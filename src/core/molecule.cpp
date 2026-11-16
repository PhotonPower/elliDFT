#include "core/molecule.hpp"

#include <cmath>

namespace ellidft {

int Molecule::nuclear_charge() const {
  int z = 0;
  for (const Atom& a : atoms) z += a.Z;
  return z;
}

int Molecule::n_electrons() const { return nuclear_charge() - charge; }

double Molecule::nuclear_repulsion_energy() const {
  double energy = 0.0;
  for (std::size_t a = 0; a < atoms.size(); ++a)
    for (std::size_t b = 0; b < a; ++b) {
      const double dx = atoms[a].r[0] - atoms[b].r[0];
      const double dy = atoms[a].r[1] - atoms[b].r[1];
      const double dz = atoms[a].r[2] - atoms[b].r[2];
      energy += atoms[a].Z * atoms[b].Z / std::sqrt(dx * dx + dy * dy + dz * dz);
    }
  return energy;
}

}  // namespace ellidft
