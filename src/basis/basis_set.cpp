#include "basis/basis_set.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "core/elements.hpp"

namespace ellidft {

double radial_gaussian_integral(int l, double zeta) {
  return std::tgamma(l + 1.5) / (2.0 * std::pow(zeta, l + 1.5));
}

std::vector<double> normalized_contraction(int l, const std::vector<double>& exponents,
                                           const std::vector<double>& coefficients) {
  if (exponents.size() != coefficients.size())
    throw std::invalid_argument("normalized_contraction: sizes differ");
  const std::size_t n = exponents.size();
  // Primitive k normalised: exp(-a_k r^2) / sqrt(I(l, 2 a_k)).
  std::vector<double> d(n);
  for (std::size_t k = 0; k < n; ++k)
    d[k] = coefficients[k] / std::sqrt(radial_gaussian_integral(l, 2.0 * exponents[k]));
  // Norm of the contraction: sum_kk' d_k d_k' I(l, a_k + a_k').
  double norm2 = 0.0;
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t m = 0; m < n; ++m)
      norm2 += d[k] * d[m] * radial_gaussian_integral(l, exponents[k] + exponents[m]);
  if (!(norm2 > 0.0)) throw std::invalid_argument("normalized_contraction: contraction has zero norm");
  for (double& x : d) x /= std::sqrt(norm2);
  return d;
}

BasisSet BasisSet::build(const Molecule& molecule, const BasisLibrary& library,
                         std::optional<bool> spherical) {
  BasisSet basis;
  basis.name_ = library.name();
  basis.spherical_ = spherical.value_or(library.spherical());

  for (std::size_t a = 0; a < molecule.atoms.size(); ++a) {
    const Atom& atom = molecule.atoms[a];
    if (!library.has_element(atom.Z))
      throw std::runtime_error("basis set '" + library.name() + "' has no functions for " +
                               element_symbol(atom.Z) + " (atom " + std::to_string(a + 1) + ")");
    for (const ShellDefinition& def : library.shells(atom.Z)) {
      Shell shell;
      shell.atom = static_cast<int>(a);
      shell.center = atom.r;
      shell.l = def.l;
      shell.spherical = basis.spherical_;
      shell.exponents = def.exponents;
      shell.coefficients = def.coefficients;
      shell.norm_coefficients = normalized_contraction(def.l, def.exponents, def.coefficients);

      basis.first_function_.push_back(basis.n_functions_);
      basis.n_functions_ += shell.n_functions();
      basis.max_l_ = std::max(basis.max_l_, shell.l);
      basis.shells_.push_back(std::move(shell));
    }
  }
  return basis;
}

BasisSet BasisSet::from_shells(std::string name, bool spherical, std::vector<Shell> shells) {
  BasisSet basis;
  basis.name_ = std::move(name);
  basis.spherical_ = spherical;
  for (Shell& shell : shells) {
    if (shell.l < 0 || shell.l > max_angular_momentum)
      throw std::invalid_argument("BasisSet::from_shells: unsupported angular momentum");
    shell.spherical = spherical;
    shell.norm_coefficients = normalized_contraction(shell.l, shell.exponents, shell.coefficients);
    basis.first_function_.push_back(basis.n_functions_);
    basis.n_functions_ += shell.n_functions();
    basis.max_l_ = std::max(basis.max_l_, shell.l);
    basis.shells_.push_back(std::move(shell));
  }
  return basis;
}

std::size_t BasisSet::n_primitives() const {
  std::size_t n = 0;
  for (const Shell& s : shells_) n += s.n_primitives();
  return n;
}

}  // namespace ellidft
