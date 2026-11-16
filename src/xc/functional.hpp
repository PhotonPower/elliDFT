#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace ellidft {

// An exchange-correlation method as a data object (specification, section
// 4): the density functional part, the fraction of exact exchange and the
// range separation parameters. SCF and response theory read these values
// from the same place.
struct Functional {
  // One term of the density functional part: a functional of libxc, named by
  // its libxc identifier (e.g. "LDA_X", "GGA_C_PBE"), times a coefficient.
  struct Component {
    std::string libxc_name;
    double coefficient = 1.0;
  };

  std::string name;
  std::vector<Component> components;
  bool needs_gradient = false;        // true for GGA parts

  // Exact exchange: -1/2 * (exact_exchange * K + long_range_exchange * K_LR)
  // in the Fock matrix, K_LR with the operator erf(omega r)/r.
  double exact_exchange = 0.0;
  double long_range_exchange = 0.0;
  double omega = 0.0;                 // 1/Bohr

  bool has_density_part() const { return !components.empty(); }

  // Functionals known by name (case-insensitive):
  //   "HF"   exact exchange only
  //   "LDA"  Slater exchange + VWN5 correlation   (LDA_X + LDA_C_VWN)
  //   "PBE"  PBE exchange + PBE correlation       (GGA_X_PBE + GGA_C_PBE)
  //   "PBE0"       HYB_GGA_XC_PBEH, 25 % exact exchange
  //   "B3LYP"      HYB_GGA_XC_B3LYP, 20 % exact exchange, correlation with
  //                VWN-RPA (the variant of the Gaussian program)
  //   "B3LYP5"     HYB_GGA_XC_B3LYP5, the same with VWN5
  //   "CAM-B3LYP"  HYB_GGA_XC_CAM_B3LYP, range separated: 0.19 K + 0.46 K_LR,
  //                omega = 0.33 / Bohr
  // The exact-exchange parameters of the hybrids are taken from libxc.
  // Throws std::invalid_argument for unknown names.
  static Functional by_name(std::string_view name);
};

}  // namespace ellidft
