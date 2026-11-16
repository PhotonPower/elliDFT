#include "xc/functional.hpp"

#include "xc/xc_evaluator.hpp"

#include <cctype>
#include <stdexcept>

namespace ellidft {

Functional Functional::by_name(std::string_view name) {
  std::string key;
  for (char c : name) key += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));

  Functional f;
  if (key == "HF") {
    f.name = "HF";
    f.exact_exchange = 1.0;
  } else if (key == "LDA") {
    // VWN5 is the variant libxc calls LDA_C_VWN (not VWN-RPA).
    f.name = "LDA";
    f.components = {{"LDA_X", 1.0}, {"LDA_C_VWN", 1.0}};
  } else if (key == "PBE") {
    f.name = "PBE";
    f.components = {{"GGA_X_PBE", 1.0}, {"GGA_C_PBE", 1.0}};
    f.needs_gradient = true;
  } else if (key == "PBE0" || key == "B3LYP" || key == "B3LYP5" || key == "CAM-B3LYP") {
    const char* libxc_name = key == "PBE0"     ? "HYB_GGA_XC_PBEH"
                             : key == "B3LYP"  ? "HYB_GGA_XC_B3LYP"
                             : key == "B3LYP5" ? "HYB_GGA_XC_B3LYP5"
                                               : "HYB_GGA_XC_CAM_B3LYP";
    f.name = key;
    f.components = {{libxc_name, 1.0}};
    const HybridParameters parameters = libxc_hybrid_parameters(libxc_name);
    f.needs_gradient = parameters.needs_gradient;
    f.exact_exchange = parameters.exact_exchange;
    f.long_range_exchange = parameters.long_range_exchange;
    f.omega = parameters.omega;
  } else {
    throw std::invalid_argument("unknown functional '" + std::string(name) + "'");
  }
  return f;
}

}  // namespace ellidft
