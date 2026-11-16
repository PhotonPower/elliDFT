#include "ellipsometry/stack.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ellidft::ellipsometry {

namespace {

// A medium of the stack for one polarisation: q = normal component of the
// wave vector in units of 2 pi / lambda, eta = ratio of the two tangential
// field amplitudes that are continuous across interfaces.
//
// s polarisation (electric field in the plane of the film): only the
// ordinary index acts, q = sqrt(N_o^2 - kx^2), eta = q.
//
// p polarisation, optic axis along the normal z: the extraordinary wave obeys
//   kx^2 / N_e^2 + kz^2 / N_o^2 = 1          (in units of (2 pi / lambda)^2;
//   Wikipedia, "Birefringence", wave-vector surface of a uniaxial material),
// hence q = (N_o / N_e) sqrt(N_e^2 - kx^2). From Maxwell's equations for
// H = H_y the tangential electric field is E_x = q H_y / N_o^2 (the in-plane
// permittivity is N_o^2), so eta = q / N_o^2. For N_e = N_o both reduce to
// the isotropic case. (The field relation is derived here.)
//
// The root with non-negative imaginary part is the wave that decays into the
// medium (time dependence exp(-i omega t), N = n + i k).
struct Medium {
  Complex q_s, eta_s, q_p, eta_p;
};

Complex decaying_root(Complex z) {
  Complex q = std::sqrt(z);
  if (q.imag() < 0.0 || (q.imag() == 0.0 && q.real() < 0.0)) q = -q;
  return q;
}

Medium medium(Complex N_o, Complex N_e, Complex kx) {
  Medium m;
  m.q_s = decaying_root(N_o * N_o - kx * kx);
  m.eta_s = m.q_s;
  if (N_e == N_o) {
    m.q_p = m.q_s;
  } else {
    // q^2 = N_o^2 (N_e^2 - kx^2) / N_e^2
    m.q_p = decaying_root(N_o * N_o * (N_e * N_e - kx * kx) / (N_e * N_e));
  }
  m.eta_p = m.q_p / (N_o * N_o);
  return m;
}

// Fresnel coefficient between media i (incident side) and j:
//   r = (eta_i - eta_j) / (eta_i + eta_j).
// For isotropic media this is r_s = (n_i cos_i - n_j cos_j) / (n_i cos_i + n_j cos_j) and
// r_p = (n_j cos_i - n_i cos_j) / (n_j cos_i + n_i cos_j) of Wikipedia, "Fresnel equations".
Complex fresnel(Complex eta_i, Complex eta_j) { return (eta_i - eta_j) / (eta_i + eta_j); }

struct Media {
  std::vector<Medium> m;          // ambient, layers, substrate
  std::vector<double> thickness;  // of the layers
  Complex kx;
};

Media media_of(const Stack& stack, double wavelength_nm, double angle_deg) {
  if (!stack.ambient || !stack.substrate) throw std::invalid_argument("reflect: ambient and substrate are required");
  if (!(wavelength_nm > 0.0)) throw std::invalid_argument("reflect: wavelength must be positive");
  if (angle_deg < 0.0 || angle_deg >= 90.0) throw std::invalid_argument("reflect: angle of incidence must be in [0, 90)");
  Media media;
  const Complex N0 = stack.ambient->index(wavelength_nm);
  media.kx = N0 * std::sin(angle_deg * M_PI / 180.0);
  media.m.push_back(medium(N0, N0, media.kx));
  for (const Layer& layer : stack.layers) {
    if (!layer.material) throw std::invalid_argument("reflect: layer without material");
    if (layer.thickness_nm < 0.0) throw std::invalid_argument("reflect: negative layer thickness");
    const Complex N_o = layer.material->index(wavelength_nm);
    const Complex N_e = layer.extraordinary ? layer.extraordinary->index(wavelength_nm) : N_o;
    media.m.push_back(medium(N_o, N_e, media.kx));
    media.thickness.push_back(layer.thickness_nm);
  }
  const Complex Ns = stack.substrate->index(wavelength_nm);
  media.m.push_back(medium(Ns, Ns, media.kx));
  return media;
}

// Reflection and transmission of the layers between the first and the last
// medium of (eta, q) with characteristic matrices. In every medium the
// continuous tangential field is U = a e^{i b} + b e^{-i b}, its partner
// V = eta (a e^{i b} - b e^{-i b}); across a layer of phase thickness
// beta = 2 pi d q / lambda
//   U(top) = cos(beta) U(bottom) - (i / eta) sin(beta) V(bottom)
//   V(top) = -i eta sin(beta) U(bottom) + cos(beta) V(bottom).
// Below there is only the transmitted wave, (U, V) = t (1, eta_last); above
// U = 1 + r, V = eta_first (1 - r). t is the ratio of the U amplitudes; the
// product of t for the two directions is the same for either choice of the
// field U.
struct Amplitudes {
  Complex r, t;
};
Amplitudes through(const std::vector<Complex>& eta, const std::vector<Complex>& q, const std::vector<double>& thickness,
                   double wavelength_nm) {
  Complex U = 1.0, V = eta.back();
  for (std::size_t j = eta.size() - 2; j >= 1; --j) {
    const Complex beta = 2.0 * M_PI * thickness[j - 1] * q[j] / wavelength_nm;
    const Complex c = std::cos(beta), s = std::sin(beta), i(0.0, 1.0);
    const Complex U_top = c * U - i / eta[j] * s * V;
    const Complex V_top = -i * eta[j] * s * U + c * V;
    U = U_top;
    V = V_top;
  }
  const Complex denominator = eta.front() * U + V;
  return {(eta.front() * U - V) / denominator, 2.0 * eta.front() / denominator};
}

}  // namespace

Reflection reflect(const Stack& stack, double wavelength_nm, double angle_deg) {
  const Media media = media_of(stack, wavelength_nm, angle_deg);
  const std::size_t n_layers = media.thickness.size();
  const std::vector<Medium>& m = media.m;

  // A layer between an upper interface with coefficient r_u and everything
  // below it with total coefficient R: summing the multiply reflected beams
  // (a geometric series, with t_u t_u' = 1 - r_u^2 and r_u' = -r_u for the
  // reversed interface) gives
  //   R_total = (r_u + R e^{2 i beta}) / (1 + r_u R e^{2 i beta}),
  //   beta = 2 pi d q / lambda.
  // Applied from the substrate upwards this covers any number of layers.
  // The unit tests compare this with the explicitly summed series.
  Complex Rs = fresnel(m[n_layers].eta_s, m[n_layers + 1].eta_s);
  Complex Rp = fresnel(m[n_layers].eta_p, m[n_layers + 1].eta_p);
  for (std::size_t j = n_layers; j >= 1; --j) {
    const double k0d = 2.0 * M_PI * media.thickness[j - 1] / wavelength_nm;
    const Complex phase_s = std::exp(Complex(0.0, 2.0) * (k0d * m[j].q_s));
    const Complex phase_p = std::exp(Complex(0.0, 2.0) * (k0d * m[j].q_p));
    const Complex rs = fresnel(m[j - 1].eta_s, m[j].eta_s);
    const Complex rp = fresnel(m[j - 1].eta_p, m[j].eta_p);
    Rs = (rs + Rs * phase_s) / (1.0 + rs * Rs * phase_s);
    Rp = (rp + Rp * phase_p) / (1.0 + rp * Rp * phase_p);
  }
  return {Rs, Rp};
}

Observation observe(const Stack& stack, double wavelength_nm, double angle_deg) {
  const Reflection front = reflect(stack, wavelength_nm, angle_deg);
  Complex I_ps = front.rp * std::conj(front.rs);
  double I_ss = std::norm(front.rs), I_pp = std::norm(front.rp);

  if (stack.substrate_thickness_nm > 0.0) {
    const Media media = media_of(stack, wavelength_nm, angle_deg);
    const std::size_t n = media.m.size();
    const Complex Nb = (stack.backside ? stack.backside : stack.ambient)->index(wavelength_nm);
    const Medium back = medium(Nb, Nb, media.kx);
    const Medium& substrate = media.m.back();
    // intensity factor of one pass through the plate (zero for an opaque substrate)
    const double pass = std::exp(-4.0 * M_PI * stack.substrate_thickness_nm * substrate.q_s.imag() / wavelength_nm);
    if (pass > 0.0) {
      Complex b1[2], ratio[2];   // first beam from the back face and the factor from beam to beam, for s and p
      for (int polarisation = 0; polarisation < 2; ++polarisation) {
        std::vector<Complex> eta(n), q(n);
        for (std::size_t j = 0; j < n; ++j) {
          eta[j] = polarisation == 0 ? media.m[j].eta_s : media.m[j].eta_p;
          q[j] = polarisation == 0 ? media.m[j].q_s : media.m[j].q_p;
        }
        const Amplitudes down = through(eta, q, media.thickness, wavelength_nm);
        std::vector<Complex> eta_up(eta.rbegin(), eta.rend()), q_up(q.rbegin(), q.rend());
        std::vector<double> thickness_up(media.thickness.rbegin(), media.thickness.rend());
        const Amplitudes up = through(eta_up, q_up, thickness_up, wavelength_nm);
        const Complex r_back = fresnel(eta.back(), polarisation == 0 ? back.eta_s : back.eta_p);
        b1[polarisation] = down.t * up.t * r_back;
        ratio[polarisation] = up.r * r_back;
      }
      const double A2 = pass * pass;   // down and up
      I_ss += std::norm(b1[0]) * A2 / (1.0 - std::norm(ratio[0]) * A2);
      I_pp += std::norm(b1[1]) * A2 / (1.0 - std::norm(ratio[1]) * A2);
      I_ps += b1[1] * std::conj(b1[0]) * A2 / (1.0 - ratio[1] * std::conj(ratio[0]) * A2);
    }
  }

  Observation o;
  o.Rs = I_ss;
  o.Rp = I_pp;
  const double total = I_ss + I_pp;
  o.N = (I_ss - I_pp) / total;
  o.C = 2.0 * I_ps.real() / total;
  o.S = 2.0 * I_ps.imag() / total;
  o.polarisation = std::sqrt(o.N * o.N + o.C * o.C + o.S * o.S);
  o.psi_deg = 0.5 * std::acos(std::clamp(o.N, -1.0, 1.0)) * 180.0 / M_PI;
  double delta = std::atan2(o.S, o.C) * 180.0 / M_PI;
  if (delta < 0.0) delta += 360.0;
  o.delta_deg = delta;
  return o;
}

PsiDelta psi_delta(const Reflection& reflection) {
  const Complex rho = reflection.rp / reflection.rs;
  PsiDelta result;
  result.psi_deg = std::atan(std::abs(rho)) * 180.0 / M_PI;
  double delta = std::arg(rho) * 180.0 / M_PI;
  if (delta < 0.0) delta += 360.0;
  result.delta_deg = delta;
  return result;
}

PsiDelta psi_delta(const Stack& stack, double wavelength_nm, double angle_deg) {
  // without a back face the coherent expression keeps Psi exact also near 0 and 90 degrees
  if (!(stack.substrate_thickness_nm > 0.0)) return psi_delta(reflect(stack, wavelength_nm, angle_deg));
  const Observation o = observe(stack, wavelength_nm, angle_deg);
  return {o.psi_deg, o.delta_deg};
}

}  // namespace ellidft::ellipsometry
