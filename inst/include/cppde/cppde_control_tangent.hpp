/*
 The control tangent's system: the linearised flow z' = J z, with J held in
 double as the -J a generated Jacobian writes, and the direction the tangent
 starts from. See multistepper::set_control_tangent().

 Copyright (C) 2026 Simon Beyer
 */

#ifndef CPPDE_CONTROL_TANGENT_HPP
#define CPPDE_CONTROL_TANGENT_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <cppde/cppde_ad_traits.hpp>
#include <cppde/cppde_types.hpp>

namespace cppde {
namespace tangent_detail {

// out = J z from the stored -J, dense column-major or CSC.
inline void apply_jacobian(const dense_matrix<double>& negJ,
                           const std::vector<double>& z,
                           std::vector<double>& out)
{
  const std::size_t n = z.size();
  out.assign(n, 0.0);
  const double* a = negJ.data.data();
  for (std::size_t j = 0; j < n; ++j, a += n) {
    const double zj = z[j];
    if (zj == 0.0) continue;
    for (std::size_t i = 0; i < n; ++i) out[i] -= a[i] * zj;
  }
}

inline void apply_jacobian(const csc_matrix<double>& negJ,
                           const std::vector<double>& z,
                           std::vector<double>& out)
{
  out.assign(z.size(), 0.0);
  for (int j = 0; j < negJ.n; ++j) {
    const double zj = z[static_cast<std::size_t>(j)];
    if (zj == 0.0) continue;
    for (int k = negJ.Ap[j]; k < negJ.Ap[j + 1]; ++k)
      out[static_cast<std::size_t>(negJ.Ai[k])] -= negJ.Ax[k] * zj;
  }
}

// The generated Jacobian at (x, t), written as -J into m.
template<class JacFunc>
inline void linearise(JacFunc& jf, const std::vector<double>& x, double t,
                      dense_matrix<double>& m, std::vector<double>& dfdt)
{
  const int n = static_cast<int>(x.size());
  if (m.rows() != n) m.resize(n, n);
  dfdt.resize(x.size());
  jf(x, m, t, dfdt);
}

template<class JacFunc>
inline void linearise(JacFunc& jf, const std::vector<double>& x, double t,
                      csc_matrix<double>& m, std::vector<double>& dfdt)
{
  dfdt.resize(x.size());
  jf(x, m, t, dfdt);
}

// Whether a generated Jacobian can be rebuilt on other parameters: its
// `params` and `F` members and the constructor taking both.
template<class JacFunc, class V, class = void>
struct rebuildable_jacobian : std::false_type {};

template<class JacFunc, class V>
struct rebuildable_jacobian<JacFunc, V, std::void_t<
    decltype(std::declval<JacFunc&>().params),
    decltype(JacFunc(std::declval<const std::vector<V>&>(),
                     std::declval<JacFunc&>().F))>> : std::true_type {};

// -J at a state, in double, from the model's own Jacobian. A model written in
// an AD scalar is evaluated on the values of the state and of its parameters,
// so the result is the Jacobian of the value run whatever width and content
// the tangents have.
template<class V, bool Sparse>
class linearisation {
public:
  using matrix = std::conditional_t<Sparse, csc_matrix<double>,
                                    dense_matrix<double>>;

  template<class JacFunc>
  void at(JacFunc& jf, const std::vector<V>& x, double t)
  {
    if constexpr (std::is_same_v<V, double>) {
      linearise(jf, x, t, m_J, m_dfdt);
    } else if constexpr (rebuildable_jacobian<JacFunc, V>::value) {
      m_p.resize(jf.params.size());
      for (std::size_t i = 0; i < m_p.size(); ++i)
        m_p[i] = V(ad_traits::scalar_value(jf.params[i]));
      JacFunc jv(m_p, jf.F);
      eval(jv, x, t);
    } else {
      eval(jf, x, t);
    }
  }

  const matrix& negJ() const { return m_J; }

private:
  template<class JacFunc>
  void eval(JacFunc& jf, const std::vector<V>& x, double t)
  {
    const std::size_t n = x.size();
    m_x.resize(n);
    for (std::size_t i = 0; i < n; ++i)
      m_x[i] = V(ad_traits::scalar_value(x[i]));
    m_dfdt_v.resize(n);
    if constexpr (Sparse) {
      jf(m_x, m_Jv, V(t), m_dfdt_v);
      if (!m_J.pattern_built) {
        m_J.n = m_Jv.n; m_J.nnz = m_Jv.nnz;
        m_J.Ap = m_Jv.Ap; m_J.Ai = m_Jv.Ai;
        m_J.Ax.assign(m_Jv.Ax.size(), 0.0);
        m_J.pattern_built = true;
      }
      for (std::size_t k = 0; k < m_Jv.Ax.size(); ++k)
        m_J.Ax[k] = ad_traits::scalar_value(m_Jv.Ax[k]);
    } else {
      const int ni = static_cast<int>(n);
      if (m_Jv.rows() != ni) m_Jv.resize(ni, ni);
      if (m_J.rows() != ni) m_J.resize(ni, ni);
      jf(m_x, m_Jv, V(t), m_dfdt_v);
      for (std::size_t k = 0; k < m_Jv.data.size(); ++k)
        m_J.data[k] = ad_traits::scalar_value(m_Jv.data[k]);
    }
  }

  matrix m_J;
  std::vector<double> m_dfdt;
  std::conditional_t<Sparse, csc_matrix<V>, dense_matrix<V>> m_Jv;
  std::vector<V> m_x, m_p, m_dfdt_v;
};

template<class Mat>
struct linear_rhs {
  const Mat* negJ;
  template<class T>
  void operator()(const std::vector<double>& z, std::vector<double>& dz,
                  const T& /*t*/) const
  { apply_jacobian(*negJ, z, dz); }
};

// Hands the stored -J to a stepper's own factorisation. Autonomous, so df/dt
// is zero.
template<class Mat>
struct linear_jac {
  const Mat* negJ;
  template<class T>
  void operator()(const std::vector<double>& /*z*/, Mat& W, const T& /*t*/,
                  std::vector<double>& dfdt) const
  {
    if constexpr (std::is_same_v<Mat, csc_matrix<double>>) {
      if (!W.pattern_built) W = *negJ;
      else std::copy(negJ->Ax.begin(), negJ->Ax.end(), W.Ax.begin());
    } else {
      W.data = negJ->data;
      W.n_rows = negJ->n_rows;
      W.n_cols = negJ->n_cols;
    }
    std::fill(dfdt.begin(), dfdt.end(), 0.0);
  }
};

// The direction a control tangent starts from: signs from a fixed hash of the
// index, so no mode of a structured model is missed by symmetry and a solve
// stays reproducible.
inline double start_direction(std::size_t i)
{
  std::uint32_t h = static_cast<std::uint32_t>(i) + 0x9e3779b9u;
  h ^= h >> 16; h *= 0x85ebca6bu;
  h ^= h >> 13; h *= 0xc2b2ae35u;
  h ^= h >> 16;
  return (h & 1u) ? 1.0 : -1.0;
}

} // namespace tangent_detail
} // namespace cppde

#endif // CPPDE_CONTROL_TANGENT_HPP
