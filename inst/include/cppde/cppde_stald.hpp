/*
 Stability limit detection (STALD) for the BDF family, a port of cvBDFStab and
 cvSLdet in CVODE. From the scaled derivatives of degree q-1, q and q+1 over
 the last five steps it estimates the magnitude rr of the dominant
 characteristic root of the step-to-step error propagation; rr near 1 means
 the step size sits at the edge of the stability region of order q.

 Copyright (C) 2026 Simon Beyer

 Derived from SUNDIALS/CVODE (cvBDFStab and cvSLdet in cvode.c),
 Copyright (c) 2002-2024 Lawrence Livermore National Security and
 Southern Methodist University, distributed under the BSD-3-Clause
 license.  See inst/COPYRIGHTS for the full license text.
 */

#ifndef CPPDE_STALD_HPP
#define CPPDE_STALD_HPP

#include <algorithm>
#include <array>
#include <cmath>

namespace cppde {

// The squared scaled derivative norms, ss[i][k]: i = 1 is the latest step,
// i = 5 the oldest, k = 1, 2, 3 the degrees q-1, q, q+1. Index 0 is unused,
// so the code reads as CVODE's.
struct stald_data {
  std::array<std::array<double, 4>, 6> ss{};

  void clear() { for (auto& r : ss) r.fill(0.0); }

  void push(double sqm2, double sqm1, double sq)
  {
    for (int i = 5; i >= 2; --i) ss[i] = ss[i - 1];
    ss[1][1] = sqm2 * sqm2;
    ss[1][2] = sqm1 * sqm1;
    ss[1][3] = sq * sq;
  }
};

// cvSLdet: 1 to 3 found a stable root, 4 to 6 a stability limit violation
// (rr above 0.98), negative when the data admit no estimate.
inline int stald_detect(const stald_data& d, int q)
{
  constexpr double TINY = 1e-10, rrcut = 0.98, vrrtol = 1e-4, vrrt2 = 5e-4,
                   sqtol = 1e-3, rrtol = 1e-2;
  const auto& ss = d.ss;

  double rat[5][4], rav[4], qkr[4], sigsq[4], smax[4], ssmax[4];
  double drr[4], rrc[4], sqmx[4], qjk[4][4], vrat[4], qc[6][4], qco[6][4];
  double rr = 0.0, sqmin = 0.0;
  int kflag = 0, kmin = 0;

  for (int k = 1; k <= 3; ++k) {
    double smink = ss[1][k], smaxk = 0.0;
    for (int i = 1; i <= 5; ++i) {
      smink = std::min(smink, ss[i][k]);
      smaxk = std::max(smaxk, ss[i][k]);
    }
    if (smink < TINY * smaxk) return -1;
    smax[k]  = smaxk;
    ssmax[k] = smaxk * smaxk;

    double sumrat = 0.0, sumrsq = 0.0;
    for (int i = 1; i <= 4; ++i) {
      rat[i][k] = ss[i][k] / ss[i + 1][k];
      sumrat += rat[i][k];
      sumrsq += rat[i][k] * rat[i][k];
    }
    rav[k]  = 0.25 * sumrat;
    vrat[k] = std::abs(0.25 * sumrsq - rav[k] * rav[k]);

    qc[5][k] = ss[1][k] * ss[3][k] - ss[2][k] * ss[2][k];
    qc[4][k] = ss[2][k] * ss[3][k] - ss[1][k] * ss[4][k];
    qc[3][k] = 0.0;
    qc[2][k] = ss[2][k] * ss[5][k] - ss[3][k] * ss[4][k];
    qc[1][k] = ss[4][k] * ss[4][k] - ss[3][k] * ss[5][k];
    for (int i = 1; i <= 5; ++i) qco[i][k] = qc[i][k];
  }

  const double vmin = std::min({vrat[1], vrat[2], vrat[3]});
  const double vmax = std::max({vrat[1], vrat[2], vrat[3]});

  if (vmin < vrrtol * vrrtol) {
    // Normal or nearly normal matrix: the three ratio sequences are constant.
    if (vmax > vrrt2 * vrrt2) return -2;
    rr = (rav[1] + rav[2] + rav[3]) / 3.0;
    double drrmax = 0.0;
    for (int k = 1; k <= 3; ++k) drrmax = std::max(drrmax, std::abs(rav[k] - rr));
    // Ratios that disagree admit no estimate.
    if (drrmax > vrrt2) return -3;
    kflag = 1;
  } else {
    // The common root of the three quartics, by elimination.
    if (std::abs(qco[1][1]) < TINY * ssmax[1]) return -4;
    double tem = qco[1][2] / qco[1][1];
    for (int i = 2; i <= 5; ++i) qco[i][2] -= tem * qco[i][1];
    qco[1][2] = 0.0;
    tem = qco[1][3] / qco[1][1];
    for (int i = 2; i <= 5; ++i) qco[i][3] -= tem * qco[i][1];
    qco[1][3] = 0.0;

    if (std::abs(qco[2][2]) < TINY * ssmax[2]) return -4;
    tem = qco[2][3] / qco[2][2];
    for (int i = 3; i <= 5; ++i) qco[i][3] -= tem * qco[i][2];

    if (std::abs(qco[4][3]) < TINY * ssmax[3]) return -4;
    rr = -qco[5][3] / qco[4][3];
    if (rr < TINY || rr > 100.0) return -5;

    for (int k = 1; k <= 3; ++k)
      qkr[k] = qc[5][k] + rr * (qc[4][k] + rr * rr * (qc[2][k] + rr * qc[1][k]));

    double sqmax = 0.0;
    for (int k = 1; k <= 3; ++k) sqmax = std::max(sqmax, std::abs(qkr[k]) / ssmax[k]);

    if (sqmax < sqtol) {
      kflag = 2;
    } else {
      // Newton corrections of rr, keeping the candidate that fits all three.
      for (int it = 1; it <= 3; ++it) {
        for (int k = 1; k <= 3; ++k) {
          const double qp = qc[4][k] + rr * rr * (3.0 * qc[2][k] + rr * 4.0 * qc[1][k]);
          drr[k] = (std::abs(qp) > TINY * ssmax[k]) ? -qkr[k] / qp : 0.0;
          rrc[k] = rr + drr[k];
        }
        for (int k = 1; k <= 3; ++k) {
          const double s = rrc[k];
          double sqmaxk = 0.0;
          for (int j = 1; j <= 3; ++j) {
            qjk[j][k] = qc[5][j] + s * (qc[4][j] + s * s * (qc[2][j] + s * qc[1][j]));
            sqmaxk = std::max(sqmaxk, std::abs(qjk[j][k]) / ssmax[j]);
          }
          sqmx[k] = sqmaxk;
        }
        sqmin = sqmx[1] + 1.0;
        for (int k = 1; k <= 3; ++k)
          if (sqmx[k] < sqmin) { kmin = k; sqmin = sqmx[k]; }
        rr = rrc[kmin];
        if (sqmin < sqtol) { kflag = 3; break; }
        for (int j = 1; j <= 3; ++j) qkr[j] = qjk[j][kmin];
      }
      if (sqmin > sqtol) return -6;
    }
  }

  // Given rr, the variances sigsq[k], and rr once more from them.
  for (int k = 1; k <= 3; ++k) {
    const double rsa = ss[1][k], rsb = ss[2][k] * rr,
                 rsc = ss[3][k] * rr * rr, rsd = ss[4][k] * rr * rr * rr;
    const double rd1a = rsa - rsb, rd1b = rsb - rsc, rd1c = rsc - rsd;
    const double rd2a = rd1a - rd1b, rd2b = rd1b - rd1c, rd3a = rd2a - rd2b;
    if (std::abs(rd1b) < TINY * smax[k]) return -7;
    const double cest1 = -rd3a / rd1b;
    if (cest1 < TINY || cest1 > 4.0) return -7;
    sigsq[k] = ss[3][k] + (rd2b / cest1) / (rr * rr);
  }
  if (sigsq[2] < TINY) return -8;

  const double ratp = sigsq[3] / sigsq[2], ratm = sigsq[1] / sigsq[2];
  const double qfac1 = 0.25 * (q * q - 1.0), qfac2 = 2.0 / (q - 1.0);
  const double bb  = ratp * ratm - 1.0 - qfac1 * ratp;
  const double tem = 1.0 - qfac2 * bb;
  if (std::abs(tem) < TINY) return -8;
  if (std::abs(1.0 / tem - rr) > rrtol) return -9;

  return rr > rrcut ? kflag + 3 : kflag;
}

} // namespace cppde

#endif // CPPDE_STALD_HPP
