---
title: Computing the Thermal Cross-Section for $\eta'\eta'\to\Delta\bar{\Delta}$ in a Large-N SU(N) Gauge Theory
tags: quantum-field-theory cosmology numerics
---

## Introduction
In this notebook, we wish to compute the thermal average of the
annihilation cross section for $\eta'\eta'\leftrightarrow\Delta\Delta$ in
the large $N$ limit. We will assume the amplitude describing this interaction
is
$$\begin{align}
  \mathcal{M}(\bar{\Delta}\Delta\leftrightarrow\eta'\eta') &=
  \frac{e^{-cN}}{\Lambda}\bar{v}_{\Delta}(p_{2})u_{\Delta}(p_1)
\end{align}$$
where $c$ is a dimensionless parameter describing the probability of
combining colors in the right way to form a baryon from a meson or vis-versa.
The factor of $\Lambda$ is here to correct for the mass dimensions and
$u_{\Delta}, v_{\Delta}$ are spinors. Squaring our the matrix element and
summing over spins, we find:
$$\begin{align}
  |\mathcal{M}(\bar{\Delta}\Delta\leftrightarrow\eta'\eta')|^2 &=
  \frac{e^{-2cN}}{\Lambda^2}\left(p_{1}\cdot p_{2}-m_{\Delta}^2\right)\\
\end{align}$$
Writing $p_{1}\cdot p_{2} = s/2 - m_{\Delta}^2$, with $s$ being the square
or the CM energy, we find:
$$\begin{align}
  |\mathcal{M}(\bar{\Delta}\Delta\leftrightarrow\eta'\eta')|^2 &=
  \frac{e^{-2cN}}{2\Lambda^2}\left(s-4m_{\Delta}^2\right)\\
\end{align}$$
From now on, we will focus on $\eta'\eta'\to\bar{\Delta}\Delta$, since
we can get the thermal cross section for $\bar{\Delta}\Delta\to\eta'\eta'$
from detailed balance. If we integrate over phase-space, we obtain the following
cross-sections:
$$\begin{align}
  \sigma_{\eta'\eta'\to\bar{\Delta}\Delta}(s) &=
  \frac{e^{-2cN}}{32\pi\Lambda^2}
  \frac{\sqrt{s-4m_{\eta}^2}\sqrt{s-4m_{\Delta}^2}}{s}
\end{align}$$
To compute the thermal average of this cross section, we need to evaluate:
$$\begin{align}
\langle\sigma v\rangle = \dfrac{2\pi^2T_{\eta'}\int_{4m_{\eta'}^2}^{\infty}
\sigma_{\eta'\eta'\to\bar{\Delta}\Delta}(s)(s-4m_{\eta'}^2)
\sqrt{s}K_{1}(\sqrt{s}/T_{\eta'})ds}{\left(4\pi m_{\eta'}^2 T K_{2}(m_{\eta'}/T_{\eta'})\right)^2}
\end{align}$$
Inserting in the expression for the cross section, we find that this is:
$$\begin{align}
\langle\sigma v\rangle = \frac{e^{-2cN}}{256\pi\Lambda^2sT^3x_{\eta}^4K^2_{2}(x_{\eta})}
\int_{4m_{\Delta}^2}^{\infty}
s^{-1/2}\sqrt{s-4m_{\Delta}^2}(s-4m_{\eta'}^2)^{3/2}
K_{1}(\sqrt{s}/T_{\eta'})ds
\end{align}$$
Notice that we changed the bounds of integration due to the mass threshold
and the fact that $m_{\Delta} > m_{\eta}$. It will be convinent to realize that
$m_{\eta} \sim N^{-3/2}m_{\Delta}$. If we use this result, we find that
$$\begin{align}
(s-4m_{\eta'}^2)^{3/2} = s^{3/2} -
6\sqrt{s}\left(\frac{\mu_{\eta}m_{\Delta}}{\mu_{\Delta}N^{3/2}}\right)^2 +
\mathcal{O}(N^{-4})
\end{align}$$
where we have identified $m_{\Delta} = \mu_{\Delta}\Lambda N$ and
$m_{\eta} = \mu_{\eta}\Lambda/\sqrt{N}$. If we make this expansion in the
integrand, preform the integration and expand in powers of $N$, we find that:
$$\begin{align}
\langle\sigma v\rangle &= \frac{m_{\Delta}^3e^{-2m_{\Delta}/T-2cN}}{\Lambda^2}
\left(\frac{1}{128\Lambda^2 T} +
\frac{1}{512T^2m_{\Delta}}\left(9T^2+\frac{2\Lambda^2m_{\Delta}}{N}\right)+
\mathcal{O}(N)\right)
\end{align}$$

Let's compare our result with the exact result for a couple of benchmark point.
We will choose $c=\Lambda = \mu_{\Delta} = \mu_{\eta} = 1$ and invesigate the
results for various $N$.

````julia
using QuadGK
using SpecialFunctions

mutable struct DarkSUNParams
  N::Float64
  Λ::Float64
  c::Float64
  μη::Float64
  μΔ::Float64
end

function σ_ηη_ΔΔ(s::Float64, params::DarkSUNParams)
  Λ::Float64 = params.Λ
  N::Float64 = params.N
  μη::Float64 = params.μη
  μΔ::Float64 = params.μΔ
  c::Float64 = params.c

  mη::Float64 = μη * Λ / sqrt(N)
  mΔ::Float64 = μΔ * N * Λ

  exp(-2N * c) * sqrt(s-4mη^2) * sqrt(s-4mΔ^2) / (32π * s * Λ^2)
end

function exact_thermal_cs(T::Float64, params::DarkSUNParams)
  Λ::Float64 = params.Λ
  N::Float64 = params.N
  μη::Float64 = params.μη
  μΔ::Float64 = params.μΔ
  mη::Float64 = μη * Λ / sqrt(N)
  mΔ::Float64 = μΔ * N * Λ

  pf::Float64 = (2π^2 * T) / (4π * mη^2 * T * besselk(2, mη/T))^2

  function integrand(s::Float64)
    return σ_ηη_ΔΔ(s, params) * (s-4mη^2) * sqrt(s) * besselk(1, sqrt(s)/T)
  end

  return quadgk(integrand, 4mΔ^2, Inf)[1] * pf
end

function approx_thermal_cs(T::Float64, params::DarkSUNParams)
  Λ::Float64 = params.Λ
  N::Float64 = params.N
  μη::Float64 = params.μη
  μΔ::Float64 = params.μΔ
  c::Float64 = params.c
  mΔ::Float64 = Λ * N * μΔ

  mΔ^3 / (T^3 * Λ^2) * exp(-2mΔ/T- 2c * N) * (1/128 +
    1/(512T^2 * mΔ) * (9T^3 + 2mΔ * Λ^2 * μη^2 / N))
end
````


````
approx_thermal_cs (generic function with 1 method)
````





Plot the results:
````julia
import PyPlot; const plt = PyPlot # python plotting
using LaTeXStrings

params = DarkSUNParams(5.0, 1.0, 1.0, 1.0, 1.0)

Ts = 10 .^(range(log10(params.Λ / 100), stop=log10(10 * params.Λ), length=100))
exact_tcs0 = [exact_thermal_cs(T, params) for T in Ts]
approx_tcs0 = [approx_thermal_cs(T, params) for T in Ts]

params.N = 10.0
exact_tcs1 = [exact_thermal_cs(T, params) for T in Ts]
approx_tcs1 = [approx_thermal_cs(T, params) for T in Ts]

params.N = 20.0
exact_tcs2 = [exact_thermal_cs(T, params) for T in Ts]
approx_tcs2 = [approx_thermal_cs(T, params) for T in Ts]

params.N = 30.0
exact_tcs3 = [exact_thermal_cs(T, params) for T in Ts]
approx_tcs3 = [approx_thermal_cs(T, params) for T in Ts]

plt.figure(dpi=100)
plt.title("lam = $(params.Λ), c = $(params.c), mud = $(params.μΔ), mue = $(params.μη)")
plt.plot(Ts, exact_tcs0, lw=3, alpha=0.6, label="exact N = 5")
plt.plot(Ts, exact_tcs1, lw=3, alpha=0.6, label="exact N = 10")
plt.plot(Ts, exact_tcs2, lw=3, alpha=0.6, label="exact N = 20")
plt.plot(Ts, exact_tcs3, lw=3, alpha=0.6, label="exact N = 30")
plt.plot(Ts, approx_tcs0, "k--")
plt.plot(Ts, approx_tcs1, "k--")
plt.plot(Ts, approx_tcs2, "k--")
plt.plot(Ts, approx_tcs3, "k--")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(L"$T$", fontsize=16)
plt.ylabel(L"\langle\sigma v\rangle", fontsize=16)
plt.legend()
plt.gcf()
````


![]({{site.baseurl}}/assets/img/2019-08-22-thermal-cross-section-eta-prime-delta_2_1.png)

````julia
import PyPlot; const plt = PyPlot # python plotting
using LaTeXStrings

params = DarkSUNParams(5.0, 1e-3, 1.0, 1.0, 1.0)

Ts = 10 .^(range(log10(params.Λ / 100), stop=log10(10 * params.Λ), length=100))
exact_tcs0 = [exact_thermal_cs(T, params) for T in Ts]
approx_tcs0 = [approx_thermal_cs(T, params) for T in Ts]

params.N = 10.0
exact_tcs1 = [exact_thermal_cs(T, params) for T in Ts]
approx_tcs1 = [approx_thermal_cs(T, params) for T in Ts]

params.N = 20.0
exact_tcs2 = [exact_thermal_cs(T, params) for T in Ts]
approx_tcs2 = [approx_thermal_cs(T, params) for T in Ts]

params.N = 30.0
exact_tcs3 = [exact_thermal_cs(T, params) for T in Ts]
approx_tcs3 = [approx_thermal_cs(T, params) for T in Ts]

plt.figure(dpi=100)
plt.title("lam = $(params.Λ), c = $(params.c), mud = $(params.μΔ), mue = $(params.μη)")
plt.plot(Ts, exact_tcs0, lw=3, alpha=0.6, label="exact N = 5")
plt.plot(Ts, exact_tcs1, lw=3, alpha=0.6, label="exact N = 10")
plt.plot(Ts, exact_tcs2, lw=3, alpha=0.6, label="exact N = 20")
plt.plot(Ts, exact_tcs3, lw=3, alpha=0.6, label="exact N = 30")
plt.plot(Ts, approx_tcs0, "k--")
plt.plot(Ts, approx_tcs1, "k--")
plt.plot(Ts, approx_tcs2, "k--")
plt.plot(Ts, approx_tcs3, "k--")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(L"$T$", fontsize=16)
plt.ylabel(L"\langle\sigma v\rangle", fontsize=16)
plt.legend()
plt.gcf()
````


![]({{site.baseurl}}/assets/img/2019-08-22-thermal-cross-section-eta-prime-delta_3_1.png)

````julia
import PyPlot; const plt = PyPlot # python plotting
using LaTeXStrings

params = DarkSUNParams(5.0, 1.0, 1.0, 5.0, 1.0)

Ts = 10 .^(range(log10(params.Λ / 100), stop=log10(10 * params.Λ), length=100))
exact_tcs0 = [exact_thermal_cs(T, params) for T in Ts]
approx_tcs0 = [approx_thermal_cs(T, params) for T in Ts]

params.N = 10.0
exact_tcs1 = [exact_thermal_cs(T, params) for T in Ts]
approx_tcs1 = [approx_thermal_cs(T, params) for T in Ts]

params.N = 20.0
exact_tcs2 = [exact_thermal_cs(T, params) for T in Ts]
approx_tcs2 = [approx_thermal_cs(T, params) for T in Ts]

params.N = 30.0
exact_tcs3 = [exact_thermal_cs(T, params) for T in Ts]
approx_tcs3 = [approx_thermal_cs(T, params) for T in Ts]

plt.figure(dpi=100)
plt.title("lam = $(params.Λ), c = $(params.c), mud = $(params.μΔ), mue = $(params.μη)")
plt.plot(Ts, exact_tcs0, lw=3, alpha=0.6, label="exact N = 5")
plt.plot(Ts, exact_tcs1, lw=3, alpha=0.6, label="exact N = 10")
plt.plot(Ts, exact_tcs2, lw=3, alpha=0.6, label="exact N = 20")
plt.plot(Ts, exact_tcs3, lw=3, alpha=0.6, label="exact N = 30")
plt.plot(Ts, approx_tcs0, "k--")
plt.plot(Ts, approx_tcs1, "k--")
plt.plot(Ts, approx_tcs2, "k--")
plt.plot(Ts, approx_tcs3, "k--")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(L"$T$", fontsize=16)
plt.ylabel(L"\langle\sigma v\rangle", fontsize=16)
plt.legend()
plt.gcf()
````


![]({{site.baseurl}}/assets/img/2019-08-22-thermal-cross-section-eta-prime-delta_4_1.png)

````julia
plt.close_figs()
````
