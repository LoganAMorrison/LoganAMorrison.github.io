---
title: Computing the Thermal Cross-Section for $4\eta'\to2\eta'$ in a Large-N SU(N) Gauge Theory
tags: quantum-field-theory cosmology numerics
---

## Introduction
Our goal of this notebook is to compute the thermal cross section for
$2\eta'\to4\eta'$ scattering. The number changing interactions for the
$\eta'$ are discribed by the Chiral Lagrangain. When expanding in powers of
$p^2/\Lambda^2$ and $N$ ($p$ being the momentum of the $\eta'$,
$\Lambda$ being the confinement scale and $N$ the number of colors),
the Lagrangian can be written as:
$$
\mathcal{L} = \dfrac{1}{2}(\partial_{\mu}\eta')^2 -\dfrac{1}{2}m_{\eta'}(\eta')^2 + \underbrace{\dfrac{\tilde{L}_{1}}{\left(f_{\eta'}\Lambda\right)^2}\left(\partial_{\mu}\eta'\right)^4}_{\mathcal{O}(p^4)}+\underbrace{\dfrac{\tilde{L}_{2}}{\left(f_{\eta'}\Lambda\right)^4}\left(\partial_{\mu}\eta'\right)^6}_{\mathcal{O}(p^6)} + \underbrace{\cdots}_{\mathcal{O}(p^8)}
$$
where we've underlined the order of the terms in the chiral expansion. When
considering $2\eta'\to4\eta'$, both the $\mathcal{O}(p^4)$ and
$\mathcal{O}(p^6)$ terms contribute. From $\mathcal{O}(p^4)$, we find
4-point interactions which can be chained together via $\eta'$ propagators
to achieve $2\eta'\to4\eta'$. From $\mathcal{O}(p^6)$, we find 6-point
contact interaction which directly contributes to $2\eta'\to4\eta'$.

Ultimately, we would like to compute $\langle\sigma(2\eta'\to4\eta') v\rangle$.
This is comuted from the following one-dimensional integrate over energies:
$$
\langle\sigma v\rangle = \dfrac{2\pi^2T_{\eta'}\int_{4m_{\eta'}^2}^{\infty}\sigma_{2\eta'\to4\eta'}(s)(s-4m_{\eta'}^2)\sqrt{s}K_{1}(\sqrt{s}/T_{\eta'})ds}{\left(4\pi m_{\eta'}^2 T K_{2}(m_{\eta'}/T_{\eta'})\right)^2}
$$
where $s$ is the squared CM energy and $T_{\eta'}$ is the $\eta'$ temperature.
In order to compute this thermal average, we will need to compute
$\sigma_{2\eta'\to4\eta'}(s)$: the zero temperature cross section.
This is given by:
$$
\sigma_{2\eta'\to4\eta'}(s) = \dfrac{1}{2\sqrt{s}\sqrt{s-4m^2}}\int\prod_{i=1}^{4}\dfrac{d^3p_{i}}{(2\pi)^2 2E_{i}}|\mathcal{M}|^2(2\pi)^4\delta^{4}\left(p_{i,1} + p_{i,2} - \sum_{i=1}^{4}p_{i})\right)
$$
with $p_{i,1}, p_{i,2}$ being the initial state $\eta'$ momenta, $p_{i}$ being
the final state $\eta'$ momenta, $E_{i}$ being the final-state $\eta'$
energies and $\mathcal{M}$ being the matrix element. The matrix element for
$2\eta'\to4\eta'$ is quite complicated. Thus, the phase-space integral is
difficult to obtain analytically. To compute the phase-space integral, we
will employ the RAMBO algorithm. To avoid computing the phase-space integral
for all values of $s$ and $m_{\eta'}$, we will rescale the phase-space
integral by the $\eta'$ mass, obtaining:
$$
\left(\dfrac{\tilde{L}^2_{1}m_{\eta'}^{8}}{f_{\eta'}^4\Lambda^4}\right)^2
f(\tilde{s}) =
\left(\dfrac{\tilde{L}^2_{1}m_{\eta'}^{8}}{f_{\eta'}^4\Lambda^4}\right)^2
\underbrace{\int\prod_{i=1}^{4}\dfrac{d^3q_{i}}{(2\pi)^2 2\epsilon_{i}}
|\tilde{\mathcal{M}}|^2(2\pi)^4\delta^{4}\left(q_{i,1} + q_{i,2} -
  \sum_{i=1}^{4}q_{i})\right)}_{\text{only function of } \tilde{s} }
$$
where $\tilde{s} = s/m_{\eta'}^2$, $q_{i} = p_{i} / m_{\eta'}$,
$\epsilon = E_{i} / m_{\eta'}$ and $\tilde{\mathcal{M}}$ is the matrix
element with all dimensions explicity factored out. In terms of
$f(\tilde{s})$, the cross section can be written as:
$$
\sigma_{2\eta'\to4\eta'}(s) = \dfrac{1}{2\sqrt{s}\sqrt{s-4m^2}}
\left(\dfrac{\tilde{L}^2_{1}m_{\eta'}^{8}}{f_{\eta'}^4\Lambda^4}\right)^2
f(\tilde{s})
$$
Our first task will be to compute $f(\tilde{s})$.

## Generate Cross section data for $2\eta'\to4\eta'$
To compute the scaleless $2\eta'\to4\eta'$ cross section ($f(\tilde{s})$
defined above), we will use `Rambo.jl`:

````julia
using Rambo
````





The matrix element for $2\eta'\to4\eta'$ is given by:

$$\begin{align}
  \mathcal{M} &= -\frac{48 L_1}{f_{\eta}^4 \Lambda^4} \bigg{[}p_{12} p_{34} p_{56}+p_{12} p_{35} p_{46} +
  p_{12} p_{36} p_{45}\\
  &\hspace{2.2cm}+p_{13} p_{24} p_{56}+p_{13} p_{25} p_{46}+
  p_{13} p_{26} p_{45}+p_{14} p_{23} p_{56}\notag\\
  &\hspace{2.2cm}+p_{14} p_{25} p_{36}+
  p_{14} p_{26} p_{35}+p_{15} (p_{23} p_{46}+p_{24} p_{36}+
  p_{26} p_{34})\notag\\
  &\hspace{2.2cm}+p_{16} (p_{23} p_{45}+p_{24} p_{35}+p_{25} p_{34})\bigg{]}\notag
\end{align}$$
where we defined $p_{ij} = p_{i}\cdot p_{j}$. In `julia`, this is (factoring
out the overall coefficient parameters, i.e. $L_{1}/(f_{\eta}\Lambda)^4$):

````julia
"""
Squared matrix element for 2η' → 4η' with all momenta rescaled by mη.

# Arguments
- `momenta::Array{FourMomentum, 1}`: four-momenta of final-state η's
"""
function scaled_msqrd(momenta)
  # Extract the center of mass energy
  Q::Float64 = sum(momenta).e
  # Compute the magnitude of the initial state eta' 3-momentum
  p::Float64 = sqrt(Q^2 / 4 - 1)
  # Chose the initial state eta's to be in the CM frame traveling
  # along z-axis
  q1 = FourMomentum(Q / 2, 0.0, 0.0, p)
  q2 = FourMomentum(Q / 2, 0.0, 0.0, -p)
  # Extract the final state eta' momenta
  q3, q4, q5, q6 = momenta

  return (528 *(
    scalar_product(q1,q4)*scalar_product(q2,q6)*scalar_product(q3,q5) +
    scalar_product(q1,q4)*scalar_product(q2,q5)*scalar_product(q3,q6) +
    scalar_product(q1,q3)*scalar_product(q2,q6)*scalar_product(q4,q5) +
    scalar_product(q1,q2)*scalar_product(q3,q6)*scalar_product(q4,q5) +
    scalar_product(q1,q6)*(scalar_product(q2,q5)*scalar_product(q3,q4) +
    scalar_product(q2,q4)*scalar_product(q3,q5) +
    scalar_product(q2,q3)*scalar_product(q4,q5)) +
    scalar_product(q1,q3)*scalar_product(q2,q5)*scalar_product(q4,q6) +
    scalar_product(q1,q2)*scalar_product(q3,q5)*scalar_product(q4,q6) +
    scalar_product(q1,q5)*(scalar_product(q2,q6)*scalar_product(q3,q4) +
    scalar_product(q2,q4)*scalar_product(q3,q6) +
    scalar_product(q2,q3)*scalar_product(q4,q6)) +
    scalar_product(q1,q4)*scalar_product(q2,q3)*scalar_product(q5,q6) +
    scalar_product(q1,q3)*scalar_product(q2,q4)*scalar_product(q5,q6) +
    scalar_product(q1,q2)*scalar_product(q3,q4)*scalar_product(q5,q6)))^2
end;
````





In order to make sure that our result using the scale-less momentum is correct,
we will also integrate the scale-full phase space so we can compare the results.
It turns out the the scale-full squared matrix element is just $m_{\eta}^6$
times the scaleless matrix element. Let's compute both the scaled and unscaled
phase space integrals. For the scaled phase-space integral, we will multiply
by $m_{\eta}^{14}$ in order to restore dimension.

````julia
import PyPlot; const plt = PyPlot # python plotting
using LaTeXStrings
# Pick a random value of η' mass
mη = 100.0
# Pick CM energies from 4*mη → 100*mη
scaled_cmes = 10 .^(range(log10(4.0 + 1e-5), stop=log10(100),length=100))
# The scaled masses are all unity
scaled_masses = [1.0, 1.0, 1.0, 1.0]
# number of events to generate for each phase space integral
nevents = 10000
# Compute the values of the scaled phase-space
scaled_ps = [integrate_phase_space(cme, scaled_masses; nevents=10000,
                                   msqrd=scaled_msqrd)[1]
             for cme in scaled_cmes]
# Compute the values of the unscaled phase-space
unscaled_ps = [integrate_phase_space(cme * mη, [mη,mη,mη,mη]; nevents=10000,
                                     msqrd=fm-> scaled_msqrd(fm))[1]
            for cme in scaled_cmes]

plt.figure(dpi=100)
plt.plot(scaled_cmes .* mη, scaled_ps .* mη^16, lw=3, alpha=0.6,
         label="scaled")
plt.plot(scaled_cmes .* mη, unscaled_ps, "k--",label="unscaled")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(L"$\sqrt{s}$", fontsize=16)
plt.ylabel("Phase Space", fontsize=16)
plt.legend()
plt.gcf()
````


![]({{site.baseurl}}/assets/img/2019-08-23-thermal-cross-section-2-to-4-eta-prime_3_1.png)


As we can see, the two agree after we restore an overall factor of
$m_{\eta}^{16}$ for the scaled phase space.

## Thermal Cross Section
Here we compute the thermal average of the cross section for
$2\eta'\to4\eta'$. We write the cross section as:
$$\begin{align}
    \sigma(s) &= \left(\dfrac{\tilde{L}^2_{1}m_{\eta'}^{8}}{f_{\eta'}^4\Lambda^4}\right)^2\dfrac{f(\sqrt{s}/m)}{2\sqrt{s}\sqrt{s-4m^2}}
\end{align}$$
where $s$ is the square of the CM energy and $m$ is the $\eta'$ mass. In the
last section, we computed $f(\sqrt{s}/m)$, which is given by:
$$\begin{align}
    f(\sqrt{s}/m) &= \int\prod_{i=1}^{4}\dfrac{d^3q_{i}}{(2\pi)^2 2\epsilon_{i}}|\tilde{\mathcal{M}}|^2(2\pi)^4\delta^{4}\left(q_{i,1} + q_{i,2} - \sum_{i=1}^{4}q_{i})\right)
\end{align}$$
here $q_{i}$ and $\epsilon_{i}$ are the momentum/energy of one of the
$\eta'$s scaled by the mass: $q_{i} = p_{i}/m$ and $\epsilon_{i} = E_{i}/m$.
We found that $n = 16$.

Given this parameterization of the cross section, we wish to compute the
thermal average, which is given by:
$$\begin{align}
    \langle\sigma v\rangle &= \dfrac{2\pi^2T\int_{4m^2}^{\infty}\sigma(s)(s-4m^2)\sqrt{s}K_{1}(\sqrt{s}/T)ds}{\left(4\pi m^2 T K_{2}(m/T)\right)^2}
\end{align}$$
To evaluate this, we will first make a change of variables from
$(s, m) \to (z, x)$. We will let:
$$\begin{align}
    s &= z^2 T^2, & m = x T.
\end{align}$$
Making the change of variables, our integral becomes:
$$\begin{align}
    \langle\sigma v\rangle &= \dfrac{1}{8x^2K_{2}^2(x)}\left(\dfrac{\tilde{L}^2_{1}m_{\eta'}^{7}}{f_{\eta'}^4\Lambda^4}\right)^2\int_{4x}^{\infty} z \ f\left(\frac{z}{x}\right)\sqrt{z^2-4x^2}K_{1}(z)dz
\end{align}$$
We can see that the integral is just some function of $x$. Let's give it a
name: $f_{T}(x)$:
$$\begin{align}
    f_{T}(x) &\equiv \dfrac{1}{8x^2K_{2}^2(x)}\int_{4x}^{\infty} z \ f\left(\frac{z}{x}\right)\sqrt{z^2-4x^2}K_{1}(z)dz
\end{align}$$
Then, the terms of $f_{T}(x)$, our thermally averaged cross section is:
$$\begin{align}
    \langle\sigma v\rangle &= \left(\dfrac{\tilde{L}^2_{1}m_{\eta'}^{7}}{f_{\eta'}^4\Lambda^4}\right)^2f_{T}(x)
\end{align}$$
Thus, our goal will be to compute $f_{T}(x)$.

First, we make an interpolating function for the zero-temperature cross
section
````julia
using Dierckx
# Generate new, dense data set for the scaled phase-space
cmes = 10 .^(range(log10(4.0 + 1e-5), stop=log10(100),length=500))
masses = [1.0, 1.0, 1.0, 1.0]
nevents = 10000
ps = [integrate_phase_space(cme, masses; nevents=10000, msqrd=scaled_msqrd)[1]
      for cme in cmes]
# Generate a linear (k = 1) spline of the data
spl = Spline1D(cmes, ps, k=1, bc="nearest", s=0.0);
````





Next, we make a fit to the large energy behavior

````julia
using LsqFit
# Model the large energy behavior as ps(cme) = cme^p₁ * 10^p₂
@. model(x, p) = x^(p[1]) * 10^(p[2])
# Take only the last twenty elements of the data to capture large energies
xdata = cmes[end-20:end]
ydata = ps[end-20:end]
fit = curve_fit(model, xdata, ydata, [16.0, -6.0]);
m, b = coef(fit)
````


````
2-element Array{Float64,1}:
 16.030261894251517
 -4.9700402079867185
````





Now we make a new function which uses the large energy behavior:

````julia
"""
Computes the cross section for 2η → 4η at zero temperature.

# Arguments
- `x::Float64`: η' mass scaled by temperature.
"""
function cs_interp(cme::Float64)
  if cme < cmes[1]
    return 0.0
  elseif cme > cmes[end]
    return cme^m * 10^b
  else
      return spl(cme)
  end
end;
````





Let's plot and see that the interpolating function is working properly:

````julia
cmes_ext = 10 .^(range(log10(4.0 + 1e-5), stop=log10(150),length=500))
plt.figure(dpi=100)
plt.plot(cmes, ps, lw=3, alpha=0.6, label="data")
plt.plot(cmes_ext, cs_interp.(cmes_ext), "k--",label="interp")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(L"$\sqrt{s}$", fontsize=16)
plt.ylabel("Phase Space", fontsize=16)
plt.legend()
plt.gcf()
````


![]({{site.baseurl}}/assets/img/2019-08-23-thermal-cross-section-2-to-4-eta-prime_7_1.png)



We are now ready to begin computing the thermal average of the cross section.
We will need a function for the integrand of the thermal average:

````julia
using SpecialFunctions
using QuadGK

"""
Thermal kernal for computing thermally averaged cross section.

# Arguments
- `z::Float64`: center of mass energy scaled by the temperature.
- `x::Float64`: η' mass scaled by temperature.
"""
function thermal_kernal(z::Float64, x::Float64)
  return z * sqrt(z^2 - 4x^2) * besselk(1, z)
end

"""
Computes the thermal cross section with couplings factored out.

# Arguments
- `x::Float64`: η' mass scaled by temperature.
"""
function thermal_cs(x)
  function integrand(z)
    return thermal_kernal(z, x) * cs_interp(z / x)
  end
  return quadgk(integrand, 4x, Inf)[1] / (8x^2 * besselk(2, x)^2)
end;
````





Now we compute and plot the results

````julia
xs = 10 .^(range(-3, stop=log10(200), length=500))
tcs = [thermal_cs(x) for x in xs]

plt.figure(dpi=100)
plt.plot(xs, tcs, "k--")
plt.yscale("log")
plt.xscale("log")
plt.ylabel(L"$\langle\sigma v\rangle\times\left(\frac{f_{\eta}^4\Lambda^4}{m_{\eta}^7\tilde{L}_{1}^2}\right)^2$", fontsize=16)
plt.xlabel(L"$x$", fontsize=16)
plt.gcf()
````


![]({{site.baseurl}}/assets/img/2019-08-23-thermal-cross-section-2-to-4-eta-prime_9_1.png)

````julia
plt.close_figs()
````
