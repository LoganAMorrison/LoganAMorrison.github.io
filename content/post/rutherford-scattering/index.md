---
title: Rutherford Scattering
subtitle: 
math: true

# Summary for listings and search engines
summary: 

# Link this post with a project
projects: []

# Date published
date: "2018-2-25"

# Date updated
lastmod: "2018-2-25"

# Is this an unpublished draft?
draft: false

# Show this page in the Featured widget?
featured: false

authors:
- Logan A. Morrison

tags:
- quantum-field-theory 

categories:
- Teaching

---

## Classical Calculation
$\require{cancel}$
$\newcommand{\bra}[1]{\left< #1 \right|}$
$\newcommand{\ket}[1]{\left| #1 \right>}$
$\newcommand{\bk}[2]{\left< #1 \middle| #2 \right>}$
$\newcommand{\bke}[3]{\left< #1 \middle| #2 \middle| #3 \right>}$
$\newcommand{\Tr}{\mathrm{Tr}}$

We consider a charged particle of mass $m$ scattering off a heavy
nucleus. Let the impact parameter be $b$ and the energy of the impinging
particle be $E$. The classical Lagrangian for a particle with charge $e$
scattering off a heavy nucleus of charge $eZ$ is:
$$\begin{align}
    \mathcal{L}= T - V = \dfrac{1}{2}m \left(\dfrac{d\vec{r}}{dt}\right)^2 - V(r)
\end{align}$$
with
$$\begin{align}
    V(r) = \dfrac{e^2Z}{4\pi r}
\end{align}$$
In the plane of scattering, we need only to consider two variables - $\|\vec{r}\| = r$ and $\theta$. We can write the derivative of $\vec{r}$ as:
$$\begin{align}
    \dfrac{d\vec{r}}{dt} = \dfrac{dr}{dt}\hat{\vec{r}} + r \dfrac{d\theta}{dt}\hat{\vec{\theta}}
\end{align}$$
Our Lagrangian then reads:
$$\begin{align}
    \mathcal{L}= \dfrac{1}{2}m \left(\left(\dfrac{dr}{dt}\right)^2 + r^2 \left(\dfrac{d\theta}{dt}\right)^2\right) - \dfrac{e^2Z}{4\pi r}
\end{align}$$
The Euler-Lagrange equation of $\theta$ reads:
$$\begin{align}
    0 = \dfrac{\partial\mathcal{L}}{d\theta} = \dfrac{d}{dt}\dfrac{\partial\mathcal{L}}{\partial\dot{\theta}} = \dfrac{d}{dt}mr^2\dot{\theta}
\end{align}$$
Hence, $mr^2\dot{\theta}$ is a constant, which is the angular momentum,
$L = b\sqrt{2mE}$. Therefore,
$$\begin{align}
    \dfrac{d\theta}{dt} = \dfrac{L}{mr^2}
\end{align}$$
By energy conservation, we can write
$$\begin{align}
    E = \dfrac{1}{2}m \left(\left(\dfrac{dr}{dt}\right)^2 + r^2 \left(\dfrac{d\theta}{dt}\right)^2\right) + \dfrac{e^2Z}{4\pi r}
\end{align}$$
We can readily solve this equation for the time derivative of $r$,
obtaining:
$$\begin{align}
    \dfrac{dr}{dr} = \dfrac{1}{mr}\sqrt{2mEr^2 - \dfrac{e^2Zm}{2\pi}r - L^2}
\end{align}$$
We can write
$$\begin{align}
    \dfrac{d\theta}{dt} = \dfrac{d\theta}{dr}\dfrac{dr}{dt} = \dfrac{L}{mr^2}
\end{align}$$
Therefore,
$$\begin{align}
    \dfrac{d\theta}{dr} = \dfrac{L}{mr^2}\left(\dfrac{dr}{dt}\right)^{-1} = \dfrac{L}{r\sqrt{2mEr^2 - \dfrac{e^2Zm}{2\pi}r - L^2}}
\end{align}$$
This equation can be integrated to obtain the final value of $\theta$.
We would want to integrate this equation over $\infty\to\infty$. We can
instead integrate this equation over $r_{0}\to\infty$ where $r_0$ is the
minimum value of $r$. The minimum value of $r$ will occur when
$dr/dt = 0$, which happens when
$$\begin{align}
    2mEr^2 - \dfrac{e^2Zm}{2\pi}r - L^2 = 0
\end{align}$$
We can solve this equation by completing the square:
$$\begin{align}
    r_0 = \dfrac{e^2Z}{8\pi E} + \sqrt{\dfrac{L^2}{2mE} + \dfrac{e^4Z^2}{64\pi^2 E^2}} = A + \sqrt{A^2 + b^2}
\end{align}$$
where
$$\begin{align}
    A & = \dfrac{e^2Z}{8\pi E}
\end{align}$$
Given these definitions,
$$\begin{align}
    \theta(\infty)-\theta(r_{0}) = \int_{r_{0}}^{\infty}dr \dfrac{b}{r\sqrt{(r-A)^2-A^2-b^2}}
\end{align}$$
Integrating this give us:
$$\begin{align}
    \theta(\infty)-\theta(r_{0}) = -\dfrac{\pi}{2} - i\log(\dfrac{-b+iA}{\sqrt{A^2+b^2}})
\end{align}$$
If we rotate our system such that $\theta(r_0) = \pi/2$ (which sends
$\theta\to\theta/2$), then we find that
$$\begin{align}
    \sin(\theta/2) = \dfrac{A}{\sqrt{A^2+b^2}}
\end{align}$$
We thus find that
$$\begin{align}
    b = A\cot(\theta/2)
\end{align}$$
This implies that the impact parameter as a function of $\theta$ is
$$\begin{align}
    b = \dfrac{e^2Z}{8\pi E\tan(\theta/2)}
\end{align}$$
Recall that the differential scattering cross-sections is
$$\begin{align}
    \dfrac{d\sigma}{d\Omega} = \dfrac{b}{\sin\theta}\dfrac{db}{d\theta}
\end{align}$$
we find that
$$\begin{align}
    \dfrac{d\sigma}{d\Omega} = \dfrac{e^4Z^2}{256\pi^2E^2\sin^4(\theta/2)} = \dfrac{\alpha^2Z^2}{16E^2\sin^4(\theta/2)} = \dfrac{\alpha^2Z^2}{4m^2v_i^2\sin^4(\theta/2)}
\end{align}$$
where $\alpha = e^2/4\pi$.

## Quantum Field Theory Calculation

In this section, we will investigate Rutherford Scattering through the
point of view of quantum field theory. We will consider the electron and
positron as dynamic fields and the photon field as static. The
Lagrangian we will consider is
$$\begin{align}
\mathcal{L}= \overline{\Psi}\left(i\cancel{\partial}-m\right)\Psi
    -eA_{\mu}\overline{\Psi}\gamma^{\mu}\Psi
\end{align}$$
Our goal will be to calculate the differential cross-section
$\dfrac{d\sigma}{d\Omega}$ for an electron scattering of the static
potential. The $T$ matrix element for this process, to first order in
the electromagnetic coupling constant, is

$$\begin{align}
\langle{p'}|iT|\rangle{p}
& = \bra{p',s'}T\exp\left[-i\int d^4x \mathcal{L}\_{I}(x)\right]\ket{p,s} \\\\
& = ie\int d^4x \bra{p',s'}A_{\mu}(x)\overline{\Psi}(x)\gamma^{\mu}\Psi(x)\ket{p,s}+
\mathcal{O}(e^2)\\\\
& = ie\overline{u}^{s'}(p')\gamma^{\mu}u^s(p)
\int d^4xA_{\mu}(x)e^{i(p'-p)x}
\end{align}$$

Now Fourier transform the photon field:

$$\begin{align}
    A\_{\mu}(x) = \int \dfrac{d^4q}{(2\pi)^4}\tilde{A}\_{\mu}(q)e^{-iqx}
\end{align}$$

Inserting this into the expression for the $T$-matrix element, we find

$$\begin{align}
    \langle{p'}|iT|\rangle{p}
      & = ie\overline{u}^{s'}(p')\gamma^{\mu}u^s(p)
    \int d^4x\int \dfrac{d^4q}{(2\pi)^4}
    \tilde{A}\_{\mu}(q)e^{i(p'-p-q)x}
\end{align}$$

Integrating over $x$ will produce a factor of $(2\pi)^4\delta^4(p'-p-q)$. We can then integrate over $q$ using this delta function, producing

$$\begin{align}
    \langle{p'}|iT|\rangle{p}
      & = ie\overline{u}^{s'}(p')\gamma^{\mu}u^s(p)\tilde{A}_{\mu}(p'-p)
\end{align}$$

Now, suppose the photon field is time-independent. Then, the fourier
transform of the photon field is

$$\begin{align}
\tilde{A}\_{\mu}(p'-p) & = \int d^3zA\_{\mu}(\vec{z}) e^{i(\vec{p}'-\vec{p})\cdot\vec{z}}\int_{-\infty}^{\infty}dt'e^{i(E_{p'}-E_{p})t'}\\\\
& = (2\pi)\delta(E_{p'}-E_{p})\int d^3zA_{\mu}(\vec{z})
e^{i(\vec{p}'-\vec{p})\cdot\vec{z}}
\end{align}$$

We now define the matrix element as

$$\begin{align}
    \langle{p'}|iT|\rangle{p} = (2\pi)i\delta(E_{p'}-E_{p})\mathcal{M}
\end{align}$$

Through this definition, we can see that $\mathcal{M}$ is

$$\begin{align}
    \mathcal{M}= e\overline{u}^{s'}(p')\gamma^{\mu}u^s(p)
    \tilde{A}_{\mu}(\vec{p}'-\vec{p})
\end{align}$$

where it is understood that $E_{p'} = E_{p}$. Let's square this matrix element and average over the initial spin and sum over final spin of the electron:

$$\begin{align}
\sum_{s,s'}\left|\mathcal{M}\right|^2
& =\dfrac{1}{2}e^2\tilde{A}\_{\mu}\tilde{A}\_{\nu}
   \Tr\left[\left(\cancel{p}'+m\right)\gamma^{\mu}
   \left(\cancel{p}+m\right)\gamma^{\nu}\right]\\\\
& =\dfrac{1}{2}e^2\tilde{A}\_{\mu}\tilde{A}\_{\nu}
   \left(\Tr\left[\cancel{p}'\gamma^{\mu}\cancel{p}\gamma^{\nu}\right]
    +m^2\Tr\left[\gamma^{\mu}\gamma^{\nu}\right]\right)\\\\
& = 2e^2\tilde{A}\_{\mu}\tilde{A}\_{\nu}
    \left(p'^{\mu}p^{\nu}+p'^{\nu}p^{\mu}-(p'\cdot p)g^{\mu\nu}
    +m^2g^{\mu\nu}\right)\\\\
& = 2e^2\left[2\left(\tilde{A}\cdot p'\right)\left(\tilde{A}\cdot p\right)
    +(m^2-p'\cdot p)\left(\tilde{A}\cdot\tilde{A}\right)\right]
\end{align}$$

Now, let's assume that the four-potential is generated by a charged
particle of charge $Ze$. Then, the vector potential is given by

$$\begin{align}
A^0(\vec{x}) = \dfrac{Ze}{4\pi r} \qquad \text{and}\qquad\vec{A}(\vec{x}) = 0
\end{align}$$

where $r = \left|\vec{x}\right|$. Now, let's compute the fourier transform of
this potential. To do so, we need to regulate the integrate. We modify
it in the following way:

$$\begin{align}
    A^0(\vec{x}) = \dfrac{Ze}{4\pi r}e^{-\lambda r}
\end{align}$$

where $\lambda$ will be taken to zero at the end of the calculation. The
Fourier transform of this vector potential is

$$\begin{align}
\tilde{A}^0(k)
&= \int d^3\vec{x}\dfrac{Ze}{4\pi r}e^{-\lambda r}e^{i k r\cos(\theta)}\\\\
&= \dfrac{Ze}{2}\int_{0}^{\infty}dr\int_{-1}^{1}d(\cos\theta)re^{-\lambda r}
  e^{i k r\cos(\theta)}\\\\
&= \dfrac{Ze}{2ik}\left[\dfrac{e^{(ik-\lambda) r}}{ik-\lambda}+
   \dfrac{e^{(-ik-\lambda) r}}{ik+\lambda}\right]_{0}^{\infty}\\\\
&= -\dfrac{Ze}{2ik}\left[\dfrac{(-ik-\lambda)+(-ik+\lambda)}{k^2+\lambda^2}\right]\\\\
& = \dfrac{Ze}{k^2+\lambda^2}
\end{align}$$

where $k = \|\vec{p}'-\vec{p}\|$. Taking $\lambda$ to zero, we obtain

$$\begin{align}
    \tilde{A}^0(\vec{p}'-\vec{p}) = \dfrac{Ze}{|\vec{p}'-\vec{p}|^2}
\end{align}$$

Let's use this to simplify the matrix element. The matrix element is

$$\begin{align}
    \sum_{s,s'}\left|\mathcal{M}\right|^2
      & = 2e^2
    \left[2 \dfrac{Z^2e^2E_iE_f}{|\vec{p}'-\vec{p}|^4}
        +(m^2-E_iE_f+\vec{p}\cdot\vec{p}')
        \dfrac{Z^2e^2}{|\vec{p}'-\vec{p}|^4}\right]\\\\
      & =  \dfrac{2Z^2e^4}{|\vec{p}'-\vec{p}|^4}
    \left[E_iE_f+m^2+\vec{p}\cdot\vec{p}'\right]
\end{align}$$

We are now ready to compute the differential cross section. The
differential cross section is given by

$$\begin{align}
d\sigma = \dfrac{1}{2E_iv_i}\dfrac{d^3p_f}{(2\pi)^3(2E_f)}
(2\pi)\delta(E_f-E_i)\left|\mathcal{M}\right|^2
\end{align}$$

Using

$$\begin{align}
    p_f^2 dp_fd\Omega = p_fE_fdE_fd\Omega
\end{align}$$

(where $p_f$ and $p_i$ are the magnitude of the initial and final momentum
respectively,) and integrating over $E_f$, we obtain $E_f=E_i=E$ and
$p_f = p_i$. Thus, the differential cross section is

$$\begin{align}
    \dfrac{d\sigma}{d\Omega}
    = \dfrac{1}{4Ev_i(2\pi)^2}
    \dfrac{2Z^2e^4}{|\vec{p}'-\vec{p}|^4}p_f
    \left[E^2+m^2+\vec{p}\cdot\vec{p}'\right]
\end{align}$$

Let $\theta$ be the angle between $\vec{p}'$ and $\vec{p}$. Then,

$$\begin{align}
|\vec{p}'-\vec{p}|^4 & = \left(p_i^2+p_f^2-2p_fp_i\cos\theta\right)^2 \\\\
                     & = \left(2p_i^2-2p_i^2\cos\theta\right)^2       \\\\
                     & = 4p_i^4\left(1-\cos\theta\right)^2            \\\\
                     & = 16p_i^4\sin^4(\theta/2)
\end{align}$$

Now, in the non-relativistic limit, the momentum is $p_i=mv_i$ and $E=m$.
Hence, the differential cross section is

$$\begin{align}
    \dfrac{d\sigma}{d\Omega}
      & = \dfrac{1}{4mv_i(2\pi)^2}
    \dfrac{2Z^2e^4mv_i}{(16)m^4v_i^4\sin^4(\theta/2)}
    \left[E^2+m^2+\vec{p}\cdot\vec{p}'\right]\\\\
      & = \dfrac{1}{8\pi^2}
    \dfrac{Z^2e^4}{(16)m^4v_i^4\sin^4(\theta/2)}
    \left[2m^2+m^2v_i^2\cos\theta\right]
\end{align}$$

Since $v_i\ll 1$, we find that

$$\begin{align}
    \dfrac{d\sigma}{d\Omega}
      & =\dfrac{Z^2\alpha^2}{4m^2v_i^4\sin^4(\theta/2)}
\end{align}$$

where we used $e^2=4\pi\alpha$.
