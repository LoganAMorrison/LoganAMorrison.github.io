---
title: Evaluating Functional Determinants in a Thermal Field Theory using the Heat-Kernel
subtitle: 
math: true

# Summary for listings and search engines
summary: 

# Link this post with a project
projects: []

# Date published
date: "2019-9-02"

# Date updated
lastmod: "2019-9-02"

# Is this an unpublished draft?
draft: false

# Show this page in the Featured widget?
featured: false

authors:
- Logan A. Morrison

tags:
- quantum-field-theory 
- finite-temperature-field-theory

categories:
- Teaching

---
## Introduction
In this post, we will discuss how to compute fluctuation determinants is a
finite temperature field theory. The general form of an operator that we will
investigate is:

$$
\begin{align}
\mathcal{O} &= -\nabla^2 + \omega_{n}^2 + m^2 + V(x)
\end{align}
$$

where $m$ is a zero-temperature mass, $\omega_{n} = 2\pi n T$ for bosons and
$2\pi(n+1/2)T$ for fermions and $V(x)$ is some space-dependent function. Often,
$V(x)$ will take the form of a background- or low-mass-field-dependent mass for
the fluctuation fields of the problem. For example, consider a $\phi^4$ theory at finite temperature with the following action:

$$
\begin{align}
S[\phi] &= \int_{0}^{\beta}d\tau\int d^{3}x \dfrac{1}{2}(\partial_{\mu}\phi)^2 -\dfrac{1}{2}m^2\phi^2
+\dfrac{\lambda}{8}\phi^4
\end{align}
$$

where $\phi = \phi(\tau, x)$. If one expands the field $\phi$ as a Fourier series in the imaginary time coordinate $\tau$,

$$
\begin{align}
\phi(\tau, x) &= \sum_{n=-\infty}^{\infty}\phi_{n}e^{i\omega_{n}\tau}
\end{align}
$$

and wishes to integrate over the heavy Matsubara modes ($n\neq0$ modes), then they will need to compute a thermal fluctuation determinant of an operator of the form:

$$
\begin{align}
\mathcal{O} &= -\nabla^2 + \omega_{n}^2 + m^2 + \dfrac{3}{2}\lambda\phi_{0}^2
\end{align}
$$

Where $\phi_{0}$ is the zero-mode ($n=0$ Matsubara mode.) In this case,
$V(x) = (3\lambda/2)\phi_{0}^2(x)$. The question we would like to answer is: how does one evaluate the following fluctuation determinant

$$
\begin{align}
\mathrm{tr}\log\left(-\nabla^2 + \omega_{n}^2 + m^2 + V(x)\right)
\end{align}
$$

We will answer this question in the limit as $T\to\infty$.

## Heat Kernal Expansion
The method we will employ to evaluate the fluctuation determinant of the operator $\mathcal{O} = -\nabla^2 + \omega_{n}^2 + m^2 + V(x)$ is the so called
heat-kernel expansion. The idea is to use the following identity:

$$
\begin{align}
\log(\lambda) = -\int_{0}^{\infty}\dfrac{ds}{s}e^{-\lambda s} + \mathrm{constant}
\end{align}
$$

where the constant we left implicit is some infinite constant which will play
no role in our analysis. To derive this expression, one can consider the
following:

$$
\begin{align}
\int_{0}^{\infty}\dfrac{ds}{s}e^{-\lambda s} &=
\int_{0}^{\infty}ds\int_{\lambda}^{\infty}d\lambda'e^{-\lambda' s}
\end{align}
$$

Switching the order of integration and integrating over $s$, we find:

$$
\begin{align}
\int_{0}^{\infty}\dfrac{ds}{s}e^{-\lambda s} &=
\int_{\lambda}^{\infty}d\lambda'\dfrac{1}{\lambda'} = \lim_{\lambda'\to\infty}\log(\lambda') - \log(\lambda)
\end{align}
$$

Thus, up to an infinite constant, the identity holds. Note that if we take
the trace of a log of a determinant, we are equivalently summing over the
log of the eigenvalues of the operator. Therefore, we can say that:

$$
\begin{align}
\mathrm{tr}\log(\mathcal{O}) = \sum_{\lambda}\log(\lambda) = -\int_{0}^{\infty}\dfrac{ds}{s}\sum_{\lambda}e^{-\lambda s} =
\int_{0}^{\infty}\dfrac{ds}{s}\mathrm{tr}e^{-\mathcal{O} s}
\end{align}
$$

We will define the factor of $e^{-\mathcal{O}s}$ as the **heat kernel**:

$$\begin{align}
K(s, x, y) \equiv \langle x|e^{-\mathcal{O}_{x} s}|y\rangle
\end{align}$$

where the subscript $x$ denotes that the derivatives and function evaluation of $\mathcal{O}$ to be evaluated with $x$. Then, the trace of $e^{-\mathcal{O} s}$ is just the integral over $x$ of $K(s,x,x)$:

$$\begin{align}
\mathrm{tr}e^{-\mathcal{O} s} = \int d^{d}x K(s,x,x)
\end{align}$$

The heat-kernel turns out to satisfy the heat equation (hence the name). This is
easy to see:

$$
\begin{align}
\dfrac{\partial}{\partial s}e^{-\mathcal{O}_{x} s} = -{\mathcal{O}\_{x}} K(s,x,y)
\quad \implies \quad \left(\dfrac{\partial}{\partial s} + {\mathcal{O}\_{x}}\right)K(s,x,y) 
= 0
\end{align}
$$

The right-most equality is simply the heat equation with a source term of
$\omega_{n}^2 + m^2 + V(x)$. The boundary condition of the heat equation can be
found by setting $s=0$:

$$
\begin{align}
K(s=0,x,y) = \langle x|y\rangle = \delta^d(x-y)
\end{align}
$$

The heat-kernel can be determined exactly if we set $V(x) = 0$. Let's call:

$$
\begin{align}
\mathcal{M}^2 = \omega_{n}^2 + m^2
\end{align}
$$

for ease of notation. The heat-kernel for the operator:

$$
\begin{align}
\mathcal{O}_{0} = -\nabla^2 + \mathcal{M}^2
\end{align}
$$

is given by the following:

$$
\begin{align}
K_{0}(s,x,y) = \dfrac{1}{(4\pi s)^{d/2}}\exp\left(-\dfrac{|x-y|^2}{4s}-s\mathcal{M}^2\right)
\end{align}
$$

This can be verified by differentiation. Additionally, one can check that this
expression satisfies the boundary condition (to do this, integrate some function
$f(x)$ by Taylor expanding the function. The results are a set of gaussian
integrals which can easily be evaluate to obtain $f(y)$.) Further more,
$K_{0}(s, x, x)$ is:

$$
\begin{align}
K_{0}(s,x,x) = \dfrac{1}{(4\pi s)^{d/2}}e^{-s\mathcal{M}^2}
\end{align}
$$

We can use this result to find the general result. Recall the Zassenhaus formula (a special case of the Baker-Campbell-Hausdorff formula):

$$
\begin{align}
e^{-t(X + Y)} &= e^{-t X}e^{-t Y}e^{-\frac{t^2}{2}[X,Y]}
e^{-\frac{t^3}{3!}(2[Y,[X,Y]]+ [X,[X,Y]])}\cdots
\end{align}
$$

where the $\cdots$ represent terms of order $t^{4}$ and higher. If we break up our general operator into:

$$
\begin{align}
\mathcal{O} = \mathcal{O}_{0} + V(x),
\end{align}
$$

then we find that the general heat-kernel is given by:

$$
\begin{align}
e^{-s\mathcal{O}} &= e^{-s\mathcal{O}_{0}}e^{-sV}e^{-\frac{s^2}{2}[\mathcal{O}\_{0},V]}
e^{-\frac{s^3}{3!}(2[V,[\mathcal{O}\_{0},V]]+ [\mathcal{O}\_{0},[\mathcal{O}\_{0},V]])}\cdots
\end{align}
$$

We can easily evaluate the various commutators by considering their action on some test function $f(x)$. The results for the commutators shown are:

$$
\begin{align}
[\mathcal{O}_{0},V] &= -\nabla^2 V\\\\
2[V,[\mathcal{O}\_{0},V]]+ [\mathcal{O}\_{0},[\mathcal{O}\_{0},V]] &=-4V\nabla^2V+\nabla^4V
\end{align}
$$

The higher order terms can easily be evaluated as well. Thus, we have that:

$$
\begin{align}
K(s,x,x) &= e^{-s(\mathcal{M}^2 + V)}e^{\frac{s^2}{2}\nabla^2V}
e^{-\frac{s^3}{3!}\left(\nabla^4V-4V\nabla^2V\right)}\cdots
\end{align}
$$

In order to make progress on this expression, we need to have control over which terms in the expression are important. Let's re-introduce our original definitions:

$$
\begin{align}
\mathcal{M}^2=\omega_{n}^2+m^2 = \tilde{\omega}_{n}^2T^2 + m^2
\end{align}
$$

where we defined:

$$
\begin{align}
\tilde{\omega}_{n}^2 = \omega\_{n}^2/T^2
\end{align}
$$

Next, we will define a scaleless variable $z=sT^2$. In terms of $z$, our expression is:

$$
\begin{align}
K(s,x,x) &= e^{-z\tilde{\omega}_{n}^2}e^{-\frac{z}{T^2}(m^2 + V)}e^{\frac{z^2}{2T^4}\nabla^2V}
e^{-\frac{z^3}{3!T^6}\left(\nabla^4V-4V\nabla^2V\right)}e^{\mathcal{O}(T^{-6})}
\end{align}
$$

Now it is clear how to proceed. We expand this expression in inverse powers of $T$. If we expand in inverse powers of $T$ and integrate over $s$, we find that, for bosons:

$$
\begin{align}
    &\dfrac{T}{2}\sum_{n=-\infty}^{\infty} \mathrm{tr}\log(-{\nabla^2}+(2\pi n T)^2 + m^2 +V({x}))\\\\
    &\hspace{1cm}= -\dfrac{T}{2}\sum_{n=\infty}^{\infty} \int d^{d}{x}\left(\dfrac{T}{4\pi}\right)^{d/2}\int_{0}^{\infty}dz\dfrac{K(z,{x},{x})}{z^{d/2+1}}\notag\\\\
    &\hspace{1cm}\underset{d\to3-2\epsilon}{=} \int d^{3}{x}\bigg{[}
        -\frac{\pi^2 T^4}{90}+\frac{T^2(m^2 + V)}{24} -\frac{\left(m^2 + V\right)^2 \left(\frac{1}{\epsilon}+\log(\frac{\mu^2}{4\pi T^2e^{-\gamma_{E}}})\right)}{64\pi^2}+ \mathrm{O}({T^{-2}})
        \notag\\\\
    &\hspace{3cm} + 
        \left(-\dfrac{T(m+V)^{3/2}}{12\pi} + \dfrac{T({\nabla}^2V)}{32\pi\sqrt{m^2+V}} + \dfrac{T(-4V{\nabla}^2V + {\nabla}^2V)}{192\pi(m^2+V)^{3/2}} + \cdots\right)
        \bigg{]}\notag
\end{align}
$$

and for fermions:

$$
\begin{align}
    &-\dfrac{T}{2}\sum_{n=-\infty}^{\infty} \mathrm{tr}\log(-{\nabla^2}+(2\pi (n+1/2) T)^2 + m^2 +V({x}))\\\\
    &\hspace{1cm}=\dfrac{T}{2}\sum_{n=\infty}^{\infty} \int d^{d}{x}\left(\dfrac{T}{4\pi}\right)^{d/2}\int_{0}^{\infty}dz\dfrac{K(z,{x},{x})}{z^{d/2+1}}\notag\\\\
    &\hspace{1cm}\underset{d\to3-2\epsilon}{=}
    \int d^{3}{x}\bigg{[} -\frac{7\pi^2 T^4}{720}
    + \frac{T^2(m^2 + V)}{48}
    + \frac{\left(m^2 + V\right)^2 \left(\frac{1}{\epsilon}+\log(\frac{4e^{-\gamma_{E}}\mu^2}{\pi T^2})\right)}{64\pi^2}+\mathrm{O}({T^{-2}})\bigg{]}\notag
\end{align}
$$

It is useful to note that, for bosons, the term proportional to $(m^2+V)^2$ can be written as:

$$
\begin{align}
&-\frac{\left(m^2 + V\right)^2 \left(\frac{1}{\epsilon}+\log(\frac{\mu^2}{4\pi T^2e^{-\gamma_{E}}})\right)}{64\pi^2}\\\\
&\qquad= \frac{\left(m^2 + V\right)^2}{64\pi^2}\left[-\left(\frac{1}{\epsilon}+\log(4\pi e^{-\gamma})\right)+2\left(\dfrac{3}{4}+\log(4\pi)-\gamma_{E}\right)+\left(\log\left(\dfrac{T^2}{\mu^2}\right)-\dfrac{3}{2}\right)\right]\notag
\end{align}
$$

where we've isolated the typical subtraction term (containing the $1/\epsilon +
\log(4\pi e^{-\gamma_{E}})$) and the Coleman-Weinberg-like potential term (the
$\log(T^2/\mu^2) - 3/2$). We can do the same for fermions, obtaining:

$$
\begin{align}
&\frac{\left(m^2 + V\right)^2 \left(\frac{1}{\epsilon}+\log(\frac{4\mu^2e^{-\gamma_{E}}}{\pi T^2})\right)}{64\pi^2}\\\\
&\qquad= \frac{\left(m^2 + V\right)^2}{64\pi^2}\left[
\left(\frac{1}{\epsilon}+\log(4\pi e^{-\gamma})\right)-2\left(\dfrac{3}{4}+\log(\pi)-\gamma_{E}\right)-\left(\log\left(\dfrac{T^2}{\mu^2}\right)-\dfrac{3}{2}\right)\right]\notag
\end{align}
$$

Even though we stopped our expansion at order $T^{0}$, one can easily see how to
include higher order term: simply cary out the expansion of the heat-kernel to
higher orders in $1/T$. For expand, the next term in the non-zero mode bosonic
expansion is:

$$
\begin{align}
\dfrac{\zeta(3)}{768\pi^2 T^2}\left[\nabla^4V+2V\nabla^2V+(m^2+V)^3\right]
\end{align}
$$

and the next term in the fermionic expansion is the same thing multiplied by 7.
