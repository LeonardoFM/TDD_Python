# TDD_Python

### Firts project for TDD learning based on seismology numerical methods


## The first note commes from Gerard Meszaros preface book - xUnit Test Patterns: refactoring  test code

"If you look at how most programmers spend their time, you’ll find that
writing code is actually a small fraction."

## The second reference from Heiner book Igel - Computational seismology: a practical introduction

"For seismologists the calculation of synthetic (or theoretical) seismograms is a key
activity on the path to a better understanding of the structure of the Earth’s
interior, or the sources of seismic energy. There are many ways of doing this,
depending in particular on the assumptions made in the geophysical model.
In the most general case—an Earth in which the properties vary in three
dimensions—analytical solutions do not exist."

#### Seismic waves:

  > f = 1Hz frequency, c = 3Km/s \approx 10.000 Km/s$ velocity

-> wavelength $\lambda$
\begin{equation}
  c = \lambda  f
\end{equation}

than the wavelength $\lambda = 3Km$

Number of points per wavelength

$$ n_{\lambda} = \frac{\lambda}{dx}$$

\subsection{Spatial scales and meshing}

For discretization we must have at least 10 grid points per wavelenth.
Now, how to discretiza the entire Earth to simulate waves with $f=1Hz$, $c = 3Km/s$?

10 grid points/$\lambda$ -> 0.3Km grid spacing

Earth:
  radius $r_E = 6.371Km$
  volume $V_E=\frac{4}{3}\pi r_E^3$ => $\frac{V_E}{(0.3Km)^3} = 4.10^{13}elements$
  8 bytes * 4 * $10^13 = 320$TBytes

Dimensionality:
  Focus: 1D

\subsection{Waves in a discrete world}

Acoustic wave equation:
\begin{equation}
  \partial^2_t p(x,t) = c(x)^2\Delta p(x,t)
  \label{wave_eq}
\end{equation}
  without sources, $c(x) = c_0$, simple solution $p(x,t)=p_0 e^{(i)kx-wt}$
When we inject $p(x,t)=p_0 e^{(i)kx-wt}$ into \ref{wave_eq} the dispersion relation come out as
\begin{equation}
  c = \frac{\omega}{k} =\frac{2\pi f}{k} = \frac{2\pi /T}{2\pi/\lambda} = \lambda f
\end{equation}

How we classify partial differential equations?
\begin{equation}
  Ap_{xx} + Bp_{xt} + Cp_{tt} + Dp_x + Ep_t + Fp = 0
\end{equation}
where A,B,C,D and E are constants, when we use Fourier trasnform
\begin{equation}
  Ax^2 + Bxt + Ct^2 + Dx + Et + F = 0
\end{equation}
arising the discriminat which depends on coefficients
\begin{equation}
  B^2 - 4AC = 0,\ (parabolic)
\end{equation}
\begin{equation}
  B^2 - 4AC < 0,\ (elliptical)
\end{equation}
\begin{equation}
  B^2 - 4AC > 0,\ (hyperbolic)
\end{equation}
than our classification is in $p_{tt} - c^2p_{xx} = 0$, so
\begin{equation}
  B^2 - 4AC = 0 - 4(-c^2)1 = 4c^2 > 0,\ (hyperbolic)
\end{equation}

With initial conditions, e.g., $p(x,t=0)$ and $\partial_t p(x,t=0)$ the solutions are fixed at all times and are wavelike.

This categories come from intersections with cones, implies the form of the partial diferential equations in the Fourier domain.

\subsection{Parallel Simulations}

Hardware architecture of supercomputers influences the chois of mathematical algorithms and their efficiences. Maicon Flinn categories:

\begin{itemize}
  \item SISD: single instruction single data (serial computer)
  \item MISD: multiple instruction single data (cryptographic decoding)
  \item SIMD: single instruction multiple data (GPU cluster, CM2)
  \item MIMD: multiple instruction multiple data (supercomputers PC cluster)
\end{itemize}

we have shared memory, distributed memory and hybrid distributed-shared memory.

Languages:
  1) compilers: fortran, c , c++
  2) Message - Passing Interface (e.g. for Python)

How efficient is the parallelization?
  $P =$ fraction of the code that can be parallelized
  $S =$ serial fraction
  $n =$ number of processes
  $$speed-up = \frac{1}{\frac{P}{n} + S}$$


Demosntrate that your algothm is efficiently scaling!

\section{A bit of wave physics}

Wave equation
\begin{equation}
  \partial^2_t p(x,t) = c(x)^2\Delta p(x,t)
\end{equation}

Analytical solutions:
  if $c(x) = c_0$ and $s(x,t) = 0$: $p(x,t=o) = p_0$ and $\partial_t p(x,t=0)$
  if $s(x,t) = \delta(x-x_0)\delta(t-t_09)$: (impulse response)-> (Green's function)

  so the solution for impulse response is a Green's function
  \begin{equation}
    \partial^2_t G(x,t;x_0,t_0) c^2\Delta G(x,t;x_0,t_0) = s(x,t) = \delta(x-x_0)\delta(t-t_09)
    \label{greens_eq}
  \end{equation}

Elastic wave equation:
  \begin{itemize}
    \item $\mathbf{u}_y$: displacement
  \end{itemize}

  \begin{equation}
    \rho \partial^2_t \mathbf{u}_y = \partial_x(\mu \partial_x \mathbf{u}_y) + f_y
    \label{ew_eq}
  \end{equation}
where $\mu$ is the shear modulus.

Is we asssume homogeneous medium, then $\rho$ and $\mu$ independent of space
\begin{equation}
  \partial^2_t \mathbf{u}_y = \frac{\mu}{\rho} \partial^2_x \mathbf{u}_y + f_y
  \label{ew_eq_homog}
\end{equation}
where $c = \sqrt{\frac{\mu}{\rho}}$

Plane wave description:
\begin{equation}
  p(\mathbf{x},t) = p_0 sin(\mathbf{kx} - \omega t)
\end{equation}
where $\mathbf{k}$ is the wavenumber vector and $\mathbf{x}$ the position vector

for elastic equation we have
\begin{equation}
  \mathbf{u}_y(\mathbf{x},t) = A_y sin(\mathbf{k_x x} - \omega t)
\end{equation}
 and $c = \omega/ |\mathbf{k}|$.


 We initialize a Gaussian function

\begin{equation}
f(x)=\frac{1}{\sqrt{2 \pi a}}e^{-\frac{(x-x_0)^2}{2a}}
\end{equation}

Note that this specific definition is a $\delta -$generating function. This means that $\int f(x)dx=1$ and in the limit $a\rightarrow 0$ the function $f(x)$ converges to a $\delta -$function.

Now let us calculate the second derivative using the finite-difference operator with three points
\begin{equation}
f^{\prime\prime}_{num}(x)=\frac{f(x+dx)-2 f(x)+f(x-dx)}{dx^2}
\end{equation}

and compare it with the analytical solution
\begin{equation}
f^{\prime\prime}(x)= \frac{1}{\sqrt{2\pi a}} ( \frac{(x-x_0)^2}{a^2}- \frac{1}{a} ) \ e^{-\frac{(x-x_0)^2}{2a}}
\end{equation}

In the cell below calculation of the first derivative with four points is provided with the following weights:
\begin{equation}
f^{\prime\prime}(x)=\frac{-\frac{1}{12}f(x-2dx)+\frac{4}{3}f(x-dx)-\frac{5}{2}f(x) +\frac{4}{3}f(x+dx)-\frac{1}{12}f(x+2dx)}{dx^2}
\end{equation}
