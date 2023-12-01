# Fundamental Principle of Monte Carlo Simulation

```math
p_{\Re} (\Re) \equiv \begin{cases}
1 & \quad 0 \leq \Re \leq 1, \\
0 & \quad \Re < 0 \cup \Re < 1.
\end{cases}
```

## The Optical Path Length

```math
r=-\frac{1}{c}\ln{\Re}
```

## Sampling Scattering direction

azimuthal angle in the plane of the scattering event relative to the direction of photons before scattering
``\varphi`` is uniformly distribute between 0 and ``2\pi`` Therefore,
```math
\varphi = 2\pi\Re
```
To find the angle between new trajectory and the direction of photons before scattering, we use the Petzold..

```@docs
phasePetzold()
```