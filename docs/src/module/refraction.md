# Refraction 

Photons refraction between air and water 

## Snell's Law 
angle of reflection [^1]
```math
\theta_{r} = \cos^{-1}|\hat{\xi}'\cdot\hat{n}|
```
When our incoming photon is directly downward: ``\hat{\xi}=\begin{bmatrix}0\\0\\1\end{bmatrix}`` and the normal vector can be defined by
``\hat{n}=\dfrac{1}{\sqrt{1+\left(\eta_{x}\right)^{2}+\left(\eta_{y}\right)^{2}}}\left(-\eta_{x}\hat{i}-\eta_{y}\hat{j}+\hat{k}\right)``
When, ``\eta_{x}`` is the partial derivative the partial derivative of ``\eta`` with respect to x: ``\dfrac{\partial\eta}{\partial x}``
and ``\eta_{y}`` is the partial derivative the partial derivative of ``\eta`` with respect to y: ``\dfrac{\partial\eta}{\partial y}``

Therefore, the angle of reflection, in this module, can be described by 
```math
\theta_{r} = \cos^{-1}\left(\frac{1}{\sqrt{1+\left(\eta_{x}\right)^{2}+\left(\eta_{y}\right)^{2}}}\right)
```

angle of transmission
```math
\theta_{t} = \sin^{-1}\left(\frac{1}{n_{w}}\sin\theta_{r}\right)
```
Substitute the ``\theta_{r}`` that we found above. 
```math
\theta_{t} = \sin^{-1}\left(\frac{1}{n_{w}}\sqrt{\frac{\left(\eta_{x}\right)^{2}+\left(\eta_{y}\right)^{2}}{1+\left(\eta_{x}\right)^{2}+\left(\eta_{y}\right)^{2}}}\right)
```


## Fresnel Reflectance 

In our package, we calculate the energy proportion of the light ray that being transmitted to the water, but first we identify the amplitude transmission coefficient or the ratio between electric field amplitude of the transmitted light ray to the intirial light ray. [^2]
```math
t_{\perp}=\left\(\frac{E_{t}}{E_{0}}\right\)_{\perp}=\frac{2\sin(\theta_{t})\cos(\theta_{r})}{\sin{\theta_{t}+\theta_{r}}}
```
```math
t_{\parallel}=\left\(\frac{E_{t}}{E_{0}}\right\)_{\parallel}=\frac{2\sin(\theta_{t})\cos(\theta_{r})}{\sin{\theta_{t}+\theta_{r}}\cos{\theta_{r}-\theta_{t}}}
```
When ``t_{\perp}`` is corresponding to the amplitude transmission coefficient of the light ray in which the electric field, that constitute the electro magnetic wave, perpendicular to the plane-of-incident , and ``t_{\parallel}`` is corresponding to the amplitude transmission coefficient of the light ray in which the electric field, that constitute the electro magnetic wave, parallels to the plane-of-incident.

The intensity of the light is proportional to the square of the amplitude, ``I \propto E^{2}``. 
All conceivable azimuths of waves that are polarized combine to form natural or unpolarized light. Each wave can be broken down into its constituent parts. Each component will have an equal amount due to symmetry. Then, half of the amplitude transmission coefficient yields the transmission coefficient of a surface in natural light. [^3]

```math
t = \frac{I_{t}}{I_{0}} = \frac{1}{2}\left\{\left[\frac{2\sin(\theta_{t})\cos(\theta_{r})}{\sin(\theta_{r}+\theta_{t})}\right]^2+\left[\frac{2\sin(\theta_{t})\cos(\theta_{r})}{\sin(\theta_{r}+\theta_{t})\cos(\theta_{r}-\theta_{t})}\right]^2\right\}
```

## result  

From the principle of least time, the path in which light travels is the path that can traveled in the least time. Hence, the incident, transmitted light ray, and the normal vector to the water surface lie in the same plane. The normal vector to the plane can be described by the cross product between incident light ray and normal vector to the water surface. When the light ray travels downward, along the z-axis, the z-axis lies on the plane, and hence, the plane is perpendicular to the xy plane. The transmitted light ray lies on the opposite side of the normal vector to the water surface. The normal vector to the water surface can be described by, ``\hat{n}_{xy}=\frac{1}{\sqrt{\eta_{x}^2+\eta_{y}^2}}\begin{bmatrix}-\eta_{x}\\-\eta_{y}\\0\end{bmatrix}```, and the angle of this vector to the x-axis is, 

```math 
\phi = \arccos{\frac{\eta_{x}}{\sqrt{(\eta_{x})^{2}+(\eta_{y})^{2}}}}
```

The range for arccos function is ``[0,\pi]``, however, our expected azimuthal angle can span from ``[0,2\pi]``. Hence, in order to find the azimuthal angle of the transmitted light ray, we need to consider every possible case of ``\eta_{x}`` and ``\eta_{y}``.

When ``\eta_{x}`` is equal to zero, the possible value for ``\phi`` is either ``\frac{\pi}{2}`` or ``\frac{3\pi}{2}``, and the value would be determined by ``\eta_{y}``.

When ``\eta_{y}`` is positive, the normal vector to the water surface ``\hat{n}`` lies on either the first or second quadrant, and the azimuthal angle of transmitted light ray would lie on either third or fourth qudrant. Therfore, 

```math 
\varphi = \pi+\arccos{\frac{\eta_{x}}{\sqrt{(\eta_{x})^{2}+(\eta_{y})^{2}}}}
```

When ``\eta_{y}`` is negative, the normal vector to the water surface ``\hat{n}`` lies on either the third or fourth quadrant, and the azimuthal angle of transmitted light ray would lie on either first or second qudrant. Therfore, 

```math 
\varphi = \pi-\arccos{\frac{\eta_{x}}{\sqrt{(\eta_{x})^{2}+(\eta_{y})^{2}}}}
```

The result azimuthal angle of transmitted ray can be described by,

```math
\varphi = \left\{
\begin{array}{ll}
\pi+\arccos\left(\frac{\eta_{x}}{\sqrt{(\eta_{x})^{2}+(\eta_{y})^{2}}}\right) & \text{for } \eta_{y} > 0 \\[10pt]
\pi-\arccos\left(\frac{\eta_{x}}{\sqrt{(\eta_{x})^{2}+(\eta_{y})^{2}}}\right) & \text{for } \eta_{y} < 0
\end{array}
\right.
```

## Reference 
[^1]: Mobley, C. (1994). Across the Surface. *Light and Water: Radiative Transfer in Natural Waters* (pp. 155-157). Academic Press. 
[^2]: Hecht, E. (2001). The Propagation of Light. *Optics* (pp. 113-115). Addison-Wesley. 
[^3]: Sears, F.W. (1949). Polarization *Optics* (pp. 173-174). Addison-Wesley