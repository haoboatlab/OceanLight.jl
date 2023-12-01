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

In our package, we calculate the energy proportion of the light ray that being transmitted to the water, transmission coefficient. [^2]
```math
t_{\perp}=\frac{2\sin(\theta_{t})\cos(\theta_{r})}{\sin{\theta_{t}+\theta_{r}}}
```
```math
t_{\parallel}=\frac{2\sin(\theta_{t})\cos(\theta_{r})}{\sin{\theta_{t}+\theta_{r}}\cos{\theta_{r}-\theta_{t}}}
```
When ``t_{\perp}`` is corresponding to the transmitted energy of the light ray in which the electric field, that constitute the electro magnetic wave, perpendicular to the plane-of-incident , and ``t_{\parallel}`` is corresponding to the transmitted energy of the light ray in which the electric field, that constitute the electro magnetic wave, parallels to the plane-of-incident.
To find the total transmitted energy ``t_{\perp} + t_{\parallel}``, we combine two equations, ``t_{\perp} + (-r_{\perp})=1`` and ``t_{\parallel}+r_{\parallel}=1``. Then, we normalize the total transmitted energy, so that all the values fall between 0 and 1. Therefore, 

```math
t = \frac{1}{2}\left\{\left[\frac{2\sin(\theta_{t})\cos(\theta')}{\sin(\theta'-\theta_{t})}\right]^2+\left[\frac{2\sin(\theta_{t})\cos(\theta')}{\sin(\theta'+\theta_{t})\cos(\theta'-\theta_{t})}\right]^2\right\}
```

## result  

We, then, transform the reflection angle into the azimuthal angle and polar angle in the spherical coordination. 

```math
temx = -\frac{\eta_{x}}{\sqrt{(\eta_{x})^{2}+(\eta_{y})^{2}}}
```
```math
temy = -\frac{\eta_{y}}{\sqrt{(\eta_{x})^{2}+(\eta_{y})^{2}}}
```

## Reference 
[^1]: Mobley, C. (1994). Across the Surface. *Light and Water: Radiative Transfer in Natural Waters* (pp. 155-157). Academic Press. 
[^2]: Hecht, E. (2001). The Propagation of Light. *Optics* (pp. 113-115). Addison-Wesley. 