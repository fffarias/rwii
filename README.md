# Reduced wavefield imaging and inversion (rwii)

<div style="text-align: justify"> 
Reproducing the algorithm proposed in the paper: Efficient snapshot-free reverse time migration and computation of multiparameter gradients in full-waveform inversion (https://doi.org/10.1190/geo2020-0606.1) using the SeisAcoustic package (https://github.com/SeismicJulia/SeisAcoustic.jl).

<br />  

In this work the authors propose an approximation to calculate the RTM migration (and the FWI gradient), without having to store in memory or write to disk the forward wavefield that must be correlated with the reverse propagated wavefield, avoiding significant input/output (I/O) cost.

Their approach boils down to the following equation:


<img src="https://latex.codecogs.com/gif.latex?\Gamma(x,y)&space;\approx&space;\frac{1}{2&space;\alpha}&space;(\Gamma_w(x,y)&space;-&space;\Gamma_u(x,y))" title="\Gamma(x,y) \approx \frac{1}{2 \alpha} (\Gamma_w(x,y) - \Gamma_u(x,y))" />

where

<img src="https://latex.codecogs.com/gif.latex?\Gamma_u(x,y)&space;=&space;\int&space;\Phi&space;(u(x,y,t),u(x,y,t))dt." title="\Gamma_u(x,y) = \int \Phi (u(x,y,t),u(x,y,t))dt." />

and

<img src="https://latex.codecogs.com/gif.latex?\Gamma_w(x,y)&space;=&space;\int&space;\Phi&space;(w(x,y,t),w(x,y,t))dt." title="\Gamma_w(x,y) = \int \Phi (w(x,y,t),w(x,y,t))dt." />

In RTM, <img src="https://latex.codecogs.com/gif.latex?u" title="u" /> is a source wavefield forward propagated using a smooth velocity model, whereas <img src="https://latex.codecogs.com/gif.latex?w" title="w" /> equals the weighted sum of the forward <img src="https://latex.codecogs.com/gif.latex?u" title="u" /> and the backward <img src="https://latex.codecogs.com/gif.latex?\lambda" title="\lambda" /> wavefields

<img src="https://latex.codecogs.com/gif.latex?w(x,y,t)&space;=&space;u(x,y,t)&space;&plus;&space;\alpha&space;\lambda(x,y,t)" title="w(x,y,t) = u(x,y,t) + \alpha \lambda(x,y,t)." />

Whe <img src="https://latex.codecogs.com/gif.latex?w" title="w" /> wavefield is calculated by depropagating the forward wavefield <img src="https://latex.codecogs.com/gif.latex?u" title="u" /> which must be carefully stored on the boundaries of the model at each time step and together with the correct initial conditions, it is possible to reconstruct the direct wavefield within  machine precision. When this happens, the remaining term in equation above (<img src="https://latex.codecogs.com/gif.latex?\alpha&space;\lambda" title="\alpha \lambda" />), the residual propagation, gives the approximation to the RTM after scaled by <img src="https://latex.codecogs.com/gif.latex?\alpha" title="\alpha" />. According to Figure 5 of the paper, the best alpha value depends on the image condition used.

In the examples shown here the following image condition was used:

<img src="https://latex.codecogs.com/gif.latex?\Phi&space;(x,y,t)&space;=&space;u(x,y,t)&space;\lambda(x,y,t)" title="\Phi (x,y,t) = u(x,y,t) \lambda(x,y,t)" />

The difference between conventional RTM and RWII migration can be seen in the example for a <a href="https://github.com/fffarias/rwii/blob/main/examples/layers.ipynb">layered model</a> for double precision and for the <a href="https://github.com/fffarias/rwii/blob/main/examples/marmousi.ipynb">marmousi model</a> for single precision.

</div>