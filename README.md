# Power-Corrections-For-Top-Quark-Pair-Transverse-Momentum-Distributions
Here we provide the power expansion of the One-Loop amplitudes for tt~ production. The regularization sheme is conventional dimensional regularization. The Renormalization schemes are the MSbar scheme for the strong coupling and the OS scheme for the top-quark mass and external wave function renormalization constants. 
The results are given in a mathematica script.

## Running the script
After cloning the repository one can put the file PowerCorrectionsTT which will load all the other packages. 

## Amplitudes
The results are 2* Re(Conjugate(M0l) * M1l) where M0l and M1l are the respective tree level and one loop amplitudes

### q(p1) + q~(p2) -> t(p3) + t~(p4)

Tree level can be called through the command

```
qQtT0l[s, t, mt]
```

OneLoop can be called through the command

```
qQtT1l[s,t,mt]
```
Phase Space point is given by 

```
num
```
