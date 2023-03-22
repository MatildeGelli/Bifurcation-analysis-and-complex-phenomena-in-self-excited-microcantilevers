# Bifurcation-analysis-and-complex-phenomena-in-self-excited-microcantilevers
### Authors: Matilde Gelli, Joao Mouro, Paolo Paoletti, Bruno Tiribilli and Michele Basso.
### Email: for any issue with the code contact matilde.gelli@unifi.it (or matilde.gelli@gmail.com)

Online supplementary material of the paper titled:

**Bifurcation analysis and complex phenomena in self-excited microcantilevers**

by M. Gelli[^1], J. Mouro[^2], P. Paoletti[^3], B. Tiribilli[^2], M. Basso[^1]

[^1]: Department of Information Engineering, University of Florence
[^2]: Institute for Complex Systems, National Research Council (ISC-CNR)
[^3]: School of Engineering, University of  Liverpool

Under review

MATLAB version used 2020a

### How to use it

The present work exploits DDE-Biftool, a third party package able to make continuation of bifurcations of DDEs in one and two parameters (link for the download at https://sourceforge.net/projects/ddebiftool/). Be sure to add DDE-Biftool to the path. 

### Code
* **cantilever_demo** is the main script where the low and high frequency branches are computed. Refer to this demo for getting figure 2, 4a, 5 and 7.
* **demo_Fig8** contains the code for reproducing the bifurcation diagrams in tau-Q plane. For the sake of simplicity, we provide the initial branches (the .mat files) from which the bifurcations in the two parameters are evaluated. 
