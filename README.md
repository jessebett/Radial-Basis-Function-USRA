### About the USRA* Project
This project explores the use of Radial Basis Functions (RBFs) in the interpolation of scattered data in N-dimensions. It was completed Summer 2014 by Jesse Bettencourt as an NSERC-USRA student under the supervision of Dr. Kevlahan in the Department of Mathematics and Statistics at McMaster University, Hamilton, Ontario, Canada. 

*[Undergraduate Student Research Awards (USRA)](http://www.nserc-crsng.gc.ca/students-etudiants/ug-pc/usra-brpc_eng.asp) are granted by the Natural Sciences and Engineering Research Council of Canada to 'stimulate interest in research in the natural sciences and engineering' and to encourage graduate studies and the pursuit of research careers in these fields.

### About this Repository 
This repository contains resources and working documents associated with the project. The early stages of the project focused on reviewing the published literature on RBF interpolation. This was summarized in a presentation given at the [Canadian Undergraduate Mathematics Conference (CUMC) 2014 at Carelton University][1]. The pdf, LaTeX, and figures as well as the python script to generate figures from this presentation can be found in the [CUMC Presentation][2] folder. Following the presentation, the project shifted focus to demonstrating the [SciPy implementation of RBF interpolation][3]. The [Interpolation Demonstration][4] folder contains the python files associated with this exploration. Of note, the file [SphericalHarmonicInterpolation.py][5] demonstrates how RBFs can be used to interpolate spherical harmonics given data sites and measurements on the surface of a sphere. This folder also contains iPython notebooks from early experimentation with SciPy's RBF and a Mathematica notebook from preliminary assessment of Mathematica implementation of RBF interpolation. 

#Radial Basis Function Interpolation
###What is interpolation?
Given a set of **measurements** $\{f_i\}_{i=1}^N$ taken at corresponding **data sites** $\{x_i\}_{i=1}^N$ we want to find an **interpolation function** $s(x)$ that informs us on our system at locations different from our data sites. 
![enter image description here][6]
Further, we want our function, $s(x)$ to satisfy what's called the **interpolation condition** which is that we want to our interpolation function to exactly match our measurements at our data sites. 
> Interpolation Condition: $s(x_i)=f_i$ $\forall i\in\{0 ... N \}$

This is how interpolation differs from approximation, where approximation does not necessitate that our function exactly equals our measurements at the data sites. This can be achieved through different methods, e.g., Least Squares approximation. Sometimes, when accuracy at data sites is not necessary, approximation is preferred over interpolation because it can provide a 'nicer' function which could better illustrate the relationship among the data sites and measurements. For instance, approximation is heavily utilized in experimental science where measurements can contain a measurement error associated with experimental procedures. In this environment, the interpolation condition may be undesirable because it forces the interpolation to match exactly with potential measurement error, where approximation may alleviate error influence and illustrate measured correlations better. 
![enter image description here][7]
For the purposes of this project, we focus on interpolation only. 

###Interpolation Assumption
Many interpolation methods rely on the convenient assumption that our interpolation function, $s(x)$, can be found through a linear combination of **basis functions**, $\psi_i(x)$.

>Linear Combination Assumption: $s(x)=\sum_{i=1}^N \lambda_i \psi_i$

This assumption is convenient as it allows us to utilize solving methods for systems of linear equations from linear algebra to find our interpolation function. As such, we can express our interpolation problem as a linear system.

>Interpolation as Linear System: $A\boldsymbol{\lambda}=\boldsymbol{f}$

Where $\boldsymbol{f}$ is the vector of datasite measurements $\left[ f_1, ..., f_N \right]^T$ and $\boldsymbol{\lambda}$ is the vector of linear combination coefficients $\left[ \lambda_1, ..., \lambda_N \right]^T$.

For a system with N measurement data sites,  $A$ is an NxN-matrix called the **interpolation matrix** or **Vandermonde Matrix** The elements of A are given by the basis functions, $\psi_j$  evaluated at each data site, $x_i$. 
>Elements of $A$: $a_{ij}=\psi_j(x_i)$

By using numerical methods and solving this linear system, we will have our interpolation function as a linear combination of our basis functions. 
####Familiar Example of Interpolation Basis
A choice of basis functions, $\psi_i$, which may familiar to undergraduate students is the basis of (N-1)-degree polynomials. If we wish to find a 1-Dimensional interpolation function from N distinct data sites, we can find an (N-1)-degree polynomial which goes exactly through all sites. In other words, by choosing our basis functions to be successive powers of x up to (N-1), we can solve our interpolation system for our function.
> Polynomial Interpolation Basis: $\psi_{i=1}^N=\{1,x,x^2,x^3, ..., x^{N-1}\}$

An example of this interpolation with 6 data sites can be seen in the figure below. Here the interpolation function, as a linear combination, is $s(x)=-0.02988 x^5 + 0.417 x^4 - 2.018 x^3 + 3.694 x^2 - 1.722 x - 5.511e^{-14}$
![enter image description here][8]
However, while polynomial basis is simple for 1-Dimensional interpolation, this method is not ideal for higher dimensions. To accommodate higher dimension interpolation, we must choose our basis differently.

###Well-Posedness in Higher Dimensional Interpolation
When defining our linear system we must consider whether our system is **well-posed**. That is, does there exist a solution to our interpolation problem, and if so is that solution unique? 
>Well-Posedness in Linear Systems: Our system will be well-posed if and only if $A$ is non-singular, i.e. $\det(A)\neq0$

For 1-D interpolation, many choices in basis functions will guarantee a well-posed system. In our example of polynomial interpolation, for instance, it was guaranteed that for N-distinct data sites a unique (N-1)-degree polynomial will interpolate the measurements. So without predetermining any information about our data sites (other than that they are distinct from each other), or their measurements, we can define our basis functions independently of our data and expect a unique, well-posed solution. 

However, for n-Dimensions where $n\geq2$ this is never guaranteed! That is, no matter what we choose for our set of basis functions, there will always be data sites which produce ill-posed systems. The implication of this is that we can not define our basis functions independently of our data and expect a well-posed system. This results from the Haar-Mairhuber-Curtis Theorem.

####Haar-Mairhuber-Curtis Theorem
![enter image description here][9]


  [1]: http://cumc.math.ca/2014/
  [2]: https://github.com/jessebett/USRA/tree/master/CUMC%20Presentation
  [3]: http://scipy.org/docs/scipy/reference/generated/scipy.interpolate.Rbf.html#scipy.interpolate.Rbf
  [4]: https://github.com/jessebett/USRA/tree/master/Interpolation%20Demonstration
  [5]: https://github.com/jessebett/USRA/blob/master/Interpolation%20Demonstration/SphericalHarmonicInterpolation.py
  [6]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/interpdef.png
  [7]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/interpvsapprox.png
  [8]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/polyinterp.png
  [9]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/HMC.png
