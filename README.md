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
> Interpolation Condition: 
\begin{equation*}
s(x_i)=f_i 
\end{equation*}
$\forall i\in\{0 ... N \}$


This is how interpolation differs from approximation, where approximation does not necessitate that our function exactly equals our measurements at the data sites. This can be achieved through different methods, e.g., Least Squares approximation. Sometimes, when accuracy at data sites is not necessary, approximation is preferred over interpolation because it can provide a 'nicer' function which could better illustrate the relationship among the data sites and measurements. For instance, approximation is heavily utilized in experimental science where measurements can contain a measurement error associated with experimental procedures. In this environment, the interpolation condition may be undesirable because it forces the interpolation to match exactly with potential measurement error, where approximation may alleviate error influence and illustrate measured correlations better. 
![enter image description here][7]
For the purposes of this project, we focus on interpolation only. 

###Interpolation Assumption
Many interpolation methods rely on the convenient assumption that our interpolation function, $s(x)$, can be found through a linear combination of **basis functions**, $\psi_i(x)$.

>Linear Combination Assumption: 
\begin{equation*}
s(x)=\sum_{i=1}^N \lambda_i \psi_i
\end{equation*}


This assumption is convenient as it allows us to utilize solving methods for systems of linear equations from linear algebra to find our interpolation function. As such, we can express our interpolation problem as a linear system.

>Interpolation as Linear System: 
\begin{equation*}
A\boldsymbol{\lambda}=\boldsymbol{f}
\end{equation*}


Where $\boldsymbol{f}$ is the vector of datasite measurements $\left[ f_1, ..., f_N \right]^T$ and $\boldsymbol{\lambda}$ is the vector of linear combination coefficients $\left[ \lambda_1, ..., \lambda_N \right]^T$.

For a system with N measurement data sites,  $A$ is an NxN-matrix called the **interpolation matrix** or **Vandermonde Matrix** The elements of A are given by the basis functions, $\psi_j$  evaluated at each data site, $x_i$. 
>Elements of $A$: 
\begin{equation*}
a_{ij}=\psi_j(x_i)
\end{equation*}

By using numerical methods and solving this linear system, we will have our interpolation function as a linear combination of our basis functions. 
####Familiar Example of Interpolation Basis
A choice of basis functions, $\psi_i$, which may familiar to undergraduate students is the basis of (N-1)-degree polynomials. If we wish to find a 1-Dimensional interpolation function from N distinct data sites, we can find an (N-1)-degree polynomial which goes exactly through all sites. In other words, by choosing our basis functions to be successive powers of x up to (N-1), we can solve our interpolation system for our function.
> Polynomial Interpolation Basis: 
\begin{equation*}
\psi_{i=1}^N=\{1,x,x^2,x^3, ..., x^{N-1}\}
\end{equation*}

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

Through the work of Alfréd Haar and his description of **Haar Spaces** we gain the negative result that well-posedness is not guaranteed in higher dimensional linear systems with independently chosen basis functions. To state the theorem we first define Haar spaces.

>Definition of Haar Space:
Let $\Omega \subset \mathbb{R}^N$ be a set with at least $N$ sites in it. Let $V \subset C(\Omega)$ be an $N$-dimensional subspace of continuous functions. Then, we say that $V$ is a **Haar Space** if for any collection of $N$ sites $\{x_1,...,x_N\}$ with any corresponding set of values $\{f_1,...,f_N\}$, we can find a unique function $s \in V$ such that $s(x_k)=f_k$.

From this definition we have the following lemma.
>Lemma:
Let $\Omega \subset \mathbb{R}^N$ be a set with at least $N$ sites in it and $V \subset C(\Omega)$ be a subspace. 
Then, $V$ is an $N$-dimensional Haar Space if and only if for any distinct sites $\{x_1,...,x_N\} \in \Omega$ and any basis of functions $\{\psi_1,...,\psi_N\} \in V$, we have $\det(\psi_j(x_i))\neq0$.

In other words:
>$V$ is a Haar Space $\iff$ any set of basis function produce well-posed system for any set of distinct data sites.

For the purposes of interpolation, then, interpolating within a Haar space is ideal, because then we can choose our basis independently of our data and, as per the lemma, we are guaranteed a well-posed system and a unique solution.

However, by the negative result of the Mairhuber-Curtis Theorem, there can be no Haar Spaces in $N$-Dimensions for $N \geq 2$
>Mairhuber-Curtis Theorem:
Let $\Omega \subset \mathbb{R}^N$, $N \geq 2$ contain an interior site. Then, there is no Haar space of dimension $N \geq 2$ for $\Omega$.

So, if we are interpolating scattered data in higher dimensions, by the Haar-Mairhuber-Curtis Theorem we cannot choose our basis functions independent from our data sites. However, this does not mean we cannot interpolate in higher dimensions using our interpolation assumption.
>If we can't guarantee well-posedness with independently chosen basis functions, we must choose our basis functions depending on our data sites. 

###Basis Functions for Higher Dimension Interpolation
One method for defining basis functions depending on our data sites is to take a single function and translate it for each site. That is, our basis functions will be **translates** of a single function for each data site. 

####Translates of the Basic Function
If our basis functions are translates of a function, which function should we translate? By answering this question we will arrive at the definition of Radial Basis Functions, but first let's consider a preliminary function: the basic function.
>The Basic Function: 
\begin{equation*}
\psi_i(x)=||x-x_i||
\end{equation*}


Pictured below, the basic function is the absolute valued function given by the Euclidean distance from a **center point** $x_i \in \mathbb{R}^N$. The basic function has the feature that it is radially symmetric about this center point. 

![enter image description here][10]

We can define our set of basis functions, $\{\psi_i(x)\}_{i=1}^N$, as translates of our basic function such that the center points are located at our data sites. 
>Set of Basis Functions: 
\begin{equation*}
\{\psi_i(x)=||x-x_i||\}_{i=1}^N
\end{equation*}

In other words, our set of basis functions is composed of basic functions centered each of our data sites. We can visualize one of these centered basic functions in the figure below.

![enter image description here][11]

Now that we have chosen our basis functions, we can look at the linear system which it produces. For instance, our interpolation matrix, $A$, now becomes:
\begin{equation*}
A=
\begin{bmatrix}
||x_1-x_1|| & ||x_1-x_2|| & \cdots & ||x_1-x_N||\\
||x_2-x_1|| & ||x_2-x_2||& \cdots & ||x_2-x_N||\\
\vdots & \vdots & \ddots & \vdots\\
||x_N-x_1|| & ||x_N-x_2||& \cdots & ||x_N-x_N||
\end{bmatrix}
\end{equation*}

This matrix is known as the **distance matrix** with Euclidean distance. 

But, if we're not interpolating in 1-D, then we know we're not in a Haar Space. How do we know that our linear interpolation system with the distance matrix is well-posed? 

>Lemma from Linear Algebra: 
Distance Matricies with Euclidean distance, for distinct points in $\mathbb{R}^n$ are always non-singular.

From the above lemma we know that our interpolation matrix is non-singular. Therefore, we know our system is well-posed and that there exists a unique interpolation function! 

However, the choice of $\psi_i(x)=||x-x_i||$ as our basic function is not ideal. As we can see from our above plot of the basic function centered at $x_i$, the first derivative of the basic function is discontinuous at our center point, $x_i$. This has the consequence that, at each of our data sites, the first derivative of our interpolation function will be discontinuous. This is problematic because, ideally, we would like to have a $C^\infty$ smooth interpolation function so we can use methods from calculus to analyze our function. 

How can we remedy our derivative discontinuities in our interpolation function?

####Building a Better Basic Function
In 1968, R.L. Hardy suggested that by using a $C^\infty$ smooth function as our basic function, we can produce smooth interpolation functions. These functions are called **Kernels**. The kernel suggested by Hardy was the **Multiquadric Kernel**.

>Hardy's Multiquadric Kernel:
\begin{align*}
\psi(x)=\sqrt{c^2 + x^2}
\end{align*}
where $c \neq 0$.

Notice that if we allow $c=0$ in the multiquadric kernel then we are actually describing the basic function used above. So, in other words, Hardy's multiquadric kernel is like the basic function but smoothed with a parameter $c$. By looking at a plot of the multiquadric kernel, we can see that the discontinuity from the basic function has been addressed. In fact, the multiquadric function is, as desired, $C^\infty$ smooth.

![enter image description here][12]

As before, we will define our basis functions, $\{\psi_i\}_{i=1}^{N}$, as a set of multiquadric kernels translated such that they are centered at our data sites, $x_i$.

>Basis of Multiquadric Kernels:
\begin{equation*}
\{\psi_i(x)=\sqrt{c^2 + (||x-x_i||)^2}\}_{i=1}^N
\end{equation*}

We can visualize one of these translated multiquadric kernels in the figure below. 

![enter image description here][13]

####Radial Basis Function Kernels
Notice that the multiquadric kernel is also radially symmetric about its center, $x_i$. Because of this radial symmetry, the multiquadric kernel can be described as a **Radial Basis Function**. In other words, it is a basis function which depends only on the radial distance from its center. Since our basis functions $\psi_i(x)$ depend only on distance, we can re-express them as such.
>Radial Basis Functions:
\begin{equation*}
\psi(||x-x_i||)= \phi(r)
\end{equation*}
where $r=||x-x_i||$

With our interpolation assumption, we can express our interpolation function as a linear combination of these functions, as before:
>Interpolation as Linear Combination of Radial Basis Functions:
\begin{equation*}
s(x)=\sum_{i=1}^N \lambda_i \psi(||x-x_i||)=\sum_{i=1}^N \lambda_i \phi(r) 
\end{equation*}

There are a few commonly used radial basis function kernels:

 - Multiquadric: $\phi(r)=\sqrt{1+(\epsilon r)^2}$![enter image description here][14]
 - Inverse Multiquadric: $\phi(r)=\frac{1}{\sqrt{1+(\epsilon r)^2}}$ ![enter image description here][15]
 - Inverse Quadratic: $\phi(r)=\frac{1}{1+(\epsilon r)^2}$ ![enter image description here][16]
 - Gaussian:  $\phi(r)=e^{-(\epsilon r)^2}$ ![enter image description here][17]

As before, we can use translates of these functions centered on our data sites as basis for our interpolation linear system. Further, notice that the multiquadric kernel has been rearranged to replace $c$ with a **shape parameter**, $\epsilon$ consistent with the other kernels.

However, by using the radial basis kernels as our basis, we change the interpolation matrix so that it is no longer the distance matrix as before.

>Interpolation matrix with RBF kernels:
\begin{equation*}
A=
\begin{bmatrix}
\phi_1(r_1) & \phi_1(r_2) & \cdots & \phi_1(r_N)\\
\phi_2(r_1) & \phi_2(r_2)& \cdots & \phi_2(r_N)\\
\vdots & \vdots & \ddots & \vdots\\
\phi_N(r_1) & \phi_N(r_2)& \cdots & \phi_N(r_N)
\end{bmatrix}
\end{equation*}

If $N\geq2$ then we are still not interpolating in a Haar Space, and since we are no longer using a distance matrix, can we expect well-posedness?

To answer this question we determine if our interpolation matrix is **positive-definite**.
>A matrix, $A$, is positive-definite if
\begin{align*}
& t^TAt>0 & \forall t=\left[ t_1, t_2, ..., t_n\right]\neq 0 \in \mathbb{R}^n
\end{align*}

Using this definition we have the following condition:
>If interpolation matrix, $A$, is symmetric, positive-definite , then $A$ is nonsingular and our system is well-posed.

So we can guarantee the existence of a unique solution if we choose our kernels such that $A$ will be positive-definite. In fact, we can produce positive-definite interpolation matrices by using positive-definite kernels.

>A function, $\phi: \mathbb{R}^n\times \mathbb{R}^n \rightarrow \mathbb{R}$, is said to be positive definite if
:
\begin{align*}
&\sum_{i=1}^N \sum_{j=1}^N \phi(||x-x_i||)t_i\bar{t_j}>0 &\forall t=\left[ t_1, t_2, ..., t_n\right]\neq 0 \in \mathbb{C}^n
\end{align*}

Of the common RBF kernels described above, all are positive-definite except Hardy's Multiquadric kernel. However, the multiquadric kernel is guaranteed to produce well-posed systems for other, similar reasons (that it is conditionally negative-definite). With the exception of Hardy's multiquadric kernel, by using positive-definite kernels we can produce positive-definite interpolation matrices which guarantee well-posed systems! 

So, by using Radial Basis Kernels for interpolation, we have shown that there exists a unique interpolation function $s(x)$ which interpolates scattered data in N-dimensions.
>Radial Basis Interpolation
\begin{equation*}
s(x)=\sum_{i=1}^N \lambda_i \psi(||x-x_i||)=\sum_{i=1}^N \lambda_i \phi(r) 
\end{equation*}

####Well-posed v.s. Well-conditioned
In the discussion above we have shown that radial basis interpolation is well-posed, so there exists a unique solution for the interpolation problem. However, because these systems are solved using numerical methods on computers, they are subject to computational limitations. By using computational methods we introduce a complication, just because a solution exists, doesn't mean that it is accessible through numerical methods. A common example of the limitations that can cause a solution to be inaccessible is the accumulation of rounding errors. If our solution exists, and the system behaves 'nicely' with computational solving methods, then we say the solution is **well-conditioned**.

Radial basis interpolation problems, although well-posed, have the propensity to be very ill-conditioned. This is in part due the choice **shape parameter**, $\epsilon$. For some systems, small changes in $\epsilon$ may have potentially significant influences on the system. 

In the two figures below we can see how increasing the value of epsilon will change the shape of the individual kernel basis functions.

For $\epsilon=0.4$
![enter image description here][18]

For $\epsilon=1$
![enter image description here][19]

In the three figures below, we can see how increasing the value of epsilon will cause the interpolation system to become ill-conditioned. Keep in mind that the interpolation solution for each espilon value still exists, but the computation methods create noise and are unable to find the function.

![enter image description here][20]

![enter image description here][21]

![enter image description here][22]


  [1]: http://cumc.math.ca/2014/
  [2]: https://github.com/jessebett/USRA/tree/master/CUMC%20Presentation
  [3]: http://scipy.org/docs/scipy/reference/generated/scipy.interpolate.Rbf.html#scipy.interpolate.Rbf
  [4]: https://github.com/jessebett/USRA/tree/master/Interpolation%20Demonstration
  [5]: https://github.com/jessebett/USRA/blob/master/Interpolation%20Demonstration/SphericalHarmonicInterpolation.py
  [6]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/interpdef.png
  [7]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/interpvsapprox.png
  [8]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/polyinterp.png
  [9]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/HMC.png
  [10]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/basicfunxi.png
  [11]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/basicbasis.png
  [12]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/kernelfun.png
  [13]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/kernelbasis.png
  [14]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/multiquadric.png
  [15]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/inversemultiquadric.png
  [16]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/inversequadratic.png
  [17]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/gaussian.png
  [18]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/basisgaus1.png
  [19]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/basisgaus2.png
  [20]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/conditioned.png
  [21]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/illconditioned.png
  [22]: https://github.com/jessebett/USRA/blob/master/CUMC%20Presentation/Figures/veryillconditioned.png
