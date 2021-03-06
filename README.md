# [Project Site](http://jessebett.com/Radial-Basis-Function-USRA/)
# [Canadian Undergraduate Mathematics Conference Presentation](https://docs.google.com/viewer?url=https://github.com/jessebett/Radial-Basis-Function-USRA/raw/master/CUMC%20Presentation/presentation.pdf)

## About the USRA* Project
This project explores the use of Radial Basis Functions (RBFs) in the interpolation of scattered data in N-dimensions. It was completed Summer 2014 by Jesse Bettencourt as an NSERC-USRA student under the supervision of Dr. Kevlahan in the Department of Mathematics and Statistics at McMaster University, Hamilton, Ontario, Canada.

*[Undergraduate Student Research Awards (USRA)](http://www.nserc-crsng.gc.ca/students-etudiants/ug-pc/usra-brpc_eng.asp) are granted by the Natural Sciences and Engineering Research Council of Canada to 'stimulate interest in research in the natural sciences and engineering' and to encourage graduate studies and the pursuit of research careers in these fields.

## About this Repository
This repository contains resources and working documents associated with the project. The early stages of the project focused on reviewing the published literature on RBF interpolation. This was summarized in a presentation given at the [Canadian Undergraduate Mathematics Conference (CUMC) 2014 at Carelton University][1]. The pdf, LaTeX, and figures as well as the python script to generate figures from this presentation can be found in the [CUMC Presentation][2] folder. Following the presentation, the project shifted focus to demonstrating the [SciPy implementation of RBF interpolation][3]. The [Interpolation Demonstration][4] folder contains the python files associated with this exploration. Of note, the file [SphericalHarmonicInterpolation.py][5] demonstrates how RBFs can be used to interpolate spherical harmonics given data sites and measurements on the surface of a sphere. This folder also contains iPython notebooks from early experimentation with SciPy's RBF and a Mathematica notebook from preliminary assessment of Mathematica implementation of RBF interpolation.


[1]: http://cumc.math.ca/
[2]: https://raw.githubusercontent.com/jessebett/Radial-Basis-Function-USRA/tree/master/CUMC%20Presentation
[3]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.Rbf.html#scipy.interpolate.Rbf
[4]: https://raw.githubusercontent.com/jessebett/Radial-Basis-Function-USRA/tree/master/Interpolation%20Demonstration
[5]: https://raw.githubusercontent.com/jessebett/Radial-Basis-Function-USRA/master/Interpolation%20Demonstration/SphericalHarmonicInterpolation.py
