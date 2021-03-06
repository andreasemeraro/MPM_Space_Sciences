
## MPM_Space_Sciences ##

## INTRODUCTION ##
This repository is made as a result of a project during the [Master in Mathematical and Physical Methods for Space Sciences](https://mpmss.i-learn.unito.it/). 
In the following part a brief description of each code will be made:
   - [Celestial mechanics and astrodynamics](https://github.com/andreasemeraro/MPM_Space_Sciences/tree/main/Celestial%20mechanics%20and%20astrodynamics):
      - [orbit.py](https://github.com/andreasemeraro/MPM_Space_Sciences/blob/main/Celestial%20mechanics%20and%20astrodynamics/Orbit.py) defines the class that                evaluates orbit parameters (converting from cartesian to keplerian or viceversa); 
      - [congiunzione_orbite.py](https://github.com/andreasemeraro/MPM_Space_Sciences/blob/main/Celestial%20mechanics%20and%20astrodynamics/congiunzione_orbite.py)           is able to propagate the equations of motions for two objects, in a given time and a given number of times along the orbit. Given a           threshold,             the code computes the points where the two objects are sufficiently close (in 2ver, it determines only the exact point at which the objects come                     close enough, warns for a conjunctions under the threshold and breaks the loop). There are also plots, and two ways to interpolate space points, in order to        determine distance and time of closest encounter for a limited set of points chosen by the user. The code is a preliminary study, so it is not supposed to           work without a bit of manipulating.
   - [Constellation](https://github.com/andreasemeraro/MPM_Space_Sciences/tree/main/Costellation):
   After an intial definition of the ranges of the constelltion          parameters N,P,h,i,F; the algorithm calculates the objective functions J1, J3, J5 and J7 s 4-D matrices. Furthermore the global function J is calculated as the      mean of all objective functions. Finally the algorithm finds the minimum value of J and gives the optimal parameters of the constellation.
   
   
   - [Data_Analysis](https://github.com/andreasemeraro/MPM_Space_Sciences/tree/main/Data_Analysis): 
      - [BONINO](https://github.com/andreasemeraro/MPM_Space_Sciences/tree/main/Data_Analysis/Bonino): File [Esame1.R](https://github.com/andreasemeraro/MPM_Space_Sciences/blob/main/Data_Analysis/Bonino/Esame1.R) performs a descriptive analysis of  flux1 of the Sun observation data set. It works on flux1, flux1diff (called ''s1diff''), checks the residual features. File [regression_esame1.R](https://github.com/andreasemeraro/MPM_Space_Sciences/blob/main/Data_Analysis/Bonino/regression_esame1.R) performs a scatter plot of the two fluxes to verify an eventual relationship               amog them.


       - [Fenu_Maldera](https://github.com/andreasemeraro/MPM_Space_Sciences/tree/main/Data_Analysis/Fenu_Maldera): In [Esercizio_1.ipynb](https://github.com/andreasemeraro/MPM_Space_Sciences/blob/main/Data_Analysis/Fenu_Maldera/Esercizio_1.ipynb) the effective area of a detector while in [Esercizio_2.ipynb](https://github.com/andreasemeraro/MPM_Space_Sciences/blob/main/Data_Analysis/Fenu_Maldera/Esercizio_2.ipynb) the forward folding procedures for a dataset is performed.



   - [Launcher and propulsion system](https://github.com/andreasemeraro/MPM_Space_Sciences/tree/main/Launcher): 
   File [Lanciatore.xlsx](https://github.com/andreasemeraro/MPM_Space_Sciences/blob/main/Launcher/Lanciatore.xlsx) performs the stage              optimization of the launcher through Newton Raphson method. File [Propulsion system_.xlsx](https://github.com/andreasemeraro/MPM_Space_Sciences/blob/main/Launcher/Lanciatore.xlsx) calculates the total amount of propellant necessary for the trip from LEO to        LMO, using Tsiolkovsky???s rocket equation.
   
   
   - [Machine_Learning](https://github.com/andreasemeraro/MPM_Space_Sciences/tree/main/Machine_Learning): 
   The algorithm uses a given dataset of images of Earth's        surface and train a CNN to recognize if the image is a forest, cloudy, sea or a desert. We used different CNN and epochs in this algorithm. 
   
   
   - [Orbit perturbation](https://github.com/andreasemeraro/MPM_Space_Sciences/tree/main/Orbit%20perturbation):
      - matlab file [Ex_10_6_J2_final.m](https://github.com/andreasemeraro/MPM_Space_Sciences/blob/main/Orbit%20perturbation/Ex_10_6_J2_final.m) computes the            effects of J2 mars gravitational field and provides an average value of the norm of a_J2 (calculated during 700 days on Mars orbit), that is the perturbing          acceleration of the due to J2 effects.
      - File [DeltaV_Drag.m](https://github.com/andreasemeraro/MPM_Space_Sciences/blob/main/Orbit%20perturbation/DeltaV_Drag.m) calculates the total DeltaV that the drag perturbation determines, drag perturbation calculated by file [Example_Mars_Drag.m](https://github.com/andreasemeraro/MPM_Space_Sciences/blob/main/Orbit%20perturbation/Example_Mars_Drag.m), that        exploits a function of Mars atmosphere.
