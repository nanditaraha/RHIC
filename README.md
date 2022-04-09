# STAR at RHIC
In this experiment we collide Au+Au nuclei at 200 GeV energy in the cm reference at the STAR detector. 

In this code our goal is to simulate the collective behaviour of large ensembles of particles. 
It is known that every individual particle behaves uniquely (so fractals behaviours are unpredictable) but the 
collective behaviour of the entire system (or ensemble) is quite predictable may obey a certain pattern – this applies to 
all large systems like predicting the weather, stock market, anomalies in diseases, human behaviour etc. and can have 
applications to machine learning, VRML, data science etc. 

This code takes the following parameters as input:
<ul>
  <li> total multiplicity (n)</li>
   <li> no. of pairs exhibiting femtoscopy (a) </li>
    <li> no. of pairs with elliptical flow (b) </li> </ul>
        as arguments 
and simulates an ideal coleration function. Ideally the pairs generated from femtoscopy are very close to each other and 
so the relative rapidity and relative azimuthal angle is nearly zero - resulting a peak near zero in the correlation function R<sub>2</sub>.
The various patterns of the anisotropic flow of ultra-RHIC can be characterized by the Fourier expansion of the invariant 
triple distribution of particle pairs. See eqn. 2 and 3 of the reference below:
https://arxiv.org/abs/1102.3010 
The conservation of momentum (represented by the first term i.e. v1 in the Fourier expansion and is sinusoidal in &phi; - &psi; i.e.cos(&phi; - &psi;)
where &phi; is the azimuthal angle and &psi; is the reaction plane angle.
Similarly the eliptical flow v<sub>2</sub> due to the Lorentz boost is quantified by v<sub>2</sub> i.e. cos2(&phi; - &psi;) which is evident in the ridge of 
the correlation function R<sub>2</sub>.</br>
NOTE: The total number of particles N < 500. (N = S +2*a +2*b) where S is the # of single particles.
Code Usage:
make clean
make
./m -n 40 -a 10 -b 1
Output: Running like the command above creates a pdf called 'm_40_10_1_0.100.pdf' in the ps directory of the current working directory
