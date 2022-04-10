# Simulations of correlation functions at STAR.
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
and simulates an ideal colerration function. </br>

## Femtoscopy
Ideally the pairs generated from femtoscopy are very close to each other and 
so the relative rapidity and relative azimuthal angle is nearly zero - resulting a peak near zero in the correlation function R<sub>2</sub>.
In an example from the code with a = 10 and b = 1, this effect of femtoscopy is prominant as shown in the peak of the correlation function R<sub>2</sub>(dy,d&phi;) below:</br>
<img width="332" alt="Screen Shot 2022-04-09 at 1 41 22 PM" src="https://user-images.githubusercontent.com/27436642/162585361-8b4bea4c-ff44-4ce1-8c8a-a1a04aeefed1.png"><img width="329" alt="Screen Shot 2022-04-09 at 1 42 14 PM" src="https://user-images.githubusercontent.com/27436642/162585394-bbd82d52-b0ca-41d2-9469-d0a20ebf8116.png">

## Elliptical flow
The various patterns of the anisotropic flow of ultra-RHIC can be characterized by the Fourier expansion of the invariant 
triple distribution of particle pairs. See eqn. 2 and 3 of the reference below:</br>
https://arxiv.org/abs/1102.3010 <br>
<img width="389" alt="Screen Shot 2022-04-09 at 1 19 36 PM" src="https://user-images.githubusercontent.com/27436642/162584653-41f07d03-7dc6-4d85-837f-007e8c9f3e05.png">
</br>The conservation of momentum (represented by the first term i.e. v<sub>1</sub> in the Fourier expansion and is sinusoidal in &phi; - &psi; i.e.cos(&phi; - &psi;)
where &phi; is the azimuthal angle and &psi; is the reaction plane angle.
Similarly the eliptical flow v<sub>2</sub> due to the Lorentz boost is quantified by v<sub>2</sub> i.e. cos2(&phi; - &psi;) which is evident in the ridge of the correlation function R<sub>2</sub>. In an example from the code with a = 0 and b = 20, only the effect of elliptical flow persists as shown in the cos2(&phi; - &psi;) variation of R<sub>2</sub>(d&phi;) - right panel and a ridge in R<sub>2</sub>(dy,d&phi;) - left panel below:</br>
<img width="697" alt="Screen Shot 2022-04-09 at 2 04 34 PM" src="https://user-images.githubusercontent.com/27436642/162586235-c5d3cb90-ae9d-4de5-abcc-c73070b558e9.png"></br>
NOTE: The total number of particles N < 500. (N = S +2*a +2*b) where S is the # of single particles.</br>

### Code Usage:</br>

#### Must have at C++ and ROOT 6.24/06  installed.
All code is saved in Multshape. Create a folder called ps in your working directory (if it is not there) and follow the simple instructions below</br>
> mkdir ps</br>
> make clean</br>
> make</br>
> ./m -n 40 -a 10 -b 1</br>
Output: 
Running like the command above creates a pdf called 'm_40_10_1_0.100.pdf' in the ps directory of the current working directory
