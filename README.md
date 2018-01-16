# shallow-waters

PURPOSE: Replicate the results obtained in:

The Seasonal Upwelling in the Gulf of Guinea Due to Remote Forcing
Journal of Physical Oceanography, Vol. 8, 1050-1060
1978
by David Adamec and James J. O'Brien

for the course Numerical Modeling of the Atmosphere of the bachelor in
Atmospheric Sciences - University of the Republic (Uruguay).

It consists of a linear model on an equatorial beta plane integrated over
a 100-days period in a basin that approximates the tropical Atlantic
Ocean.

The motion is described by the shallow water equations, discretized in a
staggered grid in space (type C). For the time derivative a leapfrog 
scheme (leapfrog.m) is considered. Every 30 iterations, the Matsuno 
scheme (matsuno.m) is used. 

Figures are made to show the results.

A stability analysis is made.
