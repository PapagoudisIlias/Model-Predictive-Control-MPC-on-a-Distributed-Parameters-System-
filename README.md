In this code I am using the MPC algorithm in order to control a distributed parameters system (DPS) with the CasADi solver. The DPS systems are systems that can be 
described with partial differential equations. In order to control this kind of systems, at fisrt I discretize the space (x-coordinate) with the use of 
the wavelets theory. In that way the PDE can be transformed at a system of ordinary ODEs. Then, I discretize the time with the Runge–Kutta 4 method. This is essential 
for the proper use of CasAdi. At last I use the MPC algorithm at the final system of differential algebraic equations.

You have to be sure that the CasaAdi is properly isntalled and that you have the correct address on the addpath command at the code.
You can download the CasAdi from https://web.casadi.org/
