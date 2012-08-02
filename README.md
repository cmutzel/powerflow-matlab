powerflow-matlab
================

This is my MATLAB power flow project for Power Systems Analysis class at Tufts.

The powerFlow script can be run on any size system ( not static 30) and requires input files similar to those found in this repo.

The power flow script uses a well-known newton raphson method to find the voltages and phase angles at each bus in a power grid.  An initial guess is made and the jacobian is updated after each iteration.  

check out the pdf for a better explanation.  