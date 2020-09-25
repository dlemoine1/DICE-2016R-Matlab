# DICE-2016R-Matlab

Created by Derek Lemoine, University of Arizona

September 2020

These files replicate DICE-2016R in Matlab.  

They also extend it to use carbon and climate models recommended by Dietz et al. (2020) and to use a damage specification that builds on the expert survey from Pindyck (2019), as implemented by Lemoine (2021).
See Appendix to "Incentivizing Negative Emissions Through Carbon Shares" (Lemoine, 2020) for model equations and details: https://eller.arizona.edu/sites/default/files/Econ-WP-20-08.pdf

All files should be placed in the same directory.  Code will create an "output" subfolder within it.

main_dice2016r.m: The code should be run from here.  Options to control the climate, carbon, and damage models are provided at the top.  The code has not necessarily been tested on the non-default computational options.  In the default approach, solver guesses trajectories for all state variables and all non-consumption controls, imposes transition equations as constraints, and uses an analytic gradient.  Outputs the user may be interested in include Welfare, abaterate, emtax_pertCO2, T (temperature, deg C wrt 1900), and emsind (industrial emissions, Gt C).

sub_parameters.m: Defines equations of the model and parameterizes the model.  Some users may want to change some of these parameters.  When changing equations, make sure to also change their derivatives.

sub_loadguesses.m: Loads an initial guess, defines upper and lower bounds, defines normalization for computation

utilityobjective.m: Policymaker's objective, with gradient

nonlcon_utilmax.m: The constraints of the model, with gradients

trajectory.m: Simulates the model with a given trajectory of controls

*.opt: Options files for Knitro solver.  Not needed if will use Matlab's built-in fmincon instead.
