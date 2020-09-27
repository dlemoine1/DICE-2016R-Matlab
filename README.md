# DICE-2016R-Matlab

Created by Derek Lemoine, University of Arizona

September 2020

Suggested citation: Lemoine, Derek.  2020.  DICE-2016R-Matlab. https://github.com/dlemoine1/DICE-2016R-Matlab

These files replicate DICE-2016R in Matlab.  They also extend it to use carbon and climate models recommended by Dietz et al. (2020) and to use a damage specification that builds on the expert survey in Pindyck (2019), as implemented by Lemoine (2021).

See Appendix to "Incentivizing Negative Emissions Through Carbon Shares" (Lemoine, 2020) for model equations and details: https://eller.arizona.edu/sites/default/files/Econ-WP-20-08.pdf.  See Nordhaus (2017) for original DICE-2016R description.

All files should be placed in the same directory.  Code will create an "output" subfolder within it.

main_dice2016r.m: The code should be run from here.  Options to control the climate, carbon, and damage models are provided at the top.  The code has not necessarily been tested on the non-default computational options.  In the default approach, solver guesses trajectories for all state variables and all non-consumption controls, imposes transition equations as constraints, and uses an analytic gradient.  Outputs the user may be interested in include Welfare, abaterate, emtax_pertCO2 (in 2010 dollars), T (temperature, deg C wrt 1900), Carbon_ppm (atmospheric carbon, in ppm), and emsind (industrial emissions, Gt C).

sub_parameters.m: Defines equations of the model and parameterizes the model.  Some users may want to change some of these parameters.  When changing equations, make sure to also change their derivatives.

sub_loadguesses.m: Loads an initial guess, defines upper and lower bounds, defines normalization for computation

utilityobjective.m: Policymaker's objective, with gradient

nonlcon_utilmax.m: The constraints of the model, with gradients

trajectory.m: Simulates the model with a given trajectory of controls

*.opt: Options files for Knitro solver.  Not needed if will use Matlab's built-in fmincon instead.

DICE2016R-091916ap.gms: Source file for DICE-2016R, from http://www.econ.yale.edu/~nordhaus/homepage/homepage/DICE2016R-091916ap.gms.  Provided only for reference.  Not needed to run code.



Reference List:

Dietz, Simon, Frederick van der Ploeg, Armon Rezai, and Frank Venmans. 2020. “Are Economists Getting Climate Dynamics Right and Does It Matter?” Working Paper 8122. CESIfo. https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3545718.

Lemoine, Derek. 2021. “The Climate Risk Premium: How Uncertainty Affects the Social Cost of Carbon.” Journal of the Association of Environmental and Resource Economists, no. forthcoming. http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2560031.

Nordhaus, William D. 2017. “Revisiting the Social Cost of Carbon.” Proceedings of the National Academy of Sciences 114 (7): 1518–23. https://doi.org/10.1073/pnas.1609244114.

Pindyck, Robert S. 2019. “The Social Cost of Carbon Revisited.” Journal of Environmental Economics and Management 94 (March): 140–60. https://doi.org/10.1016/j.jeem.2019.02.003.

