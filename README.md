# Optimal_Steering_Control
A look into performance of LQR, LQR+SbW compensation and Linear MPC

Add 6DoF_plant_functions, classes and init_files folders to your Matlab Path

Run SbWAdaptiveControl file for simulation.

The timestep, road condition and operating/desired velocity can all be manipulated using the init_files/sim_params.m .

The classes are vehicle,simprops and LQRprops and MPCprops

vehicle class carries vehicle params and also contains functions that can initialize a desired linear model. Model 1 (using linmodchoice variable): y ydot psi psidot statespace. Model 2: ey eydot epsi epsidot statespace also known as error dynamics state space. Model 3: psidot and beta statespace also known as sideslip model.

The 6DoF_plant_functions folder contains all the files required to simulate the high DoF plant that the controller is tested against. It includes 6DoF chassis and pacjeka wheel-tire model. You can ask for next state, velocity states/states_dot and Forces as output.

The MPC is solved using CasaDi for Matlab: https://web.casadi.org/

The steer-by-wire compensator was designed based on the concept by Dr. Rajamani in his Vehicle Dynamics and Control 2nd edition book: R. Rajamani, Vehicle dynamics and control. New York: Springer, 2012.

Each code file is commented.
