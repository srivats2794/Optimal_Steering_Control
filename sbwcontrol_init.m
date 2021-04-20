%% Reference generator
reference_generator;

%% LQR init
veh.linmodchoice=2;
veh.bm= LinearStateSpace(veh,data); %(bicycle model)
lqr_c= LQRprops;
%% TUNE THESE FOR LQR
lqr_c.Q=diag([0.5,0.75,0.5,0.75]);
lqr_c.R= 0.25 ;
%%
lqr_init(lqr_c,veh);

%% Steer-By-Wire compensator
Ksbw= -0.0001; %Tune this for compensator performance

delta_sbw(1)=0; theta_vf=0;

%% Plotter choice
plot_states=1; % Do you want to plot vehicle responses?