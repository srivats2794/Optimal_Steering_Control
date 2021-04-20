clc;clear;

run('C:\Users\srivats\Desktop\vehicle_control\reference_generator.m');
mpc_c= MPCprops;

addpath('C:\Users\srivats\Documents\MATLAB\casadi-v3.5.1')
import casadi.*

veh.linmodchoice=2;
veh.bm= LinearStateSpace(veh,data); %(bicycle model)
mpc_c.Np=40; 
mpc_c.Nc=1;

mpc_c.Q = diag([1,1,1,1]);
mpc_c.R = 1;
mpc_c.st_ul= [0.01;0.01;1;1];
mpc_c.st_ll= -mpc_c.st_ul;
mpc_c.in_ul= deg2rad(10);
mpc_c.in_ll= -mpc_c.in_ul;

e_y = SX.sym('e_y'); e_ydot = SX.sym('e_ydot'); ...
e_psi = SX.sym('e_psi'); e_psidot = SX.sym('e_psidot'); 
delta = SX.sym('delta');
psidot_des = SX.sym('psidot_des');

states = [e_y;e_ydot;e_psi;e_psidot]; n_states = length(states);
controls = delta; n_controls = length(controls);
pars=psidot_des;

%% f(x,u)
rhs = veh.bm.Ad*states+veh.bm.Bd*controls+veh.bm.Kd_psidotdes*pars; 
f = Function('f',{states,controls,pars},{rhs}); 

%% MPC setup
U = SX.sym('U',n_controls,mpc_c.Np); 
P = SX.sym('P',n_states+mpc_c.Np);
X = SX.sym('X',n_states,(mpc_c.Np+1));
% A vector that represents the states over the optimization problem.

%% Objectives over horizon

obj = 0; % Objective function
g = [];  % constraints vector

st  = X(:,1); % initial state
g = [g;st-P(1:4)]; % initial condition constraints
for k = 1:mpc_c.Np
    u_ref=0;
    st = X(:,k);  con = U(:,k); par= P(k+4);
    obj = obj+(st)'*mpc_c.Q*(st) + (con-u_ref)'*mpc_c.R*(con-u_ref); % calculate obj
    st_next_prop= f(st,con,par);
    st_next=X(:,k+1);
    g = [g;st_next-st_next_prop]; % compute constraints
end
% make the decision variable one column  vector
OPT_variables = [reshape(X,4*(mpc_c.Np+1),1);reshape(U,mpc_c.Np,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;

args.lbg(1:4*(mpc_c.Np+1)) = 0;  % -1e-20  % Equality constraints
args.ubg(1:4*(mpc_c.Np+1)) = 0;  % 1e-20   % Equality constraints

args.lbx(1:4:4*(mpc_c.Np+1),1) = mpc_c.st_ll(1); %error y lower bound 
args.ubx(1:4:4*(mpc_c.Np+1),1) = mpc_c.st_ul(1); %error y upper bound  
args.lbx(2:4:4*(mpc_c.Np+1),1) = mpc_c.st_ll(2); %error ydot lower bound
args.ubx(2:4:4*(mpc_c.Np+1),1) = mpc_c.st_ul(2); %error ydot upper bound
args.lbx(3:4:4*(mpc_c.Np+1),1) = mpc_c.st_ll(3); %error psi lower bound
args.ubx(3:4:4*(mpc_c.Np+1),1) = mpc_c.st_ul(3); %error psi upper bound
args.lbx(4:4:4*(mpc_c.Np+1),1) = mpc_c.st_ll(4); %error psidot lower bound 
args.ubx(4:4:4*(mpc_c.Np+1),1) = mpc_c.st_ul(4); %error psidot upper bound 

args.lbx(4*(mpc_c.Np+1)+1:4*(mpc_c.Np+1)+mpc_c.Np,1) = mpc_c.in_ll; % delta lower bound
args.ubx(4*(mpc_c.Np+1)+1:4*(mpc_c.Np+1)+mpc_c.Np,1) = mpc_c.in_ul; % delta upper bound

%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP

% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
states_curr = [data.inits.Px0;data.inits.Py0; data.inits.psi0; data.inits.Vx0;data.inits.Vy0; ...
               data.inits.psidot0; data.inits.theta0; data.inits.thetadot0; data.inits.phi0; ...
               data.inits.phidot0; data.inits.omega0; data.inits.omega0; data.inits.omega0; ...
               data.inits.omega0; data.inits.alpha0; data.inits.alpha0; data.inits.alpha0; data.inits.alpha0];

x0 = [0 ; 0 ;0 ; 0];    % initial condition.

u0 = zeros(mpc_c.Np,1);        % two control inputs for each robot
X0 = repmat(x0,1,mpc_c.Np+1)'; % initialization of the states decision variables

sim_tim = data.Tsim+7; % Maximum simulation time

data.Tvec= 0:data.Ts:(sim_tim); % Time vector
states_history= zeros(18,data.N+1);
states_history(:,1)= states_curr;
% Start MPC
mpciter = 0;
u_cl=[];

[ref_size,~]= size(ref_states);

u_ref=0;
main_loop = tic;
while(mpciter < sim_tim / data.Ts) % new - condition for ending the loop
    remainder= rem(mpciter,mpc_c.Nc);
    if remainder==0 || mpciter==0  
        args.p(1:4) = x0; % initial condition of the robot posture
        for k = 1:mpc_c.Np %new - set the reference to track
            psidot_ref = ref_states(mpciter+k+data.ld,7);
            args.p(k+4) = psidot_ref;
        end
        args.x0  = [reshape(X0',4*(mpc_c.Np+1),1);reshape(u0',mpc_c.Np,1)];
        sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
        u = reshape(full(sol.x(4*(mpc_c.Np+1)+1:end))',1,mpc_c.Np)'; % get controls only from the solution
        X0 = reshape(full(sol.x(1:4*(mpc_c.Np+1)))',4,mpc_c.Np+1)'; % get solution TRAJECTORY
        u_curr= u(1:mpc_c.Nc,:);
        X0 = [X0(mpc_c.Nc+1:end,:);repmat(X0(end,:),mpc_c.Nc,1)];
        u0 = [u(mpc_c.Nc+1:size(u,1),:);u(size(u,1),:)*ones(mpc_c.Nc,1)];
        count=1;
    end
    data.inputs.delta(mpciter+1)= u_curr(count,:);
    u_ref= u_curr(count,:);
    ins= [0;0;u_curr(count,:)];
    states_curr = plant_euler_forward(data.Ts,states_curr,ins,veh.params);
    rotP1= [cos(states_curr(3)) -sin(states_curr(3));
            sin(states_curr(3))  cos(states_curr(3))];
    vxy_plant=rotP1*[states_curr(4);states_curr(5)];
    y_err= states_curr(2)-ref_states(mpciter+1+data.ld,3);
    ydot_err= vxy_plant(2)-ref_states(mpciter+1+data.ld,6);
    psi_err= states_curr(3)-ref_states(mpciter+1+data.ld,4);
    psidot_err= states_curr(6)-ref_states(mpciter+1+data.ld,7);
    x0 = [y_err ; ydot_err ; psi_err; psidot_err];
    states_history(:,mpciter+2)= states_curr;
    states_history(4:5,mpciter+2)= vxy_plant;
    count=count+1;
    mpciter
    mpciter = mpciter + 1;
end

main_loop_time = toc(main_loop);
average_mpc_time = main_loop_time/(mpciter+1)
beta= atan2(states_history(5,:),states_history(4,:))-states_history(3,:);

states_history(19,:)= beta;
clearvars -except states_history ref_states data veh plant_loop_time

[~,columns]=size(states_history);
Tvec_sim= 0:data.Ts:((columns-1)*data.Ts); % Time vector

figure(1)
nexttile
plot(states_history(1,:),states_history(2,:),'-k',ref_states(:,2),ref_states(:,3),'--b','LineWidth',1);
set(gca,'fontsize',12);
xlim([-5 400]) 
ylim([-1 4])
title ('Plant Trajectory vs Reference','FontSize',10);
legend('Plant','Reference','FontSize',10);
ylabel('Lateral Displacement(m)','FontSize',15); xlabel('Longitudinal Displacement(m)','FontSize',15);

nexttile
plot(Tvec_sim,rad2deg(states_history(3,:)),'-k',Tvec_sim,rad2deg(ref_states(1:columns,4)),'--b','LineWidth',1);
set(gca,'fontsize',12);
title ('Yaw of Plant','FontSize',10);
ylabel('$\psi$ (deg)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
legend('Plant','Reference','FontSize',10);

nexttile
plot(Tvec_sim,rad2deg(states_history(6,:)),'-k',Tvec_sim,rad2deg(ref_states(1:columns,7)),'--b','LineWidth',1);
set(gca,'fontsize',12);
title ('Yaw Rate of Plant','FontSize',10);
ylabel('$\dot{\psi}$ (deg/s)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
legend('Plant','Reference','FontSize',10);

nexttile
plot(Tvec_sim,states_history(4,:),'-k',Tvec_sim,ref_states(1:columns,5),'--b','LineWidth',1);
set(gca,'fontsize',12);
title ('Longitudinal Velocity of Plant','FontSize',10);
ylabel('$V_{x}$ (m/s)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
legend('Plant','Reference','FontSize',10);

nexttile
plot(Tvec_sim,states_history(5,:),'-k',Tvec_sim,ref_states(1:columns,6),'--b','LineWidth',1);
set(gca,'fontsize',12);
title ('Lateral Velocity of Plant','FontSize',10);
ylabel('$V_{y}$ (m/s)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
legend('Plant','Reference','FontSize',10);

nexttile
plot(Tvec_sim,rad2deg(states_history(19,:)),'-k','LineWidth',1);
set(gca,'fontsize',12);
title ('Sideslip of Plant','FontSize',10);
ylabel('$\beta$ (deg)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);

sgtitle('VEHICLE RESPONSE PLOTS','FontSize',20)

figure(2)
plot(Tvec_sim(:,1:(end-1)),rad2deg(data.inputs.delta),'-r','LineWidth',1);
set(gca,'fontsize',12);
title ('Inputs by MPC','FontSize',10);
ylabel('$\delta$ (deg)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);