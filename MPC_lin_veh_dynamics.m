% Uses CasaDi to solve for the multiple shooting MPC problem
clc;clear;

run('C:\Users\srivats\Desktop\vehicle_control\reference_generator.m');

mpc_c= MPCprops;

addpath('C:\Users\srivats\Documents\MATLAB\casadi-v3.5.1')
import casadi.*

veh.linmodchoice=1;
veh.bm= LinearStateSpace(veh,data); %(bicycle model)
mpc_c.Np=40; % Prediction horizon
mpc_c.Nc=1; % Control horizon

mpc_c.Q = diag([0.25,0.75,0.5,0.75]); % State deviation penalty
mpc_c.R = 0.5; % Input deviation penalty
mpc_c.st_ul= max(ref_states(:,2:11))+0.001; % Upper bound for states
mpc_c.st_ll= min(ref_states(:,2:11))-0.001; % Lower bound for states
mpc_c.in_ul= deg2rad(5); % Upper bound for input
mpc_c.in_ll= -mpc_c.in_ul; % Lower bound for input

% State symbols: Lateral displacement, lateral vel, yaw, yawrate
% Control symbols: Steering angle - delta
y = SX.sym('y'); ydot = SX.sym('ydot'); psi = SX.sym('psi'); psidot = SX.sym('psidot'); 
delta = SX.sym('delta');
states = [y;ydot;psi;psidot]; n_states = length(states);
controls = delta; n_controls = length(controls);

%% f(x,u)
rhs = (veh.bm.Ac*[y;ydot;psi;psidot])+(veh.bm.Bc*delta); % Xdot= Ax+Bu
f = Function('f',{states,controls},{rhs}); 

%% MPC setup
U = SX.sym('U',n_controls,mpc_c.Np); 

P = SX.sym('P',n_states + mpc_c.Np*(n_states+n_controls));
% parameters (which include the initial state and the reference along the
% predicted trajectory (reference states and reference controls))

X = SX.sym('X',n_states,(mpc_c.Np+1));
% A vector that represents the states over the optimization problem.

%% Objectives over horizon

obj = 0; % Objective function
g = [];  % constraints vector


st  = X(:,1); % initial state
g = [g;st-P(1:4)]; % initial condition constraints
for k = 1:mpc_c.Np
    st = X(:,k);  con = U(:,k);
    obj = obj+(st-P(5*k:5*k+3))'*mpc_c.Q*(st-P(5*k:5*k+3)) + ...
              (con-P(5*k+4))'*mpc_c.R*(con-P(5*k+4)) ; % calculate obj
    % the number 5 is (n_states+n_controls)
    rot= [cos(st(3)) -sin(st(3));
          sin(st(3))  cos(st(3))];
    vxy_loc= rot'*[data.Vx_des; st(2)];
    st(2)= vxy_loc(2);
    st_next = X(:,k+1);
    k1 = f(st, con); 
    k2 = f(st + data.Ts/2*k1, con);
    k3 = f(st + data.Ts/2*k2, con); 
    k4 = f(st + data.Ts*k3, con); 
    
    k1(1)= data.Vx_des*sin(st(3))+k1(1)*cos(st(3));  
    k2(1)= data.Vx_des*sin(st(3))+k2(1)*cos(st(3));  
    k3(1)= data.Vx_des*sin(st(3))+k3(1)*cos(st(3));  
    k4(1)= data.Vx_des*sin(st(3))+k4(1)*cos(st(3));  
    
    st_next_prop=st +data.Ts/6*(k1 +2*k2 +2*k3 +k4); % new
         
    st_next_prop(2)= data.Vx_des*sin(st_next_prop(3))+st_next_prop(2)*cos(st_next_prop(3));  
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

args.lbx(1:4:4*(mpc_c.Np+1),1) = mpc_c.st_ll(2); %state y lower bound 
args.ubx(1:4:4*(mpc_c.Np+1),1) = mpc_c.st_ul(2); %state y upper bound  
args.lbx(2:4:4*(mpc_c.Np+1),1) = mpc_c.st_ll(5); %state ydot lower bound
args.ubx(2:4:4*(mpc_c.Np+1),1) = mpc_c.st_ul(5); %state ydot upper bound
args.lbx(3:4:4*(mpc_c.Np+1),1) = mpc_c.st_ll(3); %state psi lower bound
args.ubx(3:4:4*(mpc_c.Np+1),1) = mpc_c.st_ul(3); %state psi upper bound
args.lbx(4:4:4*(mpc_c.Np+1),1) = mpc_c.st_ll(6); %state psidot lower bound 
args.ubx(4:4:4*(mpc_c.Np+1),1) = mpc_c.st_ul(6); %state psidot upper bound 

args.lbx(4*(mpc_c.Np+1)+1:4*(mpc_c.Np+1)+mpc_c.Np,1) = mpc_c.in_ll; % delta lower bound
args.ubx(4*(mpc_c.Np+1)+1:4*(mpc_c.Np+1)+mpc_c.Np,1) = mpc_c.in_ul; % delta upper bound
%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP

% THE SIMULATION LOOP STARTS FROM HERE
%-------------------------------------------
states_curr = [data.inits.Px0;data.inits.Py0; data.inits.psi0; data.inits.Vx0;data.inits.Vy0; ...
               data.inits.psidot0; data.inits.theta0; data.inits.thetadot0; data.inits.phi0; ...
               data.inits.phidot0; data.inits.omega0; data.inits.omega0; data.inits.omega0; ...
               data.inits.omega0; data.inits.alpha0; data.inits.alpha0; data.inits.alpha0; data.inits.alpha0];

x0 = [states_curr(2) ; states_curr(5) ;states_curr(3) ; states_curr(6)];    % initial condition.

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
        args.p(1:4) = x0; % initial condition of the vehicle posture
        for k = 1:mpc_c.Np %set the reference to track
            y_ref = ref_states(mpciter+k+data.ld,3); 
            ydot_ref = ref_states(mpciter+k+data.ld,6);
            psi_ref = ref_states(mpciter+k+data.ld,4); 
            psidot_ref = ref_states(mpciter+k+data.ld,7);
            delta_ref= u_ref;
                      
            args.p(5*k:5*k+3) = [y_ref, ydot_ref, psi_ref, psidot_ref];
            args.p(5*k+4) = delta_ref;
        end
        
        %----------------------------------------------------------------------
        % initial value of the optimization variables
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
        x0 = [states_curr(2) ; vxy_plant(2) ;states_curr(3) ; states_curr(6)];  
        states_history(:,mpciter+2)= states_curr;
        states_history(4:5,mpciter+2)= vxy_plant;
        mpciter
        mpciter = mpciter + 1;
end

main_loop_time = toc(main_loop);
average_mpc_time = (main_loop_time/(mpciter+1))
beta= atan2(states_history(5,:),states_history(4,:))-states_history(3,:);

states_history(19,:)= beta;
clearvars -except states_history ref_states data veh plant_loop_time ...
           x0_curr average_mpc_time average_plant_time

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