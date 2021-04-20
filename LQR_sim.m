% LQR works well until 60 KmpH
sbwcontrol_init; 

%% Vehicle response recorder - lqr (Rajamani cntrl off) sbw(Rajamani cntrl on) ref(linear plant)
states_history_lqr= zeros(18,data.N+1); % Pre-Allocation to speed up code.
states_history_sbw= zeros(18,data.N+1);

%% Initial state vectors

% 18 states of nonlinear plant -> Px, Py, Psi, Vx, Vy, PsiDot, Theta,
% ThetaDot, Phi, PhiDot, OmegaFR, OmegaFL, OmegaRR, OmegaRL, AlphaFR, 
% AlphaFL, AlphaRR, AlphaRL

states_plant_lqr= [data.inits.Px0;data.inits.Py0; data.inits.psi0; data.inits.Vx0;data.inits.Vy0; ...
         data.inits.psidot0; data.inits.theta0; data.inits.thetadot0; data.inits.phi0; ...
         data.inits.phidot0; data.inits.omega0; data.inits.omega0; data.inits.omega0; ...
         data.inits.omega0; data.inits.alpha0; data.inits.alpha0; data.inits.alpha0; data.inits.alpha0];
states_history_lqr(:,1)=states_plant_lqr;

% A vector that records velocities in global frame
% Cannot be used in dynamics loop. So have to make separate vector.
% Dynamics is in body coordinates
vxy_plant_lqr=[data.inits.Vx0;data.inits.Vy0]; 

for i= 1:data.N
    %% Errors computation for LQR control vector X for u=-KX (Error dynamics states)
    % LQR only plant
    e1_plant= states_plant_lqr(2)-ref_states(i+data.ld,3);
    e1dot_plant= vxy_plant_lqr(2)-ref_states(i+data.ld,6);
    e2_plant= states_plant_lqr(3)-ref_states(i+data.ld,4);
    e2dot_plant= states_plant_lqr(6)-ref_states(i+data.ld,7);
       
    x_lqr_plant = [e1_plant;e1dot_plant;e2_plant;e2dot_plant];
    
    %% Final inputs recorded ready to be used with dynamics propagation
    data.inputs.delta_lqr(i) = lqr_c.K*x_lqr_plant;     
    %% Nonlinear plant(s) forward propagation   
    % Inputs need to be of form front axle torque, rear axle torque,
    % steering angle (road-wheel)
    ins_lqr= [0;0;data.inputs.delta_lqr(i)];
    
    % Propagating dynamics. 
    % ins-> (step size,current states,inputs,parameters) 
    % outs-> next states , dot of states(optional), tire forces(optional)
    states_plant_lqr = plant_euler_forward(data.Ts,states_plant_lqr,ins_lqr,veh.params);
        
    % Rotation matrices to rotate velocities from dynamics propagation from
    % body to world
    rotP1= [cos(states_plant_lqr(3)) -sin(states_plant_lqr(3));
           sin(states_plant_lqr(3))  cos(states_plant_lqr(3))];
    
    vxy_plant_lqr=rotP1*[states_plant_lqr(4);states_plant_lqr(5)];
    
    %% Recording states in the history vector for plotting purposes
    states_history_lqr(:,i+1)= states_plant_lqr;
    states_history_lqr(4:5,i+1)= vxy_plant_lqr;
   
end

%% Calculating body sideslip beta (tan inv of Vy/Vx - psi)
beta_lqr= atan2(states_history_lqr(5,:),states_history_lqr(4,:))-states_history_lqr(3,:);

states_history_lqr(19,:)= beta_lqr;

clearvars -except states_history_lqr ref_states data veh ...
                  plot_states
[~,columns]=size(states_history_lqr);
Tvec_sim= 0:data.Ts:((columns-1)*data.Ts); % Time vector

if plot_states==1
    figure(1)
    nexttile
    plot(states_history_lqr(1,:),states_history_lqr(2,:),'-k', ...
        ref_states(1:columns,2),ref_states(1:columns,3),':b', ...
        'LineWidth',2);
    set(gca,'fontsize',12);
    title ('Trajectory','FontSize',10);
    legend('Vehicle Response','Reference','FontSize',7);
    ylabel('Lateral Displacement(m)','FontSize',15); xlabel('Longitudinal Displacement(m)','FontSize',15);
    
    nexttile
    plot(Tvec_sim,rad2deg(states_history_lqr(3,:)),'-k', ...
         Tvec_sim,rad2deg(ref_states(1:columns,4)),':b', ...
        'LineWidth',2);
    set(gca,'fontsize',12);
    title ('Yaw','FontSize',10);
    legend('Vehicle Response','Reference','FontSize',7,'Location','southwest');
    ylabel('$\psi$ (deg)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
    
    nexttile
    plot(Tvec_sim,rad2deg(states_history_lqr(6,:)),'-k', ...
         Tvec_sim,rad2deg(ref_states(1:columns,7)),':b', ...
         'LineWidth',2);
    set(gca,'fontsize',12);
    title ('Yaw Rate','FontSize',10);
    legend('Vehicle Response','Reference','FontSize',7,'Location','southwest');
    ylabel('$\dot{\psi}$ (deg)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
    
    nexttile
    plot(Tvec_sim,states_history_lqr(5,:),'-k', ...
         Tvec_sim,ref_states(1:columns,6),':b', ...
        'LineWidth',2);
    set(gca,'fontsize',12);
    title ('Lateral Velocity','FontSize',10);
    legend('Vehicle Response','Reference','FontSize',7,'Location','southwest');
    ylabel('$\dot{y}$ (m/s)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
    
    nexttile
    plot(Tvec_sim,states_history_lqr(4,:),'-k', ...
         Tvec_sim,ref_states(1:columns,5),':b', ...
         'LineWidth',2);
    set(gca,'fontsize',12);
    title ('Longitudinal Velocity','FontSize',10);
    legend('Vehicle Response','Reference','FontSize',7,'Location','southwest');
    ylabel('$\dot{x}$ (m/s)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
    
    nexttile
    plot(Tvec_sim,rad2deg(states_history_lqr(19,:)),'-k', ...
        'LineWidth',2);
    set(gca,'fontsize',12);
    title ('SideSlip','FontSize',10);
    ylabel('$\beta$ (deg)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
        
    sgtitle('VEHICLE RESPONSE PLOTS','FontSize',20)

    figure(2)
    plot(Tvec_sim(1:end-1),rad2deg(data.inputs.delta_lqr),'-k', ...
        'LineWidth',2);
    set(gca,'fontsize',12);
    title ('Control inputs','FontSize',10);
    ylabel('$\delta$ (deg)','FontSize',15,'interpreter','latex'); xlabel('Simulation time (sec)','FontSize',15);
    
end