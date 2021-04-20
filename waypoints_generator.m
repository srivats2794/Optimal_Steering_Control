clc;clear;
init;

% if you have continuous time reference for Px,Py,Psi,Vx,Vy,PsiDot already skip this and
% reference_generator files. 
load('waypoints_40cm.mat'); % CarSim recorded data for a DLC maneuver. 
% Replace with your own waypoints

clearvars -except X_cg Y_cg data veh
divisor= 0.75/0.0005;
j=1;
[rows,~]= size(X_cg.Data(:,1));
for i= 1:rows
    reminder= rem(i,divisor);
    if i==1 || reminder==0
     x_ref(j,1)=X_cg.Data(i,1); 
     y_ref(j,1)=Y_cg.Data(i,1);
     j=j+1;
    end      
end

clearvars -except x_ref y_ref data veh 
% scatter(x_ref,y_ref)
[rows,~]=size(x_ref);
t_ref(1,1)=0;

% Calculating desired time of arrival for each waypoint. 
for i=2:rows
    dist= sqrt(((x_ref(i)-x_ref(i-1))^2)+((y_ref(i)-y_ref(i-1))^2));
    t_ref(i,1)= t_ref(i-1)+dist/data.Vx_des;
end

trajectory = waypointTrajectory([x_ref, y_ref, zeros(rows,1)], ...
    'TimeOfArrival',t_ref, ...
    'SampleRate',1/data.Ts);
clearvars -except data veh trajectory 