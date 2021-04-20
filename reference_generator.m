waypoints_generator;
count=1;
while ~isDone(trajectory)
   %[currentPosition,orientationLog(count)] = trajectory();
   %plot(currentPosition(1),currentPosition(2),'bo')
   [pos(count,:),orient(count),vel(count,:),acc(count,:),angVel(count,:)] = trajectory();
   time(count)= (count-1)*data.Ts;
   %pause(trajectory.SamplesPerFrame/trajectory.SampleRate)
   count = count + 1;
end

[row,~]= size(pos);

for i=1:row
   euler(:,i)= eulerd(orient(i),'ZYX','frame');
   yaw(i)= deg2rad(euler(1,i));
end

ref_states= [time',pos(:,1:2),yaw', vel(:,1:2),angVel(:,3)];
data.Tsim= time(end); % Simulation time
data.inits= states_init(data,veh,[ref_states(1,2);ref_states(1,3);ref_states(1,4);0;ref_states(1,6);ref_states(1,7)]);
%data.inits= states_init(data,veh,[0.1;1;0;0;0;0]);
data.N = single(data.Tsim/data.Ts);
data.Tvec= 0:data.Ts:(data.Tsim); % Time vector

clearvars -except ref_states data veh

ref_states(:,8:9)= [ref_states(:,5)/veh.Rw, ref_states(:,5)/veh.Rw];

ref_states(:,10:11)= 0;

difference= ref_states(end,2)-ref_states(end-1,2);
%Adding just copy of last known reference
%so that MPC doesnt throw an error at the last N value of loop 
% Cuz it'll look N+Np steps ahead. In an actual vehicle, reference will
% keep going or vel control will bring vehicle to stop.
for i= data.N+2:data.N+1500 
    ref_states(i,2)= ref_states(i-1,2)+difference;
    ref_states(i,3:5)= ref_states(i-1,3:5);
    ref_states(i,8:9)= ref_states(i-1,8:9);
end
ref_states(:,12)= atan((ref_states(:,7)*(veh.l))./ref_states(:,5));
clearvars difference i 
%% Reference form
% Time,Px,Py,Psi,Vx,Vy,PsiDot,OmegaF,OmegaR,AlphaF,AlphaR,Delta