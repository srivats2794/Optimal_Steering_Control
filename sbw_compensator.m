function [delta_sbw,theta_vf_next]= sbw_compensator(vels,Fyf, delta_ref, theta_vf_curr,Ksbw,data, veh, deltasbw_curr)

lp= veh.Izz/(veh.m*veh.l_r);
ay_p= (veh.l/(veh.m*veh.l_r))*Fyf;

Neul= data.Ts/0.001;
    for i=1:Neul
        gfunc= (cos(theta_vf_curr)/vels(1))*(((veh.l_f-lp)*vels(4)* ...
                cos(theta_vf_curr))+((veh.l_f*vels(3)*vels(3)- ...
                vels(2))*sin(theta_vf_curr)));
        
        theta_vf_dot=  -vels(3)+((cos(theta_vf_curr)^2)/vels(1))*ay_p+gfunc;
       
        F_del_func= (Ksbw*vels(1)*delta_ref)/veh.l; %(vels(1)*delta_curr)/ ...
                    (veh.l+((Ksbw/g)*vels(1)*vels(1)));
        
        theta_vf_curr= theta_vf_dot*data.Ts+theta_vf_curr;  
    end
    deltadot= -vels(3)+gfunc+F_del_func;
    delta_sbw= deltadot*data.Ts+deltasbw_curr;
    theta_vf_next=theta_vf_curr;
end