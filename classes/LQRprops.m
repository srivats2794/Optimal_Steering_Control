classdef LQRprops < handle
    properties
       Q,R,K
    end
    
    methods
        function lqr_init(lqr_d, veh_d) 
            [~,lqr_d.K,~] = idare(veh_d.bm.Ad,veh_d.bm.Bd,lqr_d.Q,lqr_d.R,[],[]);
            lqr_d.K=-lqr_d.K;
        end
    end
end
    