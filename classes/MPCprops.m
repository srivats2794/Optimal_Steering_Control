classdef MPCprops
   
    properties
        Np, Nc % Prediction horizon and Control horizon
        st_ul, st_ll % State upper and lower limits
        in_ul, in_ll % Input upper and lower limits
        Q, R % State deviation and input deviation penality matrices
        mpc_choice
        kin_choice
    end

end

