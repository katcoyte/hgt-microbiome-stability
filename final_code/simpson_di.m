function [simps, invsimps] = simpson_di(data_vec)

    N = sum(data_vec);
    %f = data_vec/N;
    f = data_vec;
    fsq = f.^2;
    simps = sum(fsq);
    invsimps = 1/simps;
    
    %invsimps = -sum(f.*log(f));

end
