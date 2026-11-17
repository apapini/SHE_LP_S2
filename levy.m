%Allows for the simulation of different Levy random fields all using
%the same decay rate (decay is implemented in the main code)

%kap = number of samples needed for each index
%PP = Poisson Process, BM = Brownian motion, CCP = Compound PP
%Option = select the type of process
%Lambda = intensity of the processes
%ell=current frequency
%stepsize= time step size for the increment of the Lévy processes

function LP = levy(Option, Lambda, ell, stepsize, kap)

    %number of random variable components
    numberRV =2*ell+1;
    
    if Option == 1 
        % generate just independent Brownian motions
        LP = sqrt(stepsize)*randn(numberRV,kap);
    elseif Option == 2 
        % generate just a independent compound Poisson processes
        LP = poissrnd(Lambda*stepsize,numberRV,kap) - Lambda*stepsize;
    elseif Option == 3 
        % generate sum of independent Brownian motions and compound Poisson
        % processes
        LP = sqrt(stepsize)*randn(numberRV,kap)+poissrnd(Lambda*stepsize,numberRV,kap)- Lambda*stepsize;
    elseif Option == 4 
        % generate identical compound Poisson processes
        LP = poissrnd(Lambda*stepsize,1,1)*ones(numberRV,kap) - Lambda*stepsize;
    end
end
