%Allows for the simmulation of different Levy random fields all using
%the same decay rate (all under Ass 2.8 with equality)
%The decay rate is later on implemented in the main part
%Here only 2(l+1) identically distributed Levy samples are generated
%kap=truncation parameter of the spectral methods
function LP = levy(Option, Lambda, ell, stepsize, kap)
    numberRV =2*ell+1;
    if Option == 1 
        % generate just a BM
        LP = sqrt(stepsize)*randn(numberRV,kap);
    elseif Option == 2 
        %generate just a PP
        LP = poissrnd(Lambda*stepsize,numberRV,kap);
    elseif Option == 3 
        %generate BM + PP
        LP = sqrt(stepsize)*randn(numberRV,kap)+poissrnd(Lambda*stepsize,numberRV,kap);
    elseif Option == 4 
        %generate just a same PP
        LP = poissrnd(Lambda*stepsize,1,1)*ones(numberRV,kap);
        %Does not really make sense. It just adds up the spherical
        %harmonics and obtains a wierd dependence structure. If sth like
        %this reasonable, then with CPP such that direction can change...
    elseif Option == 5 
        %generate BM + same PP
        PP = poissrnd(Lambda*stepsize,1,1)*ones(numberRV,kap);
        LP = sqrt(stepsize)*randn(numberRV,kap)+PP;
    elseif Option == 6 
        %generate BM + CCP with jumps of height -1 or 1 with prob 0.5
        PP = poissrnd(Lambda*stepsize,numberRV,kap);
        for k = 1:numberRV
            PP(k) = sum((-1+2*binornd(1,0.5,PP(k),1)));
        end
        LP = sqrt(stepsize)*randn(numberRV,kap)+PP;
    elseif Option == 7 
        %generate BM + CCP with jumps of height -4 or 4 with prob 0.5
        PP = poissrnd(Lambda*stepsize,numberRV,kap);
        for k = 1:numberRV
            PP(k) = sum((-4+8*binornd(1,0.5,PP(k),1)));
        end
        LP = sqrt(stepsize)*randn(numberRV,kap)+PP;
    elseif Option == 8 
        %generate BM + CCP with jumps normally distr
        PP = poissrnd(Lambda*stepsize,numberRV,kap);
        for k = 1:numberRV
            PP(k) = sum(randn(PP(k),1));
        end
        LP = sqrt(stepsize)*randn(numberRV,kap)+PP;
    elseif Option == 9 
        %generate BM + CCP with same PP and jumps of height -1 or 1 with prob 0.5
        PP = poissrnd(Lambda*stepsize,1,1)*ones(numberRV,kap);
        for k = 1:numberRV
            PP(k) = sum((-1+2*binornd(1,0.5,PP(k),1)));
        end
        LP = sqrt(stepsize)*randn(numberRV,kap)+PP;
    elseif Option == 10 
        %generate BM + CCP with same PP and jumps of height -4 or 4 with prob 0.5
        PP = poissrnd(Lambda*stepsize,1,1)*ones(numberRV,kap);
        for k = 1:numberRV
            PP(k) = sum((-4+8*binornd(1,0.5,PP(k),1)));
        end
        LP = sqrt(stepsize)*randn(numberRV,kap)+PP;
    elseif Option == 11 
        %generate BM + CCP with same PP and jumps normally distr
        PP = poissrnd(Lambda*stepsize,1,1)*ones(numberRV,kap);
        for k = 1:numberRV
            PP(k) = sum(randn(PP(k),kap));
        end
        LP = sqrt(stepsize)*randn(numberRV,kap)+PP;
    else 
        LP = zeros(numberRV,kap);
    end
end
