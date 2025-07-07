% Strong convergence stochastic heat equation 
% with additive Lévy noise on the sphere
% Monte Carlo simulation
clear all 
close all

%set seed
randn('state',112);

M = 10; % number of MC samples

% truncation index
k = 10;
kappa = 2^k;

% time interval
time = 1;

% initial condition, decay spectrum
beta = 3;

% noise regularity
alpha = [1 2 3 4 5];
ExpLP=1;
VarLP=1;

%vector of strong error
strong_error_MC_kappa = zeros(k,length(alpha));

%variable to check the error against lower order truncation.
j = k-1;
%auxiliary vector of error computad for every regularity parameter.
err = zeros(1,length(alpha));

%loop to constract the teoretical spectral error.
for monte = 1:M
    for i=0:(kappa-1)
        % initial condition = 0
        err = err + (1 - exp(-2*(kappa-i)*(kappa-i+1)*time))*(1+2*(kappa-i))/(2*(kappa-i)*(kappa-i+1))*(kappa-i).^(-alpha)*VarLP...
                + (1 - exp(-(kappa-i)*(kappa-i+1)*time))^2 *(1+2*(kappa-i))/((kappa-i)*(kappa-i+1))^2*(kappa-i).^(-alpha)*ExpLP^2;
        %loop to stop the troncation and place the error in the right position
        if (kappa-i)==(2^j +1)
            strong_error_MC_kappa(j+1,:) = sqrt(err);
            j=j-1;
        end
    end
end

figure();
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(1)/2),'--','Color','#0072BD');
hold on;
loglog(2.^(0:1:(k-1)), strong_error_MC_kappa(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(2)/2),':','Color','#77AC30');
loglog(2.^(0:1:(k-1)), strong_error_MC_kappa(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(3)/2),'-.','Color','#D95319');
loglog(2.^(0:1:(k-1)), strong_error_MC_kappa(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(4)/2),'-','Color','#1ec8c3');
loglog(2.^(0:1:(k-1)), strong_error_MC_kappa(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(5)/2),'.-','Color','#fb8f36');
loglog(2.^(0:1:(k-1)), strong_error_MC_kappa(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
hold off;

title('Strong error spectral approximation')
xlabel('Number of series elements $\kappa$','Interpreter','latex')
ylabel('Strong error')
h_legend=legend('$O(\kappa^{1/2})$','$\alpha = 1$','$O(\kappa)$','$\alpha = 2$','$O(\kappa^{3/2})$','$\alpha = 3$','$O(\kappa^{2})$','$\alpha = 4$','$O(\kappa^{5/2})$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
print -depsc2 -r0 strong_MC_spectral.eps
%-----------------------------------------------------------------------------------------------------------------------------

%% Forward EM, real MC simulation

M = 10; % number of MC samples

time = 1; %final time
n = 14; % time grid for reference solution with 2^n steps

% space parameter
alpha = [1 2 3 4 5]; % regularity of the component of lévy noise
kappa = floor(2^(n/2)); % truncation space

tic %tic-toc for the time of the simulation
% error stored in vector for different time disretizations
strong_error_MC_FEM = zeros(n-1,length(alpha));  

%auxiliary vector to store the error
error = zeros(n-1,length(alpha));
%loop to compute the solution error and the MC loop
for monte = 1:M
    for i = 1:kappa 
        
        %selection of the noise and intensity
        option =2;
        lambda = 1;
        %selection of the timestep
        h = time/2^n;
        LP = levy(option, lambda, i, h,2^n)';
        
        %Construction of the final time  reference solution and vector of solution
        %initial value 0, only noise
        final_time_sol_reference = zeros(length(alpha),1+2*i);
        sol = zeros(length(alpha),1+2*i); 
        
        %computation of the 1*ell+1 component
        for j = 1:(2^n) 
                sol = sol + (i+1).^(-alpha.'/2).*(1-i*(i+1)*h)^(2^n-j).*LP(j, :); %fEM as reference solution'i-s+1
        end;
        
        % rescale solutions for all but X_l0.
        final_time_sol_reference = sol;

        %computation of the component for the solution at different h and k
        for m = 1:(n-1)
            %reduction of the noise factor to the right size for the smaller solution 
            B = LP(1:(length(LP(:,1))-1),:) + LP(2:length(LP(:,1)),:);
            LP = B(1:2:length(B(:,1)),:);
            
            %timing step selection
            h = time/2^(n-m);
            sol = zeros(length(alpha),1+2*i); % initial value 0, only noise

            %reducing the floor of the components and computation of solution l^2 norm
            if i <= floor(2^((n-m)/2))
                for j = 1:(2^(n-m))
                    sol = sol + (i+1).^(-alpha.'/2).*(1-i*(i+1)*h)^(2^(n-m)-j).*LP(j,:);
                end;
            end;

            %computation of the error with MC methos of L^2(Omega,L^2(S^2)) of X(t_k)-X^(k,h)(t_k)
            for j = 1:(1+2*i)
                error(n-m,:) = error(n-m,:) + (final_time_sol_reference(:,j).' - sol(:,j).').^2;
            end;
        end;
    end;
end;

%final computation of the error and toc for the timing check
strong_error_MC_FEM = sqrt(error/M);
toc

%plot using the estimate h^{-\alpha/2}
figure();
loglog(2.^(1:1:(n-2)), 2.^(-(1:1:(n-2))*alpha(1)/4),'--','Color','#0072BD');
hold on;
loglog(2.^(1:1:(n-2)), strong_error_MC_FEM(1:(n-2),1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
loglog(2.^(1:1:(n-2)), 2.^(-(1:1:(n-2))*alpha(2)/4),':','Color','#77AC30');
loglog(2.^(1:1:(n-2)), strong_error_MC_FEM(1:(n-2),2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
loglog(2.^(1:1:(n-2)), 2.^(-(1:1:(n-2))*alpha(3)/4),'-.','Color','#D95319');
loglog(2.^(1:1:(n-2)), strong_error_MC_FEM(1:(n-2),3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
loglog(2.^(1:1:(n-2)), 2.^(-(1:1:(n-2))*alpha(4)/4),'-','Color','#1ec8c3');
loglog(2.^(1:1:(n-2)), strong_error_MC_FEM(1:(n-2),4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
loglog(2.^(1:1:(n-2)), 2.^(-(1:1:(n-2))*alpha(5)/4),'.-','Color','#fb8f36');
loglog(2.^(1:1:(n-2)), strong_error_MC_FEM(1:(n-2),5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
hold off;

title('Strong error forward EM with 10 MC samples')
xlabel('Number of time steps','Interpreter','latex')
ylabel('Strong error')
h_legend=legend('$O(h^{1/4})$','$\alpha = 1$','$O(h^{1/2})$','$\alpha = 2$','$O(h^{3/4})$','$\alpha = 3$','$O(h)$','$\alpha = 4$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
print -depsc2 -r0 strong_MC_fEM.eps
%---------------------------------------------------------------------------------------------------------------------------------

%% Backward EM

M = 50; % number of MC samples

time = 1;
n = 14; % time grid for reference solution with 2^n steps

% space parameter
alpha = [1 2 3 4 5]; % decay of angular power spectrum Q-Wiener process
kappa = floor(2^(n/2)); % truncation space

% error stored in vector for different time disretizations
strong_error_MC_BEM = zeros(n-1,length(alpha));  

error = zeros(n-1,length(alpha));

    for monte = 1:M
        monte
        for i = 1:kappa 
            i
            %selection of the noise and intensity
            option =2;
            lambda = 1;
            h = time/2^n;%selection of the timestep
            LP = levy(option, lambda, i, h,2^n)';    
           
            final_time_sol_reference = zeros(length(alpha),1+2*i);
            sol = zeros(length(alpha),1+2*i); % initial value 0, only nois
            
            for j = 1:(2^n)
                sol = sol + (i+1).^(-alpha.'/2)/(1+i*(i+1)*h)^(2^n+1-j).*LP(j,:);
            end;
            final_time_sol_reference = sol;
            
            for m = 1:(n-1)
                %reduction of the noise factor to the right size for the smaller solution 
                B = LP(1:(length(LP(:,1))-1),:) + LP(2:length(LP(:,1)),:);
                LP = B(1:2:length(B(:,1)),:);
                
                h = time/2^(n-m);
                sol = zeros(length(alpha),1+2*i); % initial value 0, only noise
                
                %if i <= floor(2^((n-m)/2))
                for j = 1:(2^(n-m))
                    sol = sol + (i+1).^(-alpha.'/2)/(1+i*(i+1)*h)^(2^(n-m)+1-j)*LP(j,:);
                end;
                %end;

                for j = 1:(1+2*i)
                    error(n-m,:) = error(n-m,:) + (final_time_sol_reference(:,j).' - sol(:,j).').^2;
                end;
            end;
        end;
    end;
    
    strong_error_MC_BEM = sqrt(error/M);

figure();
loglog(2.^(1:1:(n-2)), 2.^(-(1:1:(n-2))*alpha(1)/4),'--','Color','#0072BD');
hold on;
loglog(2.^(1:1:(n-2)), strong_error_MC_BEM(1:(n-2),1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
loglog(2.^(1:1:(n-2)), 2.^(-(1:1:(n-2))*alpha(2)/4),':','Color','#77AC30');
loglog(2.^(1:1:(n-2)), strong_error_MC_BEM(1:(n-2),2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
loglog(2.^(1:1:(n-2)), 2.^(-(1:1:(n-2))*alpha(3)/4),'-.','Color','#D95319');
loglog(2.^(1:1:(n-2)), strong_error_MC_BEM(1:(n-2),3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
loglog(2.^(1:1:(n-2)), 2.^(-(1:1:(n-2))*alpha(4)/4),'-','Color','#1ec8c3');
loglog(2.^(1:1:(n-2)), strong_error_MC_BEM(1:(n-2),4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
loglog(2.^(1:1:(n-2)), 2.^(-(1:1:(n-2))*alpha(5)/4),'.-','Color','#fb8f36');
loglog(2.^(1:1:(n-2)), strong_error_MC_BEM(1:(n-2),5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
hold off;

title('Strong error backward EM with 10 MC samples')
xlabel('Number of time steps','Interpreter','latex')
ylabel('Strong error')
h_legend=legend('$O(h^{1/4})$','$\alpha = 1$','$O(h^{1/2})$','$\alpha = 2$','$O(h^{3/4})$','$\alpha = 3$','$O(h)$','$\alpha = 4$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
print -depsc2 -r0 strong_MC_bEM.eps

%figure reduction for high order of approximation
figure();
loglog(2.^(1:1:(n-5)), 2.^(-(1:1:(n-5))*alpha(1)/4),'--','Color','#0072BD');
hold on;
loglog(2.^(1:1:(n-5)), strong_error_MC_BEM(1:(n-5),1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
loglog(2.^(1:1:(n-5)), 2.^(-(1:1:(n-5))*alpha(2)/4),':','Color','#77AC30');
loglog(2.^(1:1:(n-5)), strong_error_MC_BEM(1:(n-5),2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
loglog(2.^(1:1:(n-5)), 2.^(-(1:1:(n-5))*alpha(3)/4),'-.','Color','#D95319');
loglog(2.^(1:1:(n-5)), strong_error_MC_BEM(1:(n-5),3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
loglog(2.^(1:1:(n-5)), 2.^(-(1:1:(n-5))*alpha(4)/4),'-','Color','#1ec8c3');
loglog(2.^(1:1:(n-5)), strong_error_MC_BEM(1:(n-5),4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
loglog(2.^(1:1:(n-5)), 2.^(-(1:1:(n-5))*alpha(5)/4),'.-','Color','#fb8f36');
loglog(2.^(1:1:(n-5)), strong_error_MC_BEM(1:(n-5),5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
hold off;

%title('Strong error backward EM with 10 MC samples')
xlabel('Number of time steps','Interpreter','latex')
ylabel('Strong error')
h_legend=legend('$O(h^{1/4})$','$\alpha = 1$','$O(h^{1/2})$','$\alpha = 2$','$O(h^{3/4})$','$\alpha = 3$','$O(h)$','$\alpha = 4$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
print -depsc2 -r0 strong_MC_bEM.eps
