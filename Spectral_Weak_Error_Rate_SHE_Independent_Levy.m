
%set seed
randn('state',112);

%% Backward EM

M = 20; % number of MC samples

k=7; %troncation parameter of spectral series
kappa=2^k; %max spectral sum

time = 1; %global time
n = 16; % time grid for reference solution with 2^(n) steps
h = time/2^(n);%selection of the timestep

%selection of the noise and intensity
option = 2;%see levy.m for the different noise selections.
lambda = 1;%base intensity parameter for the poisson processes.

% convergence rate for different integrability parameter, i.e. p then select 0
% convergence rate for different regularity parameter, i.e alpha then select 1
pOrAlpDec=1; 

if pOrAlpDec==0
    alpha = [5]; % decay of angular power spectrum Q-Wiener process
    % selection of the noise and intensity
    p= [2 3 4 5 6];
    % error stored in vector for different spectral disretizations
    weak_error_MC_kappa = zeros(k,length(p));  
    error = zeros(k,length(p));
elseif pOrAlpDec==1 
    alpha = [1 2 3 4 5];
    % selection of the noise and intensity
    p= [3];
    % error stored in vector for different spectral disretizations
    weak_error_MC_kappa = zeros(k,length(alpha));  
    error = zeros(k,length(alpha));
end
    
auxerror=zeros(k,length(p));
for monte = 1:M
    
    if pOrAlpDec==0
        auxerror=zeros(k,length(p));
        aux=zeros(length(p),1);
    elseif pOrAlpDec==1
        auxerror=zeros(k,length(alpha));
        aux=zeros(length(alpha),1);        
    end;
    
    for i = 1:kappa
        LP = levy(option, lambda, i, h,2^n)';

        if pOrAlpDec==0
            sol = zeros(length(p), 1+2*i); % initial value 0, only noise 
        elseif pOrAlpDec==1
            sol = zeros(length(alpha), 1+2*i); % initial value 0, only noise
        end

        for j = 1:(2^(n))
            sol = sol + (i).^(-alpha.'/2).*LP(j,:)./(1+i*(i+1)*h)^(2^n+1-j);
        end;
          
        for m = 1:(1+2*i)
            aux=aux+sol(:,m).^2;
        end;
            
        for l = 0:(k-1)
            if i==2^(l)
               auxerror(l+1,:)=aux;
            end;
        end;
                        
    end;
        
    for l = 0:(k-1)
        error(l+1,:) = error(l+1,:) + (aux.'.^((p/2))-auxerror(l+1,:).^((p/2))); %check for the deca of weak error for phi 
    end;
end;

%weak error for the case of phi
weak_error_MC_kappa = abs(error/M);

if pOrAlpDec==1
    %PICTURES for alpha
    figure();
    loglog(2.^(0:1:k-1), weak_error_MC_kappa (1,1)*2.^(-(0:1:k-1)*(alpha(1))),'--','Color','#0072BD');
    hold on;
    loglog(2.^(0:1:(k-1)), weak_error_MC_kappa (:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
    loglog(2.^(0:1:k-1), weak_error_MC_kappa (1,2)*2.^(-(0:1:k-1)*(alpha(2))),':','Color','#77AC30');
    loglog(2.^(0:1:(k-1)), weak_error_MC_kappa (:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
    loglog(2.^(0:1:k-1), weak_error_MC_kappa (1,3)*2.^(-(0:1:k-1)*(alpha(3))),'-.','Color','#D95319');
    loglog(2.^(0:1:(k-1)), weak_error_MC_kappa (:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
    loglog(2.^(0:1:k-1), weak_error_MC_kappa (1,4)*2.^(-(0:1:k-1)*(alpha(4))),'-','Color','#1ec8c3');
    loglog(2.^(0:1:(k-1)), weak_error_MC_kappa (:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
    loglog(2.^(0:1:k-1), weak_error_MC_kappa (1,5)*2.^(-(0:1:k-1)*(alpha(5))),'.-','Color','#fb8f36');
    loglog(2.^(0:1:(k-1)), weak_error_MC_kappa (:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
    hold off;
    
    %title('Weak error spectral approximation, p=3, X^0=0')
    xlabel('Number of series elements $\kappa$','Interpreter','latex')
    ylabel('Weak error')
    h_legend=legend('$O(\kappa^{-1})$','$\alpha = 1$','$O(\kappa^{-2})$','$\alpha = 2$','$O(\kappa^{-3})$','$\alpha = 3$','$O(\kappa^{-4})$','$\alpha = 4$','$O(\kappa^{-5})$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
    print -depsc2 -r0 strong_expect_spectral.eps

elseif pOrAlpDec==0
    %PICTURES for p
    figure();
    loglog(2.^(0:1:k-1), weak_error_MC_kappa (1,1)*2.^(-(0:1:k-1)*(alpha)),'--','Color','#0072BD');
    hold on;
    loglog(2.^(0:1:(k-1)), weak_error_MC_kappa (:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
    %loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha)),':','Color','#77AC30');
    loglog(2.^(0:1:(k-1)), weak_error_MC_kappa (:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
    %loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha)),'-.','Color','#D95319');
    loglog(2.^(0:1:(k-1)), weak_error_MC_kappa (:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
    %loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha)),'-','Color','#1ec8c3');
    loglog(2.^(0:1:(k-1)), weak_error_MC_kappa (:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
    %loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha)),'.-','Color','#fb8f36');
    loglog(2.^(0:1:(k-1)), weak_error_MC_kappa (:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
    hold off;
    
    %title('Weak error spectral approximation, $\varphi(X(t))=||X(t)||^p_{L^2(S^2)}$, $X^0=0$','Interpreter','latex')
    xlabel('Number of series elements $\kappa$','Interpreter','latex')
    ylabel('Weak error')
    h_legend=legend('$O(\kappa^{-5})$','$p = 2$','$p = 3$','$p = 4$','$p = 5$','$p = 6$','Location','SouthWest','Interpreter','latex');
    print -depsc2 -r0 strong_expect_spectral.eps
end
