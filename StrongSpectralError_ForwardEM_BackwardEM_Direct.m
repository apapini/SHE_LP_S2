% Strong convergence stochastic heat equation with additive noise
% on the sphere

%Fix expectation and variance of the Levy Process considered
%For the exact computation we do not need more information

%Expectation of L(1)
ExpLP = 1;
%Variance of L(1)
VarLP = 1;

% truncation index
k = 10;
kappa = 2^k;

% time interval
time = 1;

% initial condition
% eta = 1;

% noise regularity
alpha = [1,2,3,4,5];

strong_error_kappa = zeros(k,length(alpha));

    j = k-1;
    err = 0;
    for i=0:(kappa-1)
        %initial condition
        %init = exp(-2*(kappa-i)*(kappa-i+1)*time)*(1+4*(kappa-i))*(1+(kappa-i)*(kappa-i+1))^(-eta);
        
        %just brownian motion
        %err = err + (1 - exp(-2*(kappa-i)*(kappa-i+1)*time))*(1+2*(kappa-i))/(2*(kappa-i)*(kappa-i+1))*(kappa-i).^(-alpha); %(1+ 4*(kappa -i))*(kappa-i+1)^(-beta)*exp(-2*(kappa-i)*(kappa-i+1)*time);
        
        %Levy process computation
        err = err + (1 - exp(-2*(kappa-i)*(kappa-i+1)*time))*(1+2*(kappa-i))/(2*(kappa-i)*(kappa-i+1))*(kappa-i).^(-alpha)*VarLP...
                + (1 - exp(-(kappa-i)*(kappa-i+1)*time))^2 *(1+2*(kappa-i))/((kappa-i)*(kappa-i+1))^2*(kappa-i).^(-alpha)*ExpLP^2;
        
        if (kappa-i)==(2^j +1)
            if ExpLP==0
                strong_error_kappa(j+1,:) = sqrt(err).^2;
            else
                strong_error_kappa(j+1,:) = sqrt(err);
            end;
            j=j-1;
        end
    end

if VarLP==0
    figure();
    loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha(1)/2+1)),'--','Color','#0072BD');
    hold on;
    loglog(2.^(0:1:(k-1)), strong_error_kappa(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
    loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha(2)/2+1)),':','Color','#77AC30');
    loglog(2.^(0:1:(k-1)), strong_error_kappa(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
    loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha(3)/2+1)),'-.','Color','#D95319');
    loglog(2.^(0:1:(k-1)), strong_error_kappa(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
    loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha(4)/2+1)),'-','Color','#1ec8c3');
    loglog(2.^(0:1:(k-1)), strong_error_kappa(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
    loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha(5)/2+1)),'.-','Color','#fb8f36');
    loglog(2.^(0:1:(k-1)), strong_error_kappa(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
    hold off;
    
    %title('Convergence of the mean')
    xlabel('Number of series elements $\kappa$','Interpreter','latex')
    ylabel('Mean error')
    h_legend=legend('$O(\kappa^{-1/2-1})$','$\alpha = 1$','$O(\kappa^{-2/2-1})$','$\alpha = 2$','$O(\kappa^{-3/2-1})$','$\alpha = 3$','$O(\kappa^{-4/2-1})$','$\alpha = 4$','$O(\kappa^{-5/2-1})$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
    print -depsc2 -r0 strong_expect_spectral.eps
else
    if ExpLP==0
        figure();
        loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(1)),'--','Color','#0072BD');
        hold on;
        loglog(2.^(0:1:(k-1)), strong_error_kappa(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
        loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(2)),':','Color','#77AC30');
        loglog(2.^(0:1:(k-1)), strong_error_kappa(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
        loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(3)),'-.','Color','#D95319');
        loglog(2.^(0:1:(k-1)), strong_error_kappa(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
        loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(4)),'-','Color','#1ec8c3');
        loglog(2.^(0:1:(k-1)), strong_error_kappa(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
        loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(5)),'.-','Color','#fb8f36');
        loglog(2.^(0:1:(k-1)), strong_error_kappa(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
        hold off;
        
        %title('Strong error spectral approximation')
        xlabel('Number of series elements $\kappa$','Interpreter','latex')
        ylabel('Second moment error')
        h_legend=legend('$O(\kappa^{-1})$','$\alpha = 1$','$O(\kappa^{-2})$','$\alpha = 2$','$O(\kappa^{-3})$','$\alpha = 3$','$O(\kappa^{-4})$','$\alpha = 4$','$O(\kappa^{-5})$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
        print -depsc2 -r0 strong_expect_spectral.eps
    else
        figure();
        loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(1)/2),'--','Color','#0072BD');
        hold on;
        loglog(2.^(0:1:(k-1)), strong_error_kappa(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
        loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(2)/2),':','Color','#77AC30');
        loglog(2.^(0:1:(k-1)), strong_error_kappa(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
        loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(3)/2),'-.','Color','#D95319');
        loglog(2.^(0:1:(k-1)), strong_error_kappa(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
        loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(4)/2),'-','Color','#1ec8c3');
        loglog(2.^(0:1:(k-1)), strong_error_kappa(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
        loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(5)/2),'.-','Color','#fb8f36');
        loglog(2.^(0:1:(k-1)), strong_error_kappa(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
        hold off;
        
        %title('Strong error spectral approximation')
        xlabel('Number of series elements $\kappa$','Interpreter','latex')
        ylabel('Strong error')
        h_legend=legend('$O(\kappa^{-1/2})$','$\alpha = 1$','$O(\kappa^{-1})$','$\alpha = 2$','$O(\kappa^{-3/2})$','$\alpha = 3$','$O(\kappa^{-2})$','$\alpha = 4$','$O(\kappa^{-5/2})$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
        print -depsc2 -r0 strong_expect_spectral.eps
    end;

end;

%% EM

% time discretization
% time interval
time = 1;
n = 10; %up to 2^20 but just every 2nd step

% error stored in vector for different time disretizations
strong_error_EM = zeros(n,length(alpha));  

    for k = 1:n
        m = 2*k;
        h = time/2^m;
        err = zeros(2^k,length(alpha));
        err_intermediate = zeros(2^k,length(alpha));
        %init = (1+4*i)*(1+i*(i+1))^(-eta)*(exp(-i*(i+1)*time) - (1-i*(i+1)*h)^(2^m))^2;
        %init=0;

            for j = 1:2^m
                err = err + VarLP*(2*(1:2^k).'+1).*((1:2^k).'+1).^(-alpha) .* ...
                      ((1-exp(-2*(1:2^k).'.*((1:2^k).'+1)*h))./(2*(1:2^k).'.*((1:2^k).'+1)).*exp(-2*(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j)) ...
                        - 2*(1-exp(-(1:2^k).'.*((1:2^k).'+1)*h))./((1:2^k).'.*((1:2^k).'+1)).*exp(-(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j)).*(1-(1:2^k).'.*((1:2^k).'+1)*h).^(2^m-j)...
                        + h*(1-(1:2^k).'.*((1:2^k).'+1)*h).^(2*(2^m-j)));
                err_intermediate = err_intermediate +  sqrt((2*(1:2^k).'+1)).*((1:2^k).'+1).^(-alpha/2) .* ...
                      ((1-exp(-(1:2^k).'.*((1:2^k).'+1)*h))./((1:2^k).'.*((1:2^k).'+1)).*exp(-(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j))...
                        - h*(1-(1:2^k).'.*((1:2^k).'+1)*h).^(2^m-j));
            end
        strong_error_EM(k,:) = sqrt(sum(err)+ExpLP^2*sum(err_intermediate.^2));
    end
             

figure();
loglog(2.^(2:2:(2*n)), 2.^(-(1:2:(2*n))*alpha(1)/4)/2,'--','Color','#0072BD');
hold on;
loglog(2.^(2:2:(2*n)), strong_error_EM(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(2)/4)/2,':','Color','#77AC30');
loglog(2.^(2:2:(2*n)), strong_error_EM(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(3)/4)/2,'-.','Color','#D95319');
loglog(2.^(2:2:(2*n)), strong_error_EM(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(4)/4)/2,'-','Color','#1ec8c3');
loglog(2.^(2:2:(2*n)), strong_error_EM(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
%loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(5)/4)/2,'.-','Color','#fb8f36');
loglog(2.^(2:2:(2*n)), strong_error_EM(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
hold off;

title('Strong error forward EM')
xlabel('Number of time steps','Interpreter','latex')
ylabel('Strong error')
h_legend=legend('$O(h^{1/4})$','$\alpha = 1$','$O(h^{1/2})$','$\alpha = 2$','$O(h^{3/4})$','$\alpha = 3$','$O(h)$','$\alpha = 4$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
print -depsc2 -r0 strong_expect_fEM.eps


%% Backward EM
%% computation of spectral rate for strong, mean and second moment via direct computation.

%ExpLP=1, VarLP=1, then Strong Spectral Error;
%ExpLP=1, VarLP=0, then Mean Spectral Error;
%ExpLP=0, VarLP=1, then Second Moment Spectral Error;

ExpLP = 1;
VarLP = 1;

k = 10;
kappa = 2^k;

time = 1;
alpha = [1,2,3,4,5];
n = 10; 
strong_error_bEM = zeros(n,length(alpha));  

    for k = 1:n
        
        m = 2*k;
        h = time/2^m;
        err = zeros(2^k,length(alpha));
        err_intermediate = zeros(2^k,length(alpha));
            for j = 1:2^m
                
                if ExpLP==0
                    err = err + (2*(1:2^k).'+1).*((1:2^k).'+1).^(-alpha) .*...
                        ((1-exp(-2*(1:2^k).'.*((1:2^k).'+1)*h))./(2*(1:2^k).'.*((1:2^k).'+1)).*exp(-2*(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j)) ...
                        - h./(1+(1:2^k).'.*((1:2^k).'+1)*h).^(2*(2^m-j+1)));
                    err_intermediate = err_intermediate +  (2*(1:2^k).'+1).*((1:2^k).'+1).^(-alpha) .* ...
                        ((1-exp(-(1:2^k).'.*((1:2^k).'+1)*h))./((1:2^k).'.*((1:2^k).'+1)).*exp(-(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j))...
                        - h./(1+(1:2^k).'.*((1:2^k).'+1)*h).^(2^m-j+1));
                else
                    err = err + (2*(1:2^k).'+1).*((1:2^k).'+1).^(-alpha) .*...
                       ((1-exp(-2*(1:2^k).'.*((1:2^k).'+1)*h))./(2*(1:2^k).'.*((1:2^k).'+1)).*exp(-2*(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j)) ...
                       - 2*(1-exp(-(1:2^k).'.*((1:2^k).'+1)*h))./((1:2^k).'.*((1:2^k).'+1)).*exp(-(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j))./(1+(1:2^k).'.*((1:2^k).'+1)*h).^(2^m-j+1)...
                       + h./(1+(1:2^k).'.*((1:2^k).'+1)*h).^(2*(2^m-j+1)));
                    err_intermediate = err_intermediate +  (2*(1:2^k).'+1).*((1:2^k).'+1).^(-alpha) .* ...
                      ((1-exp(-(1:2^k).'.*((1:2^k).'+1)*h))./((1:2^k).'.*((1:2^k).'+1)).*exp(-(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j))...
                        - h./(1+(1:2^k).'.*((1:2^k).'+1)*h).^(2^m-j+1));
                end;
            end
            if ExpLP==0
                strong_error_bEM(k,:) = sqrt(sum(err)+ExpLP^2*sum(err_intermediate.^2)).^2;
            else
                if VarLP==0
                    strong_error_bEM(k,:) = sqrt(ExpLP^2*sum(err_intermediate.^2));
                else
                    strong_error_bEM(k,:) = sqrt(sum(err)+ExpLP^2*sum(err_intermediate.^2));
                end;
            end;
    end             

if ExpLP==0
figure();
loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(1)/2),'--','Color','#0072BD');
hold on;
loglog(2.^(2:2:(2*n)), strong_error_bEM(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*min(alpha(2)/2,1)),':','Color','#77AC30');
loglog(2.^(2:2:(2*n)), strong_error_bEM(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
%loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*min(alpha(3)/2,1)),'-.','Color','#D95319');
loglog(2.^(2:2:(2*n)), strong_error_bEM(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
%loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*min(alpha(4)/2,1)),'-','Color','#1ec8c3');
loglog(2.^(2:2:(2*n)), strong_error_bEM(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
loglog(2.^(2:2:(2*n)), strong_error_bEM(1,3)*2.^(-(2:2:(2*n))*min(alpha(5)/4,1)),'.-','Color','#fb8f36');
loglog(2.^(2:2:(2*n)), strong_error_bEM(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
hold off;

%title('Strong error backward EM')
xlabel('Number of time steps','Interpreter','latex')
ylabel('Second moment error')
h_legend=legend('$O(h^{1/2})$','$\alpha = 1$','$O(h)$','$\alpha = 2$','$\alpha = 3$','$\alpha = 4$','$O(h)$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
print -depsc2 -r0 strong_expect_bEM.eps
else
    if VarLP==0
        figure();
        loglog(2.^(2:2:(2*n)), 2.^(-(1:2:(2*n))*1)/4,'--','Color','#0072BD');
        hold on;
        loglog(2.^(2:2:(2*n)), strong_error_bEM(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
        %loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*1)/4,':','Color','#77AC30');
        loglog(2.^(2:2:(2*n)), strong_error_bEM(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
        %loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*1)/4,'-.','Color','#D95319');
        loglog(2.^(2:2:(2*n)), strong_error_bEM(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
        %loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*1)/4,'-','Color','#1ec8c3');
        loglog(2.^(2:2:(2*n)), strong_error_bEM(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
        %loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(5)/4)/2,'.-','Color','#fb8f36');
        loglog(2.^(2:2:(2*n)), strong_error_bEM(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
        hold off;
        
        %title('Strong error backward EM')
        xlabel('Number of time steps','Interpreter','latex')
        ylabel('Mean error')
        h_legend=legend('$O(h)$','$\alpha = 1$','$\alpha = 2$','$\alpha = 3$','$\alpha = 4$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
        print -depsc2 -r0 strong_expect_bEM.eps
        
        else
        figure();
        loglog(2.^(2:2:(2*n)), 2.^(-(1:2:(2*n))*alpha(1)/4)/4,'--','Color','#0072BD');
        hold on;
        loglog(2.^(2:2:(2*n)), strong_error_bEM(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
        loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(2)/4)/4,':','Color','#77AC30');
        loglog(2.^(2:2:(2*n)), strong_error_bEM(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
        loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(3)/4)/4,'-.','Color','#D95319');
        loglog(2.^(2:2:(2*n)), strong_error_bEM(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
        loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(4)/4)/4,'-','Color','#1ec8c3');
        loglog(2.^(2:2:(2*n)), strong_error_bEM(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
        %loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(5)/4)/2,'.-','Color','#fb8f36');
        loglog(2.^(2:2:(2*n)), strong_error_bEM(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
        hold off;
        
        %title('Strong error backward EM')
        xlabel('Number of time steps','Interpreter','latex')
        ylabel('Strong error')
        h_legend=legend('$O(h^{1/4})$','$\alpha = 1$','$O(h^{1/2})$','$\alpha = 2$','$O(h^{3/4})$','$\alpha = 3$','$O(h)$','$\alpha = 4$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
        print -depsc2 -r0 strong_expect_bEM.eps
    end;
end;
