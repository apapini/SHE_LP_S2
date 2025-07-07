% Strong, mean, second moment, Spectral convergence stochastic heat equation 
% with additive Lévy noise on the sphere
clear all 
close all

%set seed
randn('state',112);

% truncation index
k = 10;
kappa = 2^k;

% time interval
time = 1;

% noise regularity
alpha = [1 2 3 4 5];

%ExpLP=1, VarLP=1, then Strong Spectral Error;
%ExpLP=1, VarLP=0, then Mean Spectral Error;
%ExpLP=0, VarLP=1, then Second Moment Spectral Error;
ExpLP=1;
VarLP=1;

%vector of strong error
strong_error_kappa = zeros(k,length(alpha));

%variable to check the error against lower order truncation.
j = k-1;

%auxiliary vector of error computad for every regularity parameter.
err = zeros(1,length(alpha));

%loop to constract the teoretical spectral error.
    for i=0:(kappa-1)

        % initial condition = 0
        err = err + (1 - exp(-2*(kappa-i)*(kappa-i+1)*time))*(1+2*(kappa-i))/(2*(kappa-i)*(kappa-i+1))*(kappa-i).^(-alpha)*VarLP...
                + (1 - exp(-(kappa-i)*(kappa-i+1)*time))^2 *(1+2*(kappa-i))/((kappa-i)*(kappa-i+1))^2*(kappa-i).^(-alpha)*ExpLP^2;
        
        %loop to stop the troncation and place the error in the right position
        if (kappa-i)==(2^j +1)
            if ExpLP==0
                strong_error_kappa(j+1,:) = sqrt(err).^2;
            else
                %if VarLP==0
                %    strong_error_kappa(j+1,:) = sqrt(err);
                %else
                    strong_error_kappa(j+1,:) = sqrt(err);
                %end;
            end;
            j=j-1;
        end
    end

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
    
    title('Strong error spectral approximation')
    xlabel('Number of series elements $\kappa$','Interpreter','latex')
    ylabel('Second moment error')
    h_legend=legend('$O(\kappa^{1})$','$\alpha = 1$','$O(\kappa^2)$','$\alpha = 2$','$O(\kappa^{3})$','$\alpha = 3$','$O(\kappa^{4})$','$\alpha = 4$','$O(\kappa^{5})$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
    print -depsc2 -r0 strong_SmE_spectral.eps
else
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
        
        title('Strong error spectral approximation')
        xlabel('Number of series elements $\kappa$','Interpreter','latex')
        ylabel('Mean error')
        h_legend=legend('$O(\kappa^{1/2+1})$','$\alpha = 1$','$O(\kappa^{1+1})$','$\alpha = 2$','$O(\kappa^{3/2+1})$','$\alpha = 3$','$O(\kappa^{2+1})$','$\alpha = 4$','$O(\kappa^{5/2+1})$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
        print -depsc2 -r0 strong_ME_spectral.eps
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
        
        title('Strong error spectral approximation')
        xlabel('Number of series elements $\kappa$','Interpreter','latex')
        ylabel('Strong error')
        h_legend=legend('$O(\kappa^{1/2})$','$\alpha = 1$','$O(\kappa)$','$\alpha = 2$','$O(\kappa^{3/2})$','$\alpha = 3$','$O(\kappa^{2})$','$\alpha = 4$','$O(\kappa^{5/2})$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');
        print -depsc2 -r0 strong_spectral.eps
    end;
end;
