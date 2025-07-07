
%% Backward EM

ExpLP = 0;
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
                err = err + (2*(1:2^k).'+1).*((1:2^k).'+1).^(-alpha) .*...
                       ((1-exp(-2*(1:2^k).'.*((1:2^k).'+1)*h))./(2*(1:2^k).'.*((1:2^k).'+1)).*exp(-2*(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j)) ...
                       - h./(1+(1:2^k).'.*((1:2^k).'+1)*h).^(2*(2^m-j+1)));
                err_intermediate = err_intermediate +  (2*(1:2^k).'+1).*((1:2^k).'+1).^(-alpha) .* ...
                      ((1-exp(-(1:2^k).'.*((1:2^k).'+1)*h))./((1:2^k).'.*((1:2^k).'+1)).*exp(-(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j))...
                        - h./(1+(1:2^k).'.*((1:2^k).'+1)*h).^(2^m-j+1));
                %err = err + (2*(1:2^k).'+1).*((1:2^k).'+1).^(-alpha) .*...
                       %((1-exp(-2*(1:2^k).'.*((1:2^k).'+1)*h))./(2*(1:2^k).'.*((1:2^k).'+1)).*exp(-2*(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j)) ...
                       %- 2*(1-exp(-(1:2^k).'.*((1:2^k).'+1)*h))./((1:2^k).'.*((1:2^k).'+1)).*exp(-(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j))./(1+(1:2^k).'.*((1:2^k).'+1)*h).^(2^m-j+1)...
                       %+ h./(1+(1:2^k).'.*((1:2^k).'+1)*h).^(2*(2^m-j+1)));
                %err_intermediate = err_intermediate +  (2*(1:2^k).'+1).*((1:2^k).'+1).^(-alpha) .* ...
                      %((1-exp(-(1:2^k).'.*((1:2^k).'+1)*h))./((1:2^k).'.*((1:2^k).'+1)).*exp(-(1:2^k).'.*((1:2^k).'+1)*h*(2^m-j))...
                        %- h./(1+(1:2^k).'.*((1:2^k).'+1)*h).^(2^m-j+1));
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


