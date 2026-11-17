%% Generates convergence plots for the temporal approximation 
%% using the direct simulation method for the backward Euler-Maruyama scheme
%% for the strong error, the mean error and the second moment error
%
% Assume the Levy process is defined as a series expansion with respect to the
% real-valued spherical harmonics Y_{l,m}, where the components are
% identically distributed real-valued Levy processe multiplied by the index
% l to the power of half the decay rate.
%
% Based on the expectation and variance of the real-valued Levy process the
% errors are computed exactly using the formulas derived at the beginning
% of the proofs of the respective theorems in the underlying paper.

clear all 
close all

%% Parameters

% expectation and variance of the considered real Levy process
ExpLP=1;
VarLP=1;

% time interval
time = 1;

% decay rate
alpha = [1 2 3 4 5];

% reference solution has 2^(2n) steps
n = 10;

% initialize error vectors for every decay rate
temporal_strong_error = zeros(n,length(alpha));
temporal_mean_error = zeros(n,length(alpha));
temporal_second_moment_error = zeros(n,length(alpha));

% iteration over different number of time discretization intervals
for k = 1:n
    
    % balancing of time and space discretization to ensure stability
    m = 2*k;
    h = time/2^m;
    
    % calculate the temporal approximation errors 
    % this step heavily exploits the knowledge about the structure of the
    % levy noise and calculations simplifying the exact error terms in that
    % case.
    temporal_strong_error(k,:) = sqrt(abs(sum( VarLP*(2*(1:2^k).'+1).*((1:2^k).').^(-alpha) .*(-2*(exp(-(1:2^k).'.*((1:2^k).'+1)*time).*(1+(1:2^k).'.*((1:2^k).'+1)*h).^(-2^m)-1)./...
        ((1:2^k).'.*((1:2^k).'+1)) .*...
        (exp((1:2^k).'.*((1:2^k).'+1)*h)-1)./...
        (1- ((1+(1:2^k).'.*((1:2^k).'+1)*h).*exp((1:2^k).'.*((1:2^k).'+1)*h))) +...
        (1-exp(-2*(1:2^k).'.*((1:2^k).'+1)*time))./(2*(1:2^k).'.*((1:2^k).'+1)) +...
        ((1+(1:2^k).'.*((1:2^k).'+1)*h).^(-2*2^m)-1)*h./(1-(1+(1:2^k).'.*((1:2^k).'+1)*h).^2))+...
        ExpLP^2* (2*(1:2^k).'+1).*((1:2^k).').^(-alpha).*(((1+(1:2^k).'.*((1:2^k).'+1)*h).^(-2^m)-exp(-(1:2^k).'.*((1:2^k).'+1)*time))./((1:2^k).'.*((1:2^k).'+1))).^2))); 

    temporal_mean_error(k,:) = sqrt(sum( ExpLP^2* (2*(1:2^k).'+1).*((1:2^k).').^(-alpha).*(((1+(1:2^k).'.*((1:2^k).'+1)*h).^(-2^m)-exp(-(1:2^k).'.*((1:2^k).'+1)*time))./...
        ((1:2^k).'.*((1:2^k).'+1))).^2)); 

    temporal_second_moment_error(k,:) = abs(sum( VarLP*(2*(1:2^k).'+1).*((1:2^k).').^(-alpha) .*...
        ((1-exp(-2*(1:2^k).'.*((1:2^k).'+1)*time))./(2*(1:2^k).'.*((1:2^k).'+1)) -...
        (1-(1+(1:2^k).'.*((1:2^k).'+1)*h).^(-2*2^m))*h./((1+(1:2^k).'.*((1:2^k).'+1)*h).^2 -1))+...
        ExpLP^2* (2*(1:2^k).'+1).*((1:2^k).').^(-alpha).*...
        (((1-exp(-(1:2^k).'.*((1:2^k).'+1)*time))./((1:2^k).'.*((1:2^k).'+1))).^2-...
        ((1-(1+(1:2^k).'.*((1:2^k).'+1)*h).^(-2^m))./((1:2^k).'.*((1:2^k).'+1))).^2))); 
end   
 

%% Generate plots

% generate temporal strong error plot

fig_strong = figure('Units', 'centimeters', 'Position', [2 2 12 10]);
loglog(2.^(2:2:(2*n)), 2.^(-(1:2:(2*n))*alpha(1)/4)/4,'--','Color','#0072BD');
hold on;
loglog(2.^(2:2:(2*n)), temporal_strong_error(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(2)/4)/4,':','Color','#77AC30');
loglog(2.^(2:2:(2*n)), temporal_strong_error(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(3)/4)/2,'-.','Color','#D95319');
loglog(2.^(2:2:(2*n)), temporal_strong_error(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(4)/4),'-','Color','#1ec8c3');
loglog(2.^(2:2:(2*n)), temporal_strong_error(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
%loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(5)/4)/2,'.-','Color','#fb8f36');
loglog(2.^(2:2:(2*n)), temporal_strong_error(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
hold off;

xlabel('Number of time steps')
ylabel('Strong error')
h_legend=legend('$O(h^{1/4})$','$\alpha = 1$','$O(h^{1/2})$','$\alpha = 2$','$O(h^{3/4})$','$\alpha = 3$','$O(h)$','$\alpha = 4$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');

% save spectral strong error data and plot in different formats
data_strong.n = n;
data_strong.alpha = alpha;
data_strong.temporal_strong_error = temporal_strong_error;
save('temporal_backward_strong_error_data.mat', 'data_strong');
savefig(fig_strong, 'temporal_backward_strong_error.fig');
exportgraphics(fig_strong, 'temporal_backward_strong_error.pdf', 'ContentType', 'vector');
exportgraphics(fig_strong, 'temporal_backward_strong_error.eps')


% generate temporal mean error plot

fig_mean = figure('Units', 'centimeters', 'Position', [2 2 12 10]);
loglog(2.^(2:2:(2*n)), 2.^(-(1:2:(2*n))*1)/4,'--','Color','#0072BD');
hold on;
loglog(2.^(2:2:(2*n)), temporal_mean_error(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
%loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*1)/4,':','Color','#77AC30');
loglog(2.^(2:2:(2*n)), temporal_mean_error(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
%loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*1)/4,'-.','Color','#D95319');
loglog(2.^(2:2:(2*n)), temporal_mean_error(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
%loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*1)/4,'-','Color','#1ec8c3');
loglog(2.^(2:2:(2*n)), temporal_mean_error(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
%loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(5)/4)/2,'.-','Color','#fb8f36');
loglog(2.^(2:2:(2*n)), temporal_mean_error(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
hold off;

xlabel('Number of time steps')
ylabel('Mean error')
h_legend=legend('$O(h)$','$\alpha = 1$','$\alpha = 2$','$\alpha = 3$','$\alpha = 4$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');

% save spectral mean error data and plot in different formats
data_mean.n = n;
data_mean.alpha = alpha;
data_mean.temporal_mean_error = temporal_mean_error;
save('temporal_backward_mean_error_data.mat', 'data_mean');
savefig(fig_mean, 'temporal_backward_mean_error.fig');
exportgraphics(fig_mean, 'temporal_backward_mean_error.pdf', 'ContentType', 'vector');
exportgraphics(fig_mean, 'temporal_backward_mean_error.eps')

% generate temporal second moment error plot

fig_second_moment = figure('Units', 'centimeters', 'Position', [2 2 12 10]);
loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*alpha(1)/2),'--','Color','#0072BD');
hold on;
loglog(2.^(2:2:(2*n)), temporal_second_moment_error(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*min(alpha(2)/2,1)),':','Color','#77AC30');
loglog(2.^(2:2:(2*n)), temporal_second_moment_error(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
%loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*min(alpha(3)/2,1)),'-.','Color','#D95319');
loglog(2.^(2:2:(2*n)), temporal_second_moment_error(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
%loglog(2.^(2:2:(2*n)), 2.^(-(2:2:(2*n))*min(alpha(4)/2,1)),'-','Color','#1ec8c3');
loglog(2.^(2:2:(2*n)), temporal_second_moment_error(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
%loglog(2.^(2:2:(2*n)), temporal_second_moment_error(1,3)*2.^(-(2:2:(2*n))*min(alpha(5)/4,1)),'.-','Color','#fb8f36');
loglog(2.^(2:2:(2*n)), temporal_second_moment_error(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
hold off;

xlabel('Number of time steps')
ylabel('Second moment error')
h_legend=legend('$O(h^{1/2})$','$\alpha = 1$','$O(h)$','$\alpha = 2$','$\alpha = 3$','$\alpha = 4$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');

% save spectral second moment error data and plot in different formats
data_second_moment.n = n;
data_second_moment.alpha = alpha;
data_second_moment.temporal_second_moment_error = temporal_second_moment_error;
save('temporal_backward_second_moment_error_data.mat', 'data_second_moment');
savefig(fig_second_moment, 'temporal_backward_second_moment_error.fig');
exportgraphics(fig_second_moment, 'temporal_backward_second_moment_error.pdf', 'ContentType', 'vector');
exportgraphics(fig_second_moment, 'temporal_backward_second_moment_error.eps')

