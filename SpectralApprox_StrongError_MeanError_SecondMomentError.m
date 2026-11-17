%% Generates convergence plots for the spectral approximation 
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

% expectation and variance of the considered real Levy process
ExpLP=1;
VarLP=1;

% truncation index
k = 10;
kappa = 2^k;

% time interval
time = 1;

% decay rate
alpha = [1 2 3 4 5];

% initialize error vectors for every decay rate
spectral_strong_error = zeros(k,length(alpha));
spectral_mean_error = zeros(k,length(alpha));
spectral_second_moment_error = zeros(k,length(alpha));

% initialize variable for lower order truncation levels
j = k-1;

% initialize auxiliary error vectors for every decay rate
err_1 = zeros(1,length(alpha)); % for strong and second moment error
err_2 = zeros(1,length(alpha)); % for mean error

% compute the spectral error iterating over all frequencies
for i=0:(kappa-1)

    err_1 = err_1 + (1 - exp(-2*(kappa-i)*(kappa-i+1)*time))*(1+2*(kappa-i))/(2*(kappa-i)*(kappa-i+1))*(kappa-i).^(-alpha)*VarLP...
            + (1 - exp(-(kappa-i)*(kappa-i+1)*time))^2 *(1+2*(kappa-i))/((kappa-i)*(kappa-i+1))^2*(kappa-i).^(-alpha)*ExpLP^2;
    err_2 = err_2 + (1 - exp(-(kappa-i)*(kappa-i+1)*time))^2 *(1+2*(kappa-i))/((kappa-i)*(kappa-i+1))^2*(kappa-i).^(-alpha)*ExpLP^2;
    
    % enter errors into error vectors at truncation level
    if (kappa-i)==(2^j +1)
        spectral_strong_error(j+1,:) = sqrt(err_1);
        spectral_mean_error(j+1,:) = sqrt(err_2);
        spectral_second_moment_error(j+1,:) = err_1;
        
        % update truncation level
        j=j-1;
    end
end

%% Generate plots

% generate spectral strong error plot

fig_strong = figure('Units', 'centimeters', 'Position', [2 2 12 10]);
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(1)/2),'--','Color','#0072BD');
hold on;
loglog(2.^(0:1:(k-1)), spectral_strong_error(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(2)/2),':','Color','#77AC30');
loglog(2.^(0:1:(k-1)), spectral_strong_error(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(3)/2),'-.','Color','#D95319');
loglog(2.^(0:1:(k-1)), spectral_strong_error(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(4)/2),'-','Color','#1ec8c3');
loglog(2.^(0:1:(k-1)), spectral_strong_error(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(5)/2),'.-','Color','#fb8f36');
loglog(2.^(0:1:(k-1)), spectral_strong_error(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
hold off;

xlabel('Truncation level \kappa')
ylabel('Strong error')
h_legend=legend('$O(\kappa^{-1/2})$','$\alpha = 1$','$O(\kappa^{-1})$','$\alpha = 2$','$O(\kappa^{-3/2})$','$\alpha = 3$','$O(\kappa^{-2})$','$\alpha = 4$','$O(\kappa^{-5/2})$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');

% save spectral strong error data and plot in different formats
data_strong.k = k;
data_strong.alpha = alpha;
data_strong.spectral_strong_error = spectral_strong_error;
save('spectral_strong_error_data.mat', 'data_strong');
savefig(fig_strong, 'spectral_strong_error.fig');
exportgraphics(fig_strong, 'spectral_strong_error.pdf', 'ContentType', 'vector');
exportgraphics(fig_strong, 'spectral_strong_error.eps')

% generate spectral mean error plot

fig_mean = figure('Units', 'centimeters', 'Position', [2 2 12 10]);
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha(1)/2+1)),'--','Color','#0072BD');
hold on;
loglog(2.^(0:1:(k-1)), spectral_mean_error(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha(2)/2+1)),':','Color','#77AC30');
loglog(2.^(0:1:(k-1)), spectral_mean_error(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha(3)/2+1)),'-.','Color','#D95319');
loglog(2.^(0:1:(k-1)), spectral_mean_error(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha(4)/2+1)),'-','Color','#1ec8c3');
loglog(2.^(0:1:(k-1)), spectral_mean_error(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha(5)/2+1)),'.-','Color','#fb8f36');
loglog(2.^(0:1:(k-1)), spectral_mean_error(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
hold off;

xlabel('Truncation level \kappa')
ylabel('Mean error')
h_legend=legend('$O(\kappa^{-3/2})$','$\alpha = 1$','$O(\kappa^{-2})$','$\alpha = 2$','$O(\kappa^{-5/2})$','$\alpha = 3$','$O(\kappa^{-3})$','$\alpha = 4$','$O(\kappa^{-7/2})$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');

% save spectral mean error data and plot in different formats
data_mean.k = k;
data_mean.alpha = alpha;
data_mean.spectral_mean_error = spectral_mean_error;
save('spectral_mean_error_data.mat', 'data_mean');
savefig(fig_mean, 'spectral_mean_error.fig');
exportgraphics(fig_mean, 'spectral_mean_error.pdf', 'ContentType', 'vector');
exportgraphics(fig_mean, 'spectral_mean_error.eps')

% generate spectral second moment error plot

fig_second = figure('Units', 'centimeters', 'Position', [2 2 12 10]);
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(1)),'--','Color','#0072BD');
hold on;
loglog(2.^(0:1:(k-1)), spectral_second_moment_error(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(2)),':','Color','#77AC30');
loglog(2.^(0:1:(k-1)), spectral_second_moment_error(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(3)),'-.','Color','#D95319');
loglog(2.^(0:1:(k-1)), spectral_second_moment_error(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(4)),'-','Color','#1ec8c3');
loglog(2.^(0:1:(k-1)), spectral_second_moment_error(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*alpha(5)),'.-','Color','#fb8f36');
loglog(2.^(0:1:(k-1)), spectral_second_moment_error(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
hold off;

xlabel('Truncation level \kappa')
ylabel('Second moment error')
h_legend=legend('$O(\kappa^{-1})$','$\alpha = 1$','$O(\kappa^{-2})$','$\alpha = 2$','$O(\kappa^{-3})$','$\alpha = 3$','$O(\kappa^{-4})$','$\alpha = 4$','$O(\kappa^{-5})$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');

% save spectral second moment error data and plot in different formats
data_second.k = k;
data_second.alpha = alpha;
data_second.spectral_second_moment_error = spectral_second_moment_error;
save('spectral_second_moment_error_data.mat', 'data_second');
savefig(fig_mean, 'spectral_second_moment_error.fig');
exportgraphics(fig_mean, 'spectral_second_moment_error.pdf', 'ContentType', 'vector');
exportgraphics(fig_mean, 'spectral_second_moment_error.eps')
