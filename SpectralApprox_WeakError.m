%% Generate convergence plots for the spectral approximation
%% for the weak error using the backward Euler-Maruyama method for simulation
%% and the L^2-norm to the power p as test function
%
% Assume the Levy process is defined as a series expansion with respect to the
% real-valued spherical harmonics Y_{l,m}, where the components are
% identically distributed real-valued Levy processe generated in levy.m  
% multiplied by the index l to the power of half the decay rate.
%
% We compute weak error bounds using Monte Carlo simulation based on the 
% backward Euler--Maruyama method for test functions of the form |u(x)|^p 
% where |.| represents the L^2 norm on the sphere.
%
% To generate the plots in the paper we use the values:
% option = 1; lambda = 1; pOrAlpDec = 0;
% option = 2; lambda = 1; pOrAlpDec = 0;
% option = 2; lambda = 1; pOrAlpDec = 1;
% Plot generation is only implemented for these settings.

clear all
close all

% selection of the noise and intensity
option = 1; % choose the corresponding Levy process from levy.m
lambda = 1; % jump intensity for the underlying poisson processes

% decide which parameter to vary in convergence rate plots
% select 0 for different powers p for test functions
% select 1 for different decay rates alpha 
pOrAlpDec=0; 

% set seed
rng(112);

% number of MC samples
M = 100; 

% truncation index
k=7; 
kappa=2^k; 

%time interval
time = 1; 

% time grid
n = 16; 
N = 2^n;
h = time/N;

% initialize variables alpha and p accordingly
if pOrAlpDec==0
    % decay rate alpha
    alpha = 5; 
    % initialize vector of p values for test functions
    p= [2 3 4 5 6];
    % compute length of p vector
    len = length(p);
elseif pOrAlpDec==1 
    % initialize vector of decay rates alpha
    alpha = [1 2 3 4 5];
    % fix p values for test function
    p= 3;
    % compute length of alpha vector
    len = length(alpha);
end


% initialize weak error matrix
spectral_weak_error = zeros(k,len);

for sample = 1:M

    % initialize variable for lower order truncation levels
    l = 0;
    
    % initialize auxiliary vectors for each sample
    
    % stores the L^2-norm of the sample for different truncation levels
    auxerror=zeros(k,len);

    % computes the L^2-norm iteratively
    aux=zeros(len,1); 
    
    % iterate over different levels and compute L^2-norm of coefficients
    for i = 1:kappa

        % generate Levy noise for the level
        LP = levy(option, lambda, i, h,N)';
        
        % matrix storing backward Euler approximation of
        % components of the series expansion
        sol = (i).^(-alpha.'/2).*((1+i*(i+1)*h).^(-(N:-1:1))*LP);

        % add value of the components to the L^2-norm
        aux=aux+sum(sol.^2,2);
        
        % store L^2-norm of relevant truncation levels
        if i==2^l
           auxerror(l+1,:)=aux.';
           l = l+1;
        end
                        
    end
    
    % add error contribution of the sample to error
    spectral_weak_error = spectral_weak_error + (aux.'.^((p/2))-auxerror(1:k,:).^((p/2)));
end

% divide by number of samples and take absolut value to obtain error
spectral_weak_error = abs(spectral_weak_error/M);

%% Generate plots 

% Wiener process with varying p

if pOrAlpDec ==0 && lambda == 1 && option ==1
    
    fig_WP_p = figure('Units', 'centimeters', 'Position', [2 2 12 10]);
    loglog(2.^(0:1:k-1), spectral_weak_error(1,1)*2.^(-(0:1:k-1)*(alpha)),'--','Color','#0072BD');
    hold on;
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
    %loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha)),':','Color','#77AC30');
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
    %loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha)),'-.','Color','#D95319');
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
    %loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha)),'-','Color','#1ec8c3');
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
    %loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha)),'.-','Color','#fb8f36');
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
    hold off;
        
    xlabel('Truncation level \kappa')
    ylabel('Weak error')
    h_legend=legend('$O(\kappa^{-5})$','$p = 2$','$p = 3$','$p = 4$','$p = 5$','$p = 6$','Location','SouthWest','Interpreter','latex');

    % save spectral mean error data and plot in different formats
    data_WP_p.k = k;
    data_WP_p.p = p;
    data_WP_p.alpha = alpha;
    data_WP_p.spectral_weak_error = spectral_weak_error;
    save('spectral_weak_error_WP_p_data.mat', 'data_WP_p');
    savefig(fig_WP_p, 'spectral_weak_error_WP_p.fig');
    exportgraphics(fig_WP_p, 'spectral_weak_error_WP_p.pdf', 'ContentType', 'vector');
    exportgraphics(fig_WP_p, 'spectral_weak_error_WP_p.eps')

% Poisson process with varying p

elseif pOrAlpDec ==0 && lambda == 1 && option ==2
    
    fig_PP_p = figure('Units', 'centimeters', 'Position', [2 2 12 10]);
    loglog(2.^(0:1:k-1), spectral_weak_error(1,1)*2.^(-(0:1:k-1)*(alpha)),'--','Color','#0072BD');
    hold on;
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
    %loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha)),':','Color','#77AC30');
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
    %loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha)),'-.','Color','#D95319');
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
    %loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha)),'-','Color','#1ec8c3');
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
    %loglog(2.^(0:1:k-1), 2.^(-(0:1:k-1)*(alpha)),'.-','Color','#fb8f36');
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
    hold off;
        
    xlabel('Truncation level \kappa')
    ylabel('Weak error')
    h_legend=legend('$O(\kappa^{-5})$','$p = 2$','$p = 3$','$p = 4$','$p = 5$','$p = 6$','Location','SouthWest','Interpreter','latex');

    % save spectral mean error data and plot in different formats
    data_PP_p.k = k;
    data_PP_p.p = p;
    data_PP_p.alpha = alpha;
    data_PP_p.spectral_weak_error = spectral_weak_error;
    save('spectral_weak_error_PP_p_data.mat', 'data_PP_p');
    savefig(fig_PP_p, 'spectral_weak_error_PP_p.fig');
    exportgraphics(fig_PP_p, 'spectral_weak_error_PP_p.pdf', 'ContentType', 'vector');
    exportgraphics(fig_PP_p, 'spectral_weak_error_PP_p.eps')

% Poisson process with varying alpha

elseif pOrAlpDec ==1 && lambda == 1 && option ==2
    
    fig_PP_alpha = figure('Units', 'centimeters', 'Position', [2 2 12 10]);
    loglog(2.^(0:1:k-1), spectral_weak_error(1,1)*2.^(-(0:1:k-1)*(alpha(1))),'--','Color','#0072BD');
    hold on;
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,1),'v','Color','#0072BD','MarkerSize',10,'MarkerFaceColor','#0072BD');
    loglog(2.^(0:1:k-1), spectral_weak_error(1,2)*2.^(-(0:1:k-1)*(alpha(2))),':','Color','#77AC30');
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,2),'diamond','Color','#77AC30','MarkerSize',10,'MarkerFaceColor','#77AC30');
    loglog(2.^(0:1:k-1), spectral_weak_error(1,3)*2.^(-(0:1:k-1)*(alpha(3))),'-.','Color','#D95319');
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,3),'^','Color','#D95319','MarkerSize',10,'MarkerFaceColor','#D95319');
    loglog(2.^(0:1:k-1), spectral_weak_error(1,4)*2.^(-(0:1:k-1)*(alpha(4))),'-','Color','#1ec8c3');
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,4),'o','Color','#1ec8c3','MarkerSize',10,'MarkerFaceColor','#1ec8c3');
    loglog(2.^(0:1:k-1), spectral_weak_error(1,5)*2.^(-(0:1:k-1)*(alpha(5))),'.-','Color','#fb8f36');
    loglog(2.^(0:1:(k-1)), spectral_weak_error(:,5),'>','Color','#fb8f36','MarkerSize',10,'MarkerFaceColor','#fb8f36');
    hold off;
        
    xlabel('Truncation level \kappa')
    ylabel('Weak error')
    h_legend=legend('$O(\kappa^{-1})$','$\alpha = 1$','$O(\kappa^{-2})$','$\alpha = 2$','$O(\kappa^{-3})$','$\alpha = 3$','$O(\kappa^{-4})$','$\alpha = 4$','$O(\kappa^{-5})$','$\alpha = 5$','Location','SouthWest','Interpreter','latex');

    % save spectral mean error data and plot in different formats
    data_PP_alpha.k = k;
    data_PP_alpha.p = p;
    data_PP_alpha.alpha = alpha;
    data_PP_alpha.spectral_weak_error = spectral_weak_error;
    save('spectral_weak_error_PP_alpha_data.mat', 'data_PP_alpha');
    savefig(fig_PP_alpha, 'spectral_weak_error_PP_alpha.fig');
    exportgraphics(fig_PP_alpha, 'spectral_weak_error_PP_alpha.pdf', 'ContentType', 'vector');
    exportgraphics(fig_PP_alpha, 'spectral_weak_error_PP_alpha.eps')
end    
