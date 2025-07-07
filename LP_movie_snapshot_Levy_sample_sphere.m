%% Simulation of the stochastic wave equation (additive noise)
%
% Simulation of a surface which is a
% sphere with different log-normally distributed
% radius in every point
% 
% First a isotropic Gaussian random field
% is simulated where the expansion with respect
% to the spherical harmonics is truncated.
% Then, this field is exponentiated.
%
% plot code by Denise L. Chen  9-1-93.
% in Template "Spherical Surface Harmonic"

clear all

%set seed
randn('state',112); 

% Define the truncation of the harmonic series.
degree = 2^7; %2^7;

% Create the grid
delta = pi/480; delta = pi/160;% delta = pi/8;
%theta = 0 : delta : (pi-delta); % altitude
theta = 0 : delta : pi; % altitude
phi = 0 : 2*delta : (2*pi - 2*delta); % azimuth
%phi = 0 : 2*delta : 2*pi; % azimuth
[phi,theta] = meshgrid(phi,theta);

time = 1;
timesteps = 400;%400;
h = time/timesteps;

% Define the angular power spectrum
% which defines the covariance.
decay = 3; % strictly larger than 2
Cl = ((0:degree)+1).^(-decay);
sqCl = sqrt(2*Cl);

option =4;
if option==4
    lambda = 1;
else
    lambda=0.01;
end;

for j = option:option
    for i = 1:1 % To simulate the same setting a couple of times

        %random field on the surface
        T = zeros(size(phi,2),size(theta,1));
        sol1 = zeros(size(T,1),size(T,2),(timesteps+1)); % initial condition
        
        for t=1:timesteps
           
            T = zeros(size(phi,2),size(theta,1));
        
            % Sum over all spherical harmonics up to order degree
            for l = 0:degree
                % Generate generalized Legendre polynomials
                % of degree l and scale them
          
                Plm = legendre(l,cos(theta(:,1)),'norm');
                scal = (1).^(0:1:l)*sqrt(1/pi);
                scal(1) = scal(1)/sqrt(2);
               
                Llm = Plm.*(ones(size(theta,1),1)*scal)';
            
                %Generate random field
                LP = levy(option, lambda, l, h,1);
                LRF = zeros(l+1,size(phi,2));
                LRF(1,:) = (LP(1)*ones(1,size(phi,2)));
                if l >0
                    LRF(2:end,:) = (LP(2:(l+1))*ones(1,size(phi,2))).*cos((1:1:l)'*phi(1,:)) +...
                        (LP((l+2):(2*l+1))*ones(1,size(phi,2))).*sin((1:1:l)'*phi(1,:));
                end
                
                % Add the l-th component to the resulting
                % random field.
                % The covariance scaling sqCl is done here, since
                % an earlier scaling results in numerical
                % miscalculations due to very small values (???).
                T = T + (LRF'*Llm)*sqCl(l+1); %matrix is (phi,theta)
            end;
            sol1(:,:,t+1) = sol1(:,:,t) + T;
        
        end;
        
        %make video
        videoname = strcat('LPOption',string(j),'Num',string(i),'Color'); % enter the name as a text string
        vidobj = VideoWriter(videoname);
        vidobj.FrameRate = 10;
        open(vidobj)
        
        videoname2 = strcat('LPOption',string(j),'Num',string(i),'Deep'); % enter the name as a text string
        vidobj2 = VideoWriter(videoname2);
        vidobj2.FrameRate = 10;
        open(vidobj2)
        
        sol1 = sol1/max(abs(min(min(min(sol1)))), abs(max(max(max(sol1)))));

            x = cos([phi 2*pi*ones(size(phi,1),1)]).*sin([theta theta(:,1)]);
            y = sin([phi 2*pi*ones(size(phi,1),1)]).*sin([theta theta(:,1)]);
            z = cos([theta theta(:,1)]);
       
        %make video of surface
        for t = 1 : timesteps+1
            
            %create a picture for every time step
            surf(x,y,z,[sol1(:,:,t);sol1(1,:,t)]','EdgeColor','none')
            if t == 1
                v = caxis;
            end
            if t ~= 1
                caxis(v);
            end
            axis tight equal off
            view(40,30)
            camzoom(1.5)
            colormap jet(256)
            
            %save as video
            writeVideo(vidobj, getframe(gcf))
        end    
        close(vidobj)
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('LPOption',string(j),'Num',string(i),'ColorP'),timesteps-1);   % say, i=1
            print('-djpeg', fOut, '-r100')
        
        %make video of sphere
        figure
        
        for t = 1 : timesteps+1
             %create a picture for every time step
        
            rho = exp([sol1(:,:,t);sol1(1,:,t)]');
            r = rho.*sin([theta theta(:,1)]);
            x = r.*cos([phi 2*pi*ones(size(phi,1),1)]);
            y = r.*sin([phi 2*pi*ones(size(phi,1),1)]);
             z = rho.*cos([theta theta(:,1)]);
        
            clf
            surf(x,y,z,'EdgeColor','none')
            light
            lighting phong
            axis tight equal off
            view(40,30)
            camzoom(1.5)
        
            %save as video
            writeVideo(vidobj2, getframe(gcf))    
        end
        close(vidobj2)
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('LPOption',string(j),'Num',string(i),'DeepP'),timesteps-1);   % say, i=1
            print('-djpeg', fOut, '-r100')
    end
end
