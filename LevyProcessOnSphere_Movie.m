%% Simulation of samples of a Levy process on the sphere
%
% Simulation of multiple samples of different kinds of Levy processes on the
% sphere
%
% The Levy process is defined as a series expansion with respect to the
% real-valued spherical harmonics Y_{l,m}, where the components are
% identically distributed real-valued Levy processe generated in levy.m 
% multiplied by the index l to the power of half the decay rate.
%
% Plot code by Denise L. Chen  9-1-93.
% in Template "Spherical Surface Harmonic"
%
% for the plots in the paper we used the following samples
% option 1 seed 119
% option 2 seed 122
% option 3 seed 113
% option 4 seed 118

% define the truncation of the harmonic series.
degree = 2^7;

% create the grid
delta = pi/480;
% altitude
theta = 0 : delta : pi;
% azimuth
phi = 0 : 2*delta : (2*pi - 2*delta);
% makes the grid
[phi,theta] = meshgrid(phi,theta);

% time window and timestepping
time = 1;
timesteps = 400;
h = time/timesteps;

% define the decay rate alpha 
decay = 4; % strictly larger than 2
Cl = ((0:degree)+1).^(-decay);
sqCl = sqrt(2*Cl);

% selection of type of 1-dim Levy noise 
option =[1,2,3,4];

% jump intensity parameter
lambda =1;

% can simulate multiple option if made as range
for j = option
     % can be increased to simulate the same setting a couple of times
    for i = 1:10

        rng(112+i)

        % noise process
        npL = zeros(size(phi,2),size(theta,1),(timesteps+1));
    
        % sum over all spherical harmonics up to order degree
        for l = degree:-1:0 
            
            % generate generalized Legendre polynomials
            % of degree l and scale them
            Plm = legendre(l,cos(theta(:,1)),'norm');
            scal = (1).^(0:1:l)*sqrt(1/pi);
            scal(1) = scal(1)/sqrt(2);
    
            Llm = Plm.*(ones(size(theta,1),1)*scal)';
            
            % iteration over all time steps
            for t=1:timesteps
    
                % generation of the noise increment
                % generate random fields
    
                numberRV =2*l+1;
    
                % samples of the Levy process
                LP = levy(j, lambda, l, h,1);
    
                % initialize random fields
                LRF = zeros(l+1,size(phi,2));
    
                % use samples to generate time increment of random fields
                LRF(1,:) = (LP(1)*ones(1,size(phi,2)));

                if l >0
                    LRF(2:end,:) = (LP(2:(l+1))*ones(1,size(phi,2))).*cos((1:1:l)'*phi(1,:)) +...
                        (LP((l+2):(2*l+1))*ones(1,size(phi,2))).*sin((1:1:l)'*phi(1,:));
                end
                
                % add simulated values to the increment process of the
                % stochastic noise
                npL(:,:,t+1) = npL(:,:,t+1) + (LRF'*Llm)*sqCl(l+1);
                
                clear LRF;
                clear BRF;
                clear PRF;
            end
            
        end
        
        % transform increment of noise to noise process
        for t = 2:timesteps
            npL(:,:,t+1) = npL(:,:,t) + npL(:,:,t+1);
        end
        
        % make video
        videoname = strcat('LPOption',string(j),'Num',string(i),'Color'); % enter the name as a text string
        vidobj = VideoWriter(videoname);
        vidobj.FrameRate = 10;
        open(vidobj)
        
        videoname2 = strcat('LPOption',string(j),'Num',string(i),'Deep'); % enter the name as a text string
        vidobj2 = VideoWriter(videoname2);
        vidobj2.FrameRate = 10;
        open(vidobj2)
        
        npL = npL/max(abs(min(min(min(npL)))), abs(max(max(max(npL)))));

        x = cos([phi 2*pi*ones(size(phi,1),1)]).*sin([theta theta(:,1)]);
        y = sin([phi 2*pi*ones(size(phi,1),1)]).*sin([theta theta(:,1)]);
        z = cos([theta theta(:,1)]);
       
        % make video of surface
        for t = 1 : timesteps+1
            
            % create a picture for every time step
            surf(x,y,z,[npL(:,:,t);npL(1,:,t)]','EdgeColor','none')
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
            
            % save as video
            writeVideo(vidobj, getframe(gcf))
        end    
        close(vidobj)
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
        fOut = sprintf(strcat('LPOption',string(j),'Num',string(i),'ColorP'),timesteps-1);  
        print('-djpeg', fOut, '-r100')
        
        % make video of sphere
        figure
        
        for t = 1 : timesteps+1
            
            % create a picture for every time step
            rho = exp([npL(:,:,t);npL(1,:,t)]');
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
        
            % save as video
            writeVideo(vidobj2, getframe(gcf))    
        end
        close(vidobj2)
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
        fOut = sprintf(strcat('LPOption',string(j),'Num',string(i),'DeepP'),timesteps-1);   % say, i=1
        print('-djpeg', fOut, '-r100')
    end
end
