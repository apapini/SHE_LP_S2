%% Simulation of samples of the stochastic heat equation and the underlying noise

%set seed
%randn('state',112); 
rng('default')
rng(112);

% Define the truncation of the harmonic series.
degree = 2^7;

% Create the grid
delta = pi/480; 
theta = 0 : delta : pi; % altitude
phi = 0 : 2*delta : (2*pi - 2*delta); % azimuth
[phi,theta] = meshgrid(phi,theta);

time = 0.1;
timesteps = 200; 
h = time/timesteps;

% Define the angular power spectrum
% which defines the covariance.
decay = 3; % strictly larger than 2
Cl = ((0:degree)+1).^(-decay);
sqCl = sqrt(2*Cl);


lambda = 1;

% To simulate the different samples, in the paper the third sample is used
for i = 1:3 

    % increments of the stochastic heat equation on the surface
    %driven by Wiener process + Poisson process
    TL = ones(size(phi,2),size(theta,1));
    %driven by Wiener process
    TB = ones(size(phi,2),size(theta,1));
    %driven by Poisson process
    TP = ones(size(phi,2),size(theta,1));

    % solution of the stochastic heat equation
    %driven by Wiener process + Poisson process
    solL = zeros(size(TL,1),size(TL,2),(timesteps+1));
    %driven by Wiener process
    solB = zeros(size(TB,1),size(TB,2),(timesteps+1));
    %driven by Poisson process
    solP = zeros(size(TP,1),size(TP,2),(timesteps+1));

    % noise processes
    %Wiener process + Poisson process
    npL = zeros(size(TL,1),size(TL,2),(timesteps+1));
    %Wiener process
    npB = zeros(size(TB,1),size(TB,2),(timesteps+1));
    %Poisson process
    npP = zeros(size(TP,1),size(TP,2),(timesteps+1));

    solL(:,:,1) = solL(:,:,1) + TL;
    solB(:,:,1) = solB(:,:,1) + TB;
    solP(:,:,1) = solP(:,:,1) + TP;

    % Sum over all spherical harmonics up to order degree
    for l = degree:-1:0 
        
        % Generate generalized Legendre polynomials
        % of degree l and scale them
        Plm = legendre(l,cos(theta(:,1)),'norm');
        scal = (1).^(0:1:l)*sqrt(1/pi);
        scal(1) = scal(1)/sqrt(2);

        Llm = Plm.*(ones(size(theta,1),1)*scal)';
        
        %Generate rough random field for initial value
        
        RF = zeros(l+1,size(phi,2));
        RF(1,:) = (randn(1)*ones(1,size(phi,2)));
        if l >0
            RF(2:end,:) = (randn(l,1)*ones(1,size(phi,2))).*cos((1:1:l)'*phi(1,:)) +...
                (randn(l,1)*ones(1,size(phi,2))).*sin((1:1:l)'*phi(1,:));
        end
        beta = 2; %roughness of initial condition

        TL = (RF'*Llm)*(l+1)^(-beta/2)*1; %last number to scale noise
        TB = (RF'*Llm)*(l+1)^(-beta/2)*1; %last number to scale noise
        TP = (RF'*Llm)*(l+1)^(-beta/2)*1; %last number to scale noise
    
        solL(:,:,1) = solL(:,:,1) + TL;
        solB(:,:,1) = solB(:,:,1) + TB;
        solP(:,:,1) = solP(:,:,1) + TP;
        
        %Iteration over all time steps
        for t=1:timesteps

            % generation of the noise increment
            %Generate random fields

            numberRV =2*l+1;

            %Samples for the Wiener process
            BM = sqrt(h)*randn(numberRV,1);
            %Samples for the Poisson process
            PP = poissrnd(lambda*h,numberRV,1)-lambda*h;
            %Samples for Wiener process + Poisson process
            LP = BM+PP;

            %Initialize random fields
            %Wiener process + Poisson process
            LRF = zeros(l+1,size(phi,2));
            %Wiener process
            BRF = zeros(l+1,size(phi,2));
            %Poisson process
            PRF = zeros(l+1,size(phi,2));

            %Use samples to generate time increment of random fields
            LRF(1,:) = (LP(1)*ones(1,size(phi,2)));
            BRF(1,:) = (BM(1)*ones(1,size(phi,2)));
            PRF(1,:) = (PP(1)*ones(1,size(phi,2)));
            if l >0
                LRF(2:end,:) = (LP(2:(l+1))*ones(1,size(phi,2))).*cos((1:1:l)'*phi(1,:)) +...
                    (LP((l+2):(2*l+1))*ones(1,size(phi,2))).*sin((1:1:l)'*phi(1,:));
                BRF(2:end,:) = (BM(2:(l+1))*ones(1,size(phi,2))).*cos((1:1:l)'*phi(1,:)) +...
                    (BM((l+2):(2*l+1))*ones(1,size(phi,2))).*sin((1:1:l)'*phi(1,:));
                PRF(2:end,:) = (PP(2:(l+1))*ones(1,size(phi,2))).*cos((1:1:l)'*phi(1,:)) +...
                    (PP((l+2):(2*l+1))*ones(1,size(phi,2))).*sin((1:1:l)'*phi(1,:));
            end

            % deterministic evolution of the heat equation
            % and addition of noise increment
            if l==0
                TL = TL + (LRF'*Llm)*sqCl(l+1);
                TB = TB + (BRF'*Llm)*sqCl(l+1);
                TP = TP + (PRF'*Llm)*sqCl(l+1);
            else
                %using the backward Euler scheme
                TL = (TL+  (LRF'*Llm) *sqCl(l+1))*(1+l*(l+1)*h).^(-1); 
                TB = (TB+  (LRF'*Llm) *sqCl(l+1))*(1+l*(l+1)*h).^(-1);
                TP = (TP+  (LRF'*Llm) *sqCl(l+1))*(1+l*(l+1)*h).^(-1);
            end

            %Add simulated values to solution process of the heat equation
            solL(:,:,t+1) = solL(:,:,t+1) + TL;
            solB(:,:,t+1) = solB(:,:,t+1) + TB;
            solP(:,:,t+1) = solP(:,:,t+1) + TP;
            
            %Add simulated values to the increment process of the
            %stochastic noise
            npL(:,:,t+1) = npL(:,:,t+1) + (LRF'*Llm)*sqCl(l+1);
            npB(:,:,t+1) = npB(:,:,t+1) + (BRF'*Llm)*sqCl(l+1);
            npP(:,:,t+1) = npP(:,:,t+1) + (PRF'*Llm)*sqCl(l+1);

            clear LRF;
            clear BRF;
            clear PRF;
        end
        
    end
    
    %Transform increment of noise to noise process
    for t = 2:timesteps
        npL(:,:,t+1) = npL(:,:,t) + npL(:,:,t+1);
        npB(:,:,t+1) = npB(:,:,t) + npB(:,:,t+1);
        npP(:,:,t+1) = npP(:,:,t) + npP(:,:,t+1);
    end
    
    %Create pictures and videos for the noise and the heat equation

    %Heat equation driven by Brownian motion
    videoname2 = strcat('Num',string(i),'SHEBM'); 
    vidobj2 = VideoWriter(videoname2);
    vidobj2.FrameRate = 10;
    open(vidobj2)
    
    caxis([min(min(min(solB))) max(max(max(solB)))]);
    
    solB = solB/max(abs(min(min(min(solB)))), abs(max(max(max(solB)))));
    
    for t = 1 : timesteps+1
         %create a picture for every time step
        rho = exp([solB(:,:,t);solB(1,:,t)]');
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
        camzoom(1.4)
        
        %save as video and store pictures of the interessting time points
        writeVideo(vidobj2, getframe(gcf))  
        if t==1
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHEBMtime001'),t);  
            print('-depsc', fOut, '-r100')
        end
        if t==2
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHEBMtime002'),t);  
            print('-depsc', fOut, '-r100')
        end
        if t==6
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHEBMtime006'),t);  
            print('-depsc', fOut, '-r100')
        end
        if t==11
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHEBMtime011'),t);  
            print('-depsc', fOut, '-r100')
        end
        if t==101
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHEBMtime101'),t);  
            print('-depsc', fOut, '-r100')
        end
    end
    close(vidobj2)
   
    %Brownian motion
    videoname3 = strcat('Num',string(i),'BM'); % enter the name as a text string
    vidobj3 = VideoWriter(videoname3);
    vidobj3.FrameRate = 10;
    open(vidobj3)
    
    caxis([min(min(min(npB))) max(max(max(npB)))]);
    
    npB = npB/max(abs(min(min(min(npB)))), abs(max(max(max(npB)))));
    
    for t = 1 : timesteps+1
        %create a picture for every time step
        
        rho = exp([npB(:,:,t);npB(1,:,t)]');
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
        camzoom(1.4)
        
        %save as video and store pictures of the interessting time points
        writeVideo(vidobj3, getframe(gcf))  
        if t==1
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'BMtime001'),t);  
            print('-depsc', fOut, '-r100')
        end
        if t==2
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'BMtime002'),t);  
            print('-depsc', fOut, '-r100')
        end
        if t==6
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'BMtime006'),t); 
            print('-depsc', fOut, '-r100')
        end
        if t==11
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'BMtime011'),t);  
            print('-depsc', fOut, '-r100')
        end
        if t==101
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'BMtime101'),t); 
            print('-depsc', fOut, '-r100')
        end
    end
    close(vidobj3)

    %Heat equation driven by Poisson process
    
    videoname4 = strcat('Num',string(i),'SHEPP'); % enter the name as a text string
    vidobj4 = VideoWriter(videoname4);
    vidobj4.FrameRate = 10;
    open(vidobj4)
    
    caxis([min(min(min(solP))) max(max(max(solP)))]);
    
    solP = solP/max(abs(min(min(min(solP)))), abs(max(max(max(solP)))));
    
    for t = 1 : timesteps+1
        %create a picture for every time step
        rho = exp([solP(:,:,t);solP(1,:,t)]');
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
        camzoom(1.4)
        
        %save as video and store pictures of the interessting time points
        writeVideo(vidobj4, getframe(gcf))  
        if t==1
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHEPPtime001'),t);  
            print('-depsc', fOut, '-r100')
        end
        if t==2
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHEPPtime002'),t); 
            print('-depsc', fOut, '-r100')
        end
        if t==6
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHEPPtime006'),t); 
            print('-depsc', fOut, '-r100')
        end
        if t==11
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHEPPtime011'),t); 
            print('-depsc', fOut, '-r100')
        end
        if t==101
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHEPPtime101'),t);  
            print('-depsc', fOut, '-r100')
        end
    end
    close(vidobj4)

    %Poisson process
    
    videoname5 = strcat('Num',string(i),'PP'); % enter the name as a text string
    vidobj5 = VideoWriter(videoname5);
    vidobj5.FrameRate = 10;
    open(vidobj5)
    
    caxis([min(min(min(npP))) max(max(max(npP)))]);
    
    npP = npP/max(abs(min(min(min(npP)))), abs(max(max(max(npP)))));
    
    for t = 1 : timesteps+1
        %create a picture for every time step
        rho = exp([npP(:,:,t);npP(1,:,t)]');
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
        camzoom(1.4)
        
        %save as video and store pictures of the interessting time points
        writeVideo(vidobj5, getframe(gcf))  
        if t==1
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'PPtime001'),t);  
            print('-depsc', fOut, '-r100')
        end
        if t==2
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'PPtime002'),t); 
            print('-depsc', fOut, '-r100')
        end
        if t==6
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'PPtime006'),t); 
            print('-depsc', fOut, '-r100')
        end
        if t==11
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'PPtime011'),t); 
            print('-depsc', fOut, '-r100')
        end
        if t==101
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'PPtime101'),t); 
            print('-depsc', fOut, '-r100')
        end
    end
    close(vidobj5)


    %Heat equation driven by Brownian motion + Poisson process
    
    videoname6 = strcat('Num',string(i),'SHELP'); % enter the name as a text string
    vidobj6 = VideoWriter(videoname6);
    vidobj6.FrameRate = 10;
    open(vidobj6)
    
    caxis([min(min(min(solL))) max(max(max(solL)))]);
    
    solL = solL/max(abs(min(min(min(solL)))), abs(max(max(max(solL)))));
    
    for t = 1 : timesteps+1
        %create a picture for every time step
        rho = exp([solL(:,:,t);solL(1,:,t)]');
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
        camzoom(1.4)
        
        %save as video and store pictures of the interessting time points
        writeVideo(vidobj6, getframe(gcf))  
        if t==1
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHELPtime001'),t); 
            print('-depsc', fOut, '-r100')
        end
        if t==2
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHELPtime002'),t); 
            print('-depsc', fOut, '-r100')
        end
        if t==6
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHELPtime006'),t); 
            print('-depsc', fOut, '-r100')
        end
        if t==11
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHELPtime011'),t);
            print('-depsc', fOut, '-r100')
        end
        if t==101
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'SHELPtime101'),t); 
            print('-depsc', fOut, '-r100')
        end
    end
    close(vidobj6)
    
    %Wiener process + Poisson process
    
    videoname7 = strcat('Num',string(i),'LP'); % enter the name as a text string
    vidobj7 = VideoWriter(videoname7);
    vidobj7.FrameRate = 10;
    open(vidobj7)
    
    caxis([min(min(min(npL))) max(max(max(npL)))]);
    
    npL = npL/max(abs(min(min(min(npL)))), abs(max(max(max(npL)))));
    
    for t = 1 : timesteps+1
        %create a picture for every time step
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
        camzoom(1.4)
        
        %save as video and store pictures of the interessting time points
        writeVideo(vidobj7, getframe(gcf))  
        if t==1
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'LPtime001'),t);
            print('-depsc', fOut, '-r100')
        end
        if t==2
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'LPtime002'),t);  
            print('-depsc', fOut, '-r100')
        end
        if t==6
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'LPtime006'),t); 
            print('-depsc', fOut, '-r100')
        end
        if t==11
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'LPtime011'),t); 
            print('-depsc', fOut, '-r100')
        end
        if t==101
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7])
            fOut = sprintf(strcat('Num',string(i),'LPtime101'),t); 
            print('-depsc', fOut, '-r100')
        end
    end
    close(vidobj3)

end
