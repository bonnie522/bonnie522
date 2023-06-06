% Code to run the Kuramoto model with different options:
% - different frequency distributions
% - location dependent coupling ('Chimera model')
% - add a white-noise forcing term (Monte Carlo method)
%
% A model can be chosen by changing 'modeltype' Either the parameter K
% (model 'k') or the parameters A and beta (model 'c') can then be changed
%
% The other parameters are the same for each model. 
%
% Change for example the number of oscillators 'N' or the runtime 'T'. The
% initial phases are should be initialized as an N by 1 array 'X0'. Bi- and
% unimodal discrete and Lorenzian frequency distributions can be specified 
% by 'w0','w0s', 'sd' and 'distribution'
%
% When adding a white-noise term (By changing 'D'), determine the number
% of trajectories by 'Na' and the number of trajectories to add per time
% by 'Nadd'. When you want to add more trajectories without resetting the 
% variables change 'reset' to 0. Change 'reset' back to 1 when wanting to 
% run a new experiment.
%
% The total number of trajectories produced (Monte Carlo) and the variances
% of the 3 different order parameters are displayed. Several figures are
% produced and saved under time specific names when the lines under "save
% figures" are commented out.
%
% Initializing 'plotpartial' to 1 will output extra plots for r1, r2, z1,
% z2 and their time-frequncy plots.


%%% ----- Parameters
D = 0;                      % noise strenght
T=100;                      % end time
N = 20;                     % number of oscillators (even)

%%% ----- Type of Model ----- %%%
modeltype = 'c';    % 'c': Chimera model, 'k': Kuramoto model 
plotpartial = 1;     % if 1, plots for partial order parameters 
                    % & time-frequency will be made

%%% Original Kuramoto 'k'
K = 4;            

%%% Chimera model 'c'
A = 0.2;%0.2;
beta = 0.1; %pi/2;%0.01;

%%% ----- Monte Carlo ----- %%% 
% Trajectories
Na = 10;              % #trajectories
Nadd = 10;        % number of trajectories to add each loop
reset = 1;              % to add new trajectories: set to 0

% Timestep
t0=0;                   % initial time
dt = 0.1;               % time step

%%% ----- Initialising oscillator phases ---- %%%
%%% Completely uniform: (~incoherent solution)
X0 = 2*pi.*rand(N,1);     

%%% close to chimera
X0(1:N/2)=zeros(N/2,1);
X0(N/2+1:N)=0.6*randn(N/2,1)+0.4;

X0 = mod(X0+pi,2*pi)-pi;

%X0=Xf(:,steps+1); % Continue from previous endpoint

%%% ----- Frequency distribution ----- %%%
w0 = 0.3;                 % frequency mean
w0s = 0.5;                % frequency mean 2
sd=0.025;                 % standard deviation
distribution = 'ud';    %  *types of distributions:
                        % 'bl' bimodal lorentz   'ul' unimodal lorenz
                        % 'bd' bimodal discrete  'ud' unimodal discrete
                        
if reset
   [w, Ntotal,sqordersum,sqordersum1,sqordersum2, steps, ordersum, ...
       ordersum1, ordersum2,phasesum]=initialize(N, distribution,...
       w0,w0s,sd,T,t0,dt);
end

tic
[Xf,Ntotal,Var,Var1,Var2,sqordersum,sqordersum1,sqordersum2, ordersum, ...
    ordersum1, ordersum2,traj1order,traj1order1,traj1order2,phasesum]...
    = RegularMonteCarlo(phasesum,modeltype,A,beta,Nadd,Na,dt,N,w,...
    D,K, X0, t0,T,Ntotal, ordersum, ordersum1, ordersum2,sqordersum,...
    sqordersum1,sqordersum2);
   t=0:dt:T;                       % time discretizsation
toc

disp(Ntotal)        % total number of trajectories simulated
disp([Var Var1 Var2]) 

figures = []; % to store figure handles
figures = plotswave(plotpartial,figures, ordersum,ordersum1, ...
    ordersum2,t,Ntotal,traj1order,traj1order1,traj1order2);

% %%% ----- Produce a scatterplot of the natural frequencies ----- %%% 
%     figures(8)=figure; scatter(1:length(w),sort(w))
%     
% %%% ----- Produce a scatterplot of the initial phases ----- %%%
%     figures(1) = figure; scatter(1:N,X0,'MarkerFaceColor','b')
%     ylabel('\theta_j'); xlabel('j'); ylim([-pi pi])

% %%% ----- Produce movie of the phases of the oscillators ----- %%%
% h=figure;
% 
% v=VideoWriter('phases.avi');
% open(v)
% for i = 1:length(Xf)
%     scatter(1:N,mod(Xf(:,i)+pi,2*pi)-pi,'MarkerFaceColor','b')
%     title(['r_1=',num2str(abs(ordersum1(i)),2), '     r_2=',...
%       num2str(abs(ordersum2(i)),2), '     t=',num2str(t(i)),]);
%     ylabel('\theta_j')
%     xlabel('j')
%     ylim([-pi pi])
% 
%     frame = getframe(h); 
%     writeVideo(v,frame);
% end
% close(v);    
    
% %%% ----- Saving figures ----- %%%
% ct=fix(clock);
% saveas(figures(1),[num2str(ct(3)),'0',num2str(ct(2)),'t',...
% num2str(ct(4)),num2str(ct(5)),'X0'],'epsc') %initial phases
% saveas(figures(2),[num2str(ct(3)),'0',num2str(ct(2)),'t',...
% num2str(ct(4)),num2str(ct(5)),'r'],'epsc') % abs order
% saveas(figures(3),[num2str(ct(3)),'0',num2str(ct(2)),'t',...
% num2str(ct(4)),num2str(ct(5)),'rc'],'epsc') % complex order
% saveas(figures(5),[num2str(ct(3)),'0',num2str(ct(2)),'t',...
% num2str(ct(4)),num2str(ct(5)),'rcw'],'epsc') % wave comp. order
% saveas(figures(8),[num2str(ct(3)),'0',num2str(ct(2)),'t',...
% num2str(ct(4)),num2str(ct(5)),'w'],'epsc') % freq distr
% % %%% Partial order parameters:
% saveas(figures(4),[num2str(ct(3)),'0',num2str(ct(2)),'t',...
% num2str(ct(4)),num2str(ct(5)),'r12'],'epsc') % partial abs order
% saveas(figures(6),[num2str(ct(3)),'0',num2str(ct(2)),'t',...
% num2str(ct(4)),num2str(ct(5)),'rc1w'],'epsc') % wave partial order 1
% saveas(figures(7),[num2str(ct(3)),'0',num2str(ct(2)),'t',...
% num2str(ct(4)),num2str(ct(5)),'rc2w'],'epsc') % wave partial order 2
% saveas(figures(9),[num2str(ct(3)),'0',num2str(ct(2)),'t',...
% num2str(ct(4)),num2str(ct(5)),'rc12'],'epsc') % complex partial orderparameters

function [w, Ntotal,sqordersum,sqordersum1,sqordersum2, steps, ...
    ordersum, ordersum1, ordersum2, phasesum]=...
    initialize(N, distribution,w0,w0s,sd,T,t0,dt)
    w = Distribution(distribution, N, w0,w0s, sd);
    
    Ntotal = 0;         % total number of trajectories simulated so far

    steps = (T-t0)/dt;  % number of timesteps
    
    ordersum=0;
    ordersum1=0;
    ordersum2=0;

    sqordersum = 0;
    sqordersum1 = 0;
    sqordersum2 = 0;
    
    phasesum = 0;
end

function [w] = Distribution(distribution, N, w0,w0s, sd)
    if distribution == 'bl'
        %%% bimodal lorenzian
        unif = rand(N/2,1);
        w1 = w0 + sd*tan(pi*(unif-0.5));
        w2 = w0s + sd*tan(pi*(unif-0.5));
        w = [w1 ; w2];
    elseif distribution == 'ul'    
        %%% unimodal lorenzian
        unif = rand(N,1);
        w = w0 + sd*tan(pi*(unif-0.5));
    elseif distribution == 'bd'
        %%% bimodal discrete
        w1 = w0*ones(N/2,1);
        w2 = w0s*ones(N/2,1);
        w = [w1 ; w2];  
    elseif distribution == 'ud'
        %%% unimodal discrete
        w = w0*ones(N,1);
    else
        disp('choose distribution bl, ul, bd or ud')
    end
end

function [g1] = Kuramoto(N,w,K, X)
    [~, s2]=size(X);
    g1=zeros(size(X));      % initialise
    for z=1:s2       % runs over trajectories
        [theta_j,theta_i] = meshgrid(X(:,z)); 
        g1(:,z) = w + K/N*sum(sin(theta_j-theta_i),2);
    end
end

function [g1] = KuramotoChimera(N,w,A,beta,X) % Abrams et al, 2008
    M=N/2; % number of oscillators per group
    [~, s2]=size(X);
    g1 = zeros(size(X));      % initialise
    for z=1:s2      % runs over trajectories
        [theta_j,theta_i] = meshgrid(X(:,z));
        g1(1:M,z) = w(1:M) - (1+A)/(2*M)*...
            sum(cos(theta_i(1:M,1:M)-theta_j(1:M,1:M)-beta),2)...
            -(1-A)/(2*M)*...
            sum(cos(theta_i(1:M,1:M)-theta_j(M+1:N,M+1:N)-beta),2);
        g1(M+1:N,z) = w(M+1:N) - (1+A)/(2*M)*...
            sum(cos(theta_i(M+1:N,M+1:N)-theta_j(M+1:N,M+1:N)-beta),2)...
            -(1-A)/(2*M)*...
            sum(cos(theta_i(M+1:N,M+1:N)-theta_j(1:M,1:M)-beta),2);
    end

end

function [Xf,Ntotal,Var,Var1,Var2,sqordersum,sqordersum1,...
    sqordersum2, ordersum, ordersum1, ordersum2,traj1order,...
    traj1order1,traj1order2,phasesum] = RegularMonteCarlo(phasesum,...
    modeltype,A,beta,Nadd,Na,dt,N,w,D,K, X0, t0,T,Ntotal, ...
    ordersum, ordersum1, ordersum2,sqordersum,sqordersum1,sqordersum2)
    %%% ----- Parameters ----- %%%    
    steps = (T-t0)/dt;      % number of timesteps
    
    %%% ----- Monte Carlo ----- %%%    
    while Ntotal < Na % && (Yvar>epsilon || Zvar>epsilon) 
        % initial positions of all oscillators
        Xf=ones(N, Nadd,steps+1).*X0;
        for j=1:steps
            % Brownian motion
            dW = sqrt(dt)*randn(N,Nadd);   	
            if modeltype == 'k'
                g1 = Kuramoto(N,w,K,Xf(:,:,j));	
            elseif modeltype == 'c'
                g1 = KuramotoChimera(N,w,A,beta,Xf(:,:,j)); 
            else
                error('Error. specify the type of model as c or k')
            end
            % Euler-Maruyana
            Xf(:,:,j+1) = Xf(:,:,j) + g1*dt + sqrt(2*D).*dW;  
        end
        
    %%% ----- Order parameter ----- %%%
    ett=exp(1j*Xf);
    
    %%% ----- Partial order parameters ----- %%%
    ett1 = exp(1j*Xf(1:N/2,:,:));
    ett2 = exp(1j*Xf(N/2+1:N,:,:));

%%%% ----- plotting different trajectories ---- %%%
%     ettplot = 1/N*permute(sum(ett,1),[3 1 2]);
%     figure
%     for kk=1:Nadd 
%         plot3(0:dt:T,real(ettplot(:,kk)),imag(ettplot(:,kk)));
%         hold on;
%         ylim([-1 1])
%         zlim([-1 1])
%     end
%     figure
%     for kk=1:Nadd     
%         plot(0:dt:T,abs(ettplot(:,kk)));
%         hold on;
%         ylim([0 1])
%     end
    
    %%% ---- Order parameter for 1st trajectory ----- %%%
    traj1order = 1/N*(sum(ett(:,1,:),1));
    traj1order1 = 1/(N/2)*sum(ett1(:,1,:),1);
    traj1order2 = 1/(N/2)*sum(ett2(:,1,:),1);
    
    phasesum = phasesum + sum(angle(1/N*(sum(ett(:,:,:),1))),2); 
    %%% ---- Calculate sums required for variance and mean ---- %%%
    ordersum= ordersum + sum(abs(1/N*(sum(ett(:,:,:),1))),2); 
    sqordersum= sqordersum + sum(abs(1/N*(sum(ett(:,:,steps+1),1))).^2,2);
    ordersum1 = ordersum1 + sum(abs(1/(N/2)*(sum(ett1(:,:,:),1))),2);
    ordersum2 = ordersum2 + ...
        sum(abs(1/(N/2)*(sum(ett2(:,:,:),1))),2);
    sqordersum1= sqordersum1 + ...
        sum(abs((1/(N/2)*(sum(ett1(:,:,steps+1),1)))).^2,2);
    sqordersum2= sqordersum2 + ...
        sum(abs((1/(N/2)*(sum(ett2(:,:,steps+1),1)))).^2,2);
    
    Ntotal = Ntotal + Nadd;      % update number of paths simulated   
    end
    
    %%% ---- Variances for absolute (partial) order parameters ---- %%%
    Var = abs(sum(sqordersum)/Ntotal - ...
        (sum(ordersum(1,1,steps+1))/Ntotal).^2);
    Var1 = abs(sum(sqordersum1)/Ntotal - ...
        (sum(ordersum1(1,1,steps+1))/Ntotal).^2);
    Var2 = abs(sum(sqordersum2)/Ntotal - ...
        (sum(ordersum2(1,1,steps+1))/Ntotal).^2);
    
    Xf = Xf(:,1,:);
end

function [figures]=plotswave(plotpartial, figures,ordersum,...
    ordersum1, ordersum2,t,Ntotal,traj1order,traj1order1,traj1order2) 
    
    %%% ---- Complex order parameters for 1st trajectory---- %%%
    ordercomplex = permute(traj1order,[3 2 1]);
    ordercomplex1 = permute(traj1order1,[3 2 1]);
    ordercomplex2 = permute(traj1order2,[3 2 1]);  
    
    %%% ---- Mean radius order parameter over all trajectories---- %%%
    orderabs = permute(ordersum,[3 2 1])/Ntotal;
    orderabs1 = permute(ordersum1,[3 2 1])/Ntotal;
    orderabs2 = permute(ordersum2,[3 2 1])/Ntotal;  
    
%   orderphase = permute(phasesum,[3 2 1])/Ntotal;
%   figure
%   plot3(t, real(orderabs.*exp(1j*orderphase)), ...
%   imag(orderabs.*exp(1j*orderphase)))

    %%% ----- Specifying figure sizes and positions ----- %%%
    bdwidth = 5; topbdwidth=30;
    set(0, 'Units', 'pixels');
    scnsize= get(0, 'ScreenSize');
    pos1 = [bdwidth, 2/3*scnsize(4) + bdwidth, ...
        scnsize(3)/2 - 2*bdwidth, scnsize(4)/3-(topbdwidth + bdwidth)];
    pos2 = [pos1(1) + scnsize(3)/2, pos1(2), pos1(3), pos1(4)];
    pos3 = [pos1(1), 0, pos1(3)-300, pos1(4)+250];
    pos4 = [pos3(1)+ scnsize(3)/3, pos3(2),pos3(3),pos3(4)];
    pos5 = [pos3(1)+ 2*scnsize(3)/3, pos3(2),pos3(3),pos3(4)];
    
%     %%% ----- Simple Plots ----- %%%
    % radius of order parameter
    figures(2) = figure('Position', pos1);
    plot(t,orderabs)
    xlabel('t')
    ylabel('r')
    ylim([0 1])

    % order parameter
 	  figures(3) = figure;
    plot3(t,real(ordercomplex), imag(ordercomplex))
    xlabel('t')
    ylabel('real(re^{i\psi})')
    zlabel('imag(re^{i\psi})')
    ylim([-1 1])
    zlim([-1 1])
    

    %%% ---- Wavelet Toolbox Time-Frequency analysis ---- %%%
    figures(5) = figure('Position', pos3);                  
    cwt(ordercomplex,'amor')      % 'amor' is the analytic morlet wavelet. 

    
    if plotpartial
        % partial order parameters
        figures(9) = figure;    
        plot3(t,real(ordercomplex1), imag(ordercomplex1))
        hold on;
        plot3(t,real(ordercomplex2), imag(ordercomplex2))
        xlabel('t')
        ylabel('real(re^{i\psi})')
        zlabel('imag(re^{i\psi})')
        ylim([-1 1])
        zlim([-1 1])

        % radius of partial order parameters
        figures(4) = figure('Position', pos2); 
        plot(subplot(2,1,1),t,orderabs1)
        xlabel('t')
        ylabel('r_1')
        ylim([0 1])
        hold on;
        plot(subplot(2,1,2),t,orderabs2)
        xlabel('t')
        ylabel('r_2')
        ylim([0 1])
        
        % wavelet toolbox time-frequency analysis for partial order
        % parameters
%         figures(6) = figure('Position', pos4);
%         cwt(ordercomplex1,'amor')
%         figures(7) = figure('Position', pos5);
%         cwt(ordercomplex2,'amor')
    end 
end
        


