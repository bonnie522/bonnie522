% Code to make phase plots for the Kuramoto model with location dependent
% coupling as described by Abrams, Mirello, Strogatz, Wiley (2008)

b=0.1;         % phase lag parameter
A=0.2;              % Coupling parameter
mu = (1+A)/2;       % intergroup coupling
nu = (1-A)/2;       % intragroup coupling
alpha=pi/2-b;    % phase lag

% The Kuramoto model with location dependen coupling:
f = @(t,Y) [(1-Y(1)^2)/2*(Y(1)*mu*cos(-alpha) + nu*cos(Y(2)+alpha))
     (Y(1)^2+1)/(2*Y(1))*(Y(1)*mu*sin(-alpha)-nu*sin(Y(2)+alpha))-...
     mu*sin(-alpha)-Y(1)*nu*sin(Y(2)-alpha)];

% Initialising figure, axis, title, colours
figure;
xlabel('r cos \psi','FontSize',17);
ylabel('r sin \psi','FontSize',17);
title(['A=',num2str(A),' \beta=',num2str(b)]);
hold on
colour= [ 1 0 0; 1 0.5 0; 0.6 0 0.6; 0.2 0.6 1 ]; % red orange purple blue
i=1;    % to vary among the 4 colours

% Grid in terms of r and psi
rvec = [0.001 0.1:0.1:1];
psivec = [pi 0:0.2:pi 0:-0.2:-pi];

% store startpoints and endpoints to plot on top of trajectories
startpoints.y1 = zeros(length(rvec),length(psivec)); 
startpoints.y2 = zeros(length(rvec),length(psivec));
endpoints.y1 = zeros(length(rvec),length(psivec));
endpoints.y2 = zeros(length(rvec),length(psivec));

% running over grid, calculating and plotting trajectories
for r=1:length(rvec) %0:0.25:1
    for psi=1:length(psivec) %-pi:0.3:pi
        [ts,ys] = ode45(f,[0,3000],[rvec(r);psivec(psi)]);
        plot(ys(:,1).*cos(ys(:,2)),ys(:,1).*sin(ys(:,2)),'color',...
            colour(i,1:3)); i=i+1; if i>4 i=1; end
        startpoints.y1(r,psi) = ys(1,1).*cos(ys(1,2));
        startpoints.y2(r,psi) = ys(1,1).*sin(ys(1,2));
        endpoints.y1(r,psi) = ys(length(ys),1).*cos(ys(length(ys),2));
        endpoints.y2(r,psi) = ys(length(ys),1).*sin(ys(length(ys),2));
    end
end

% running over grid, plotting start- and endpoints
for r=1:length(rvec)
    for psi=1:length(psivec)
        plot(startpoints.y1(r,psi),...
            startpoints.y2(r,psi),'bo') % starting point
        plot(endpoints.y1(r,psi),...
            endpoints.y2(r,psi),'bd','MarkerFaceColor','b') % ending point
    end
end

hold off