function [t, x, y, vx, vy] = sim()

% geometry
    L = 5e-9;
    W = 3e-9;

% physical constants
    e = 1.6e-19;
    m = 9.11e-31;
    epsilon_0 = 8.85e-12;
    konstant = -e^2 / m * 1 / (4*pi*epsilon_0);

% time grid
    tinit = 0;
    tend = 100e-15; % 10 fs

% initial conditions
    KE = rand() * e;
    theta = -pi/2 + pi * rand();
    vinit = sqrt(2*KE/m) * [cos(theta) sin(theta)]; % m/s
    xinit = -L/2;
    yinit = -W/2 + W * rand();
    rinit = [xinit yinit];
    Rinit = [rinit vinit];

% Matlab's ode45 (adaptive RK4)
    ode = @(t, R) Coloumb_Rprime(t, R, konstant, L, W);
    [t, R] = ode45(ode, [tinit tend], Rinit);
    x = R(:,1);
    y = R(:,2);
    vx = R(:,3);
    vy = R(:,4);

% plotting and animation
    clf;
    hold on;
    xlabel('x');
    ylabel('y');
    bounding_box = line(L/2*[1 -1 -1 1 1], W/2*[1 1 -1 -1 1], ...
                            'LineStyle','-', 'LineWidth',1, 'Color','k');
    proton = line(0, 0, 'color','r', 'Marker','.', 'markersize',50);
    text(0,0,'+e', 'color','w','HorizontalAlignment','center','VerticalAlignment','middle')
    electron = line(x(1), y(1), 'color','b', 'Marker','.', 'markersize',50);
    text(x(1),y(1), '-e', 'color','w','HorizontalAlignment','center','VerticalAlignment','middle')
    tra = plot(x,y,'b-', 'LineWidth',2);
    % arrow([x(1),y(1)], [x(1)+KE/e*W*cos(theta), y(1)+KE/e*W*sin(theta)]);
    if theta > 0
        text(x(1)+0.2e-9,y(1)-0.2e-9, ['KE = ', num2str(KE/e), ' eV']);
    else
        text(x(1)+0.2e-9,y(1)+0.2e-9, ['KE = ', num2str(KE/e), ' eV']);
    end
    axis(1.2*[-L/2 L/2 -W/2 W/2]);
    axis equal;
    for i = 1:length(t)
        set(tra, 'XData',x(1:i), 'YData',y(1:i))
        drawnow;
    end
end

% transformed system of first order ODE's
% Formed from Newtonian and Coloumb equations
function Rprime = Coloumb_Rprime(t, R, konstant, L, W)
    r = R(1:2);
    v = R(3:4);
    % stall the calculations if electron strays out of the bounding box
    if (abs(r(1)) > L/2) || (abs(r(2)) > W/2)
        Rprime = [zeros(size(R))];
        return;
    end
    rprime = v;
    vprime = konstant * r ./ (sqrt(sum(r.^2)))^3;
    Rprime = [rprime ; vprime];
    return;
end
