function [c_l, c_d, V, s] = VORTEX_PANEL_FUNCTION(Xb, Yb, Vinf, c, angle_of_attack)
    %% Variables declerations
    Mp = length(Xb); % Make sure there is the same number of x and y data values
    M = Mp - 1; % Next iteration of M
    angle_of_attack = (angle_of_attack * pi) /180; % Convert to radians

    %% Coordinates (x,y) of control point and panel length S are computed for each of the vortex panels
    for i = 1:M 
        ip = i + 1; % Next i
        x(i) = 0.5 * (Xb(i) + Xb(ip));
        y(i) = 0.5 * (Yb(i) + Yb(ip));
        s(i) = sqrt((Xb(ip) - Xb(i))^2 + (Yb(ip) - Yb(i))^2); 

        theta(i) = atan2((Yb(ip) - Yb(i)), (Xb(ip) - Xb(i))); % Calculate theta
        sine(i) = sin(theta(i)); % Calculate sine(theta)
        cosine(i) = cos(theta(i)); % Calculate cosine(theta)
        RHS(i) = sin(theta(i) - angle_of_attack); % RHS represents the right-hand side of Eq. (5.47)
    end
    for i = 1:M
        for j = 1:M
            if (i == j)
                CN1(i, j) = -1.0;
                CN2(i, j) = 1.0;
                CT1(i, j) = 0.5 * pi;
                CT2(i, j) = 0.5 * pi;
            else
                % Calculate the geometric constants 
                A = -(x(i) - Xb(j)) * cosine(j) - (y(i) - Yb(j)) * sine(j);
                B = (x(i) - Xb(j))^2 + (y(i) - Yb(j))^2;
                C = sin(theta(i) - theta(j));
                D = cos(theta(i) - theta(j));
                E = (x(i) - Xb(j)) * sine(j) - (y(i) - Yb(j)) * cosine(j);
                F = log(1.0 + s(j) * (s(j) + 2 .* A) / B);
                G = atan2(E * s(j), B + A * s(j));
                P = (x(i) - Xb(j)) * sin(theta(i) - 2 .* theta(j)) + (y(i) - Yb(j)) * cos(theta(i) - 2 .* theta(j));
                Q = (x(i) - Xb(j)) * cos(theta(i) - 2 .* theta(j)) - (y(i) - Yb(j)) * sin(theta(i) - 2 .* theta(j));
                CN2(i, j) = D + 0.5 * Q * F / (s(j)) - (A * C + D * E) * G / (s(j));
                CN1(i, j) = 0.5 * D * F + C * G - CN2(i, j);
                CT2(i, j) = C + 0.5 * P * F / (s(j)) + (A * D - C * E) * G / (s(j));
                CT1(i, j) = 0.5 * C * F - D * G - CT2(i, j);
            end
        end
    end
    
    %% Compute Influence Coefficients 
    for i = 1:M
       AN(i, 1)  = CN1(i, 1);
       AN(i, Mp) = CN2(i, M);
       AT(i, 1)  = CT1(i, 1);
       AT(i, Mp) = CT2(i, M);
       for j = 2:M
          AN(i, j) = CN1(i, j) + CN2(i, j - 1);
          AT(i, j) = CT1(i, j) + CT2(i, j - 1);
       end
    end
    AN(Mp, 1) = 1.0;
    AN(Mp, Mp) = 1.0;
    for j = 2:M
       AN(Mp, j) = 0.0;
    end
    RHS(Mp) = 0.0;

    %% Solve the Linear System
    Gama = inv(AN) * (RHS'); % Same as multiplying by the inverse (X = inv(A)*B)    
    for i = 1:M
        V(i) = cos(theta(i) - angle_of_attack); % Tangential velocity
        for j = 1:Mp
           V(i) = V(i) + AT(i, j) * Gama(j); % Calculate the Velocity ratio (V/Vinf) or tangential velocity
           CP(i) = 1.0 - V(i)^2; % Calculate the coeff. of pressure
        end
    end
    
    %% Calculate the Circulation
    Circ = sum(V .* s);
    
    %% Calculate the Sectional Coefficient of Lift
    %c = abs(max(Xb) - min(Xb)); % Estimate the chord, based on x
    %c_l = (2 * Circ) / (c); % Calculate the sectional coefficient of lift using Circulation
    cy = 0.0;
    cx = 0.0;
    for i=1:M
        cy = cy - CP(i) * s(i) * cosine(i); % Normal Force Coefficient
        cx = cx - CP(i) * s(i) * -sine(i); % Axial Force Coefficient
    end
    cy = cy/c;
    cx = cx/c;
    c_d  =  (cx) * (cos(angle_of_attack)) + (cy) * (sin(angle_of_attack)); 
    c_l  = -(cx) * (sin(angle_of_attack)) + (cy) * (cos(angle_of_attack));
    
    %% Summation of all Vortices
    Gama(89) = [];
    name_3 = ['The sum of all vortices is equal ', num2str(dot(Gama,s))];
    disp(name_3)
    
    %% Plot the streamlines
    [x_1, y_1] = meshgrid(-0.5:0.002:1.5, -0.5:0.001:0.5);
    psi = zeros(1001, 1001);
    for i = 1:length(x_1)
        for j = 1:length(y_1)
             psi(i, j) = Vinf .* y_1(i, j);
            for k = 1:(length(Gama) - 1)
                psi(i, j) =  psi(i, j) + (Gama(k) ./ (2 * pi)) .* (log((y_1(i, j) - y(k)).^2 + (x_1(i, j) - x(k)).^2));
            end
        end
    end
    figure;
    contour(x_1, y_1, psi, 250); % contour code to plot 250 streamlines
    hold on
    plot(Xb, Yb) % Plot airfoil based on the data
    set(gca,'FontWeight','bold')
    set(gca,'FontSize', 10)
    hold off
    angle_of_attack = (angle_of_attack * 180) / pi; % Convert back to degrees
    name_1 = strcat('Streamlines for \alpha= ', num2str(angle_of_attack));
    title(name_1)
    xlabel ('x-axis')
    ylabel ('y-axis')
    grid on
    axis image
    cleanfigure;
    
    %% Plot the coefficient of pressure Graph for the airfoil
    figure;
    plot(x./c, CP)
    set(gca,'FontWeight','bold')
    set(gca,'FontSize', 10)
    name_2 = strcat('C_{p} vs x/c for \alpha= ', num2str(angle_of_attack));
    title(name_2)
    xlabel('x/c')
    ylabel('C_{p}')
    grid on
    cleanfigure;
    
end