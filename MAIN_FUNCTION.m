clc
clear
disp('****************************************************************')
disp('***            NASA/LANGLEY LS(1)-0413MOD AIRFOIL            ***')
disp('****************************************************************')
fprintf ('\n')
disp('Coefficients of Lift/Drag values for different angles of attack:')
disp('----------------------------------------------------------------')

%% Initial declarations
% Define all constants
Vinf = 10; % Freestream velocity
c = 1.0; % Chord length 

% Call function to get X, Y locations for the airfoil
data = load('NASA_LS(1)_0413MOD_DATA_POINTS.dat');
Xb = data(:, 1);
Yb = flip(data(:, 2));

% Plot the airfoil geometry based on the data
figure;
plot(Xb, Yb) % Plot airfoil based on the data
xlabel('Xb')
ylabel('Yb')
title('NASA/LANGLEY LS(1)-0413MOD AIRFOIL - Real Shape')
axis image

%% Determine cl/cd with varying angles of attack
angle_of_attack = [0 3 6 9 12]; % 5 Different Angles of Attack
for i = 1:length(angle_of_attack)
    print = ['At "', num2str(angle_of_attack(i)), ' deg":'];
    disp(print)
    [c_l, c_d] = VORTEX_PANEL_FUNCTION(Xb, Yb, Vinf, c, angle_of_attack(i)); % Call Vortex Panel Method function
    CL(i) = c_l; % Store this value for c_l
    CD(i) = c_d; % Store this value for c_d
    printLift = ['C_L = ', num2str(CL(i), 6)];
    printDrag = ['C_D = ', num2str(CD(i), 6)];
    disp(printLift)
    disp(printDrag)
    fprintf('\n')
end

%% Plot C_L vs alpha for varying angles of attack
figure;
plot(angle_of_attack, CL, 'Linewidth', 2)
hold on
plot(angle_of_attack, CL, 'ko', 'Linewidth', 2)
title('NASA/LANGLEY LS(1)-0413MOD AIRFOIL: C_{L} vs \alpha')
xlabel('\alpha')
ylabel('C_{L}')
cl_values = {num2str(CL(1)), num2str(CL(2)), num2str(CL(3)), num2str(CL(4)), num2str(CL(5))};
text(angle_of_attack + 0.1, CL - 0.025, cl_values)
grid on
hold off

%% Plot C_D vs alpha for varying angles of attack
figure;
plot(angle_of_attack, CD, 'Linewidth', 2)
hold on
plot(angle_of_attack, CD, 'ko', 'Linewidth', 2)
title('NASA/LANGLEY LS(1)-0413MOD AIRFOIL: C_{D} vs \alpha')
xlabel('\alpha')
ylabel('C_{D}')
cl_values = {num2str(CD(1)), num2str(CD(2)), num2str(CD(3)), num2str(CD(4)), num2str(CD(5))};
text(angle_of_attack + 0.1, CD - 0.001, cl_values)
grid on
hold off
