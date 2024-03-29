%% Assignment No. 1

clc;
clear all;
close all;

set(0,'DefaultTextInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'DefaultLineLineWidth', 1.2);

%% Part 1 - Rotations.
% (a) The required function is attached below in the "Functions" Section.
% The function's name is: eul2rotmat(phi, theta, psi)

% (b) Calculation of the rotation matrix from Body to Global
psi_b = pi/7;
theta_b = pi/5;
phi_b = pi/4;
R_B2G = eul2rotmat(phi_b, theta_b, psi_b);
disp('The rotation matrix from Body to Global of the given input is:')
disp(R_B2G)

% (c) The required function is attached below in the "Functions" Section.
% The function's name is: rotmat2eul(R)


% (d) Calculation of the Euler angles in degrees from the given Rotation
% Matrix
R_d = [0.813797681, -0.440969611, 0.378522306;...
     0.46984631, 0.882564119, 0.0180283112;...
    -0.342020143, 0.163175911,  0.925416578];
[phi_d, theta_d, psi_d] = rotmat2eul(R_d); 
% Conversion to deg
phi_d = phi_d*(180/pi);
theta_d = theta_d*(180/pi);
psi_d = psi_d*(180/pi);
disp(['The Euler angles from the given rotation matrix are:'])


% Calculation othe Euler Angles from a given Rotation Matrix
R=[0.813797681, -0.440969611,  0.378522306;
    0.46984631, 0.882564119, 0.0180283112;
     -0.342020143, 0.163175911, 0.925416578];
disp('From the rotation matrix:');
disp(R);
[theta, phi, psi] = rotmat2eul(R);

disp('theta:');
disp(rad2deg(theta));
disp('phi:');
disp(rad2deg(phi));
disp('psi:');
disp(rad2deg(psi));

%% Part 2 - 3D Rigid Transformation.
% Given Rotation Matrix and translation vector in G frame
R_G2C=[0.5363, -0.844, 0;
        0.844, 0.5363, 0;
         0, 0, 1];
t_C2G = transpose([-451.2459, 257.0322, 400]);

T_G2C = rotmatandtranslationvector2posetransformation(R_G2C,R_G2C * t_C2G);

disp('The Pose Transformation Matrix is:');
disp(T_G2C);

% Now Calculate a point in 3D space in global frame to camera frame

l_G = transpose([450, 400, 50]);

l_G_Aug = [l_G; 1];

l_C_Aug = T_G2C * l_G_Aug;

disp('The 3D point in the camera frame is:');
disp(l_C_Aug(1:3));

%% Part 3 - Pose Composition
% (a) We want to find the pose transformation matrices

t_command = transpose([1,0]);
t_actual = transpose([1.01,0]);

degree_command = 0; %degrees
degree_actual = 1; %degrees

R_command = calculateRotationMatrixFromDeg(degree_command);
R_actual = calculateRotationMatrixFromDeg(degree_actual);
 
T_command_0 = poseTransformationFromAngleAndTranslationVector(degree_command,t_command); % T_c(0 -> 1)
T_command = T_command_0;
T_actual_0 = poseTransformationFromAngleAndTranslationVector(degree_actual,t_actual);% T_a(0 -> 1)
T_actual = T_actual_0;

x_0 = transpose([0,0]);
x_0_Aug = createAugmentedVector(x_0);

x_command_Aug = x_0_Aug;
x_actual_Aug = x_0_Aug;

% (b) we shall now simulate the movement for 10 steps
numberOfTimeSteps = 10;

for i=1:numberOfTimeSteps
   x_command_Aug = [x_command_Aug, createAugmentedVector(getTranslationInOriginFrame(T_command))];
   T_command = T_command * T_command_0;
    
   x_actual_Aug = [x_actual_Aug, createAugmentedVector(getTranslationInOriginFrame(T_actual))];
   T_actual = T_actual * T_actual_0;
   
end

x_command = x_command_Aug;
x_command(3,:) = [];
x_actual = x_actual_Aug;
x_actual(3,:) = [];

plotPaths(x_command, x_actual, 'robotPlot', 900);

[errorX, errorY, totalError, errorAngle] = calculateErrorAndAngle(x_command, x_actual);

disp('The error in X is:');
disp(errorX);
disp('The error in Y is:');
disp(errorY);
disp('The total error  is:');
disp(totalError);
disp('The angle error is:');
disp(errorAngle);
%% Functions
% (1.a) Euler angles (in radians) to rotation matrix.
function rot_mat = eul2rotmat(phi, theta, psi) 
    % Declaration of rotation matrices about the principal axes
    R_x = [1, 0, 0; 0, cos(phi), sin(phi); 0, -sin(phi), cos(phi)];
    R_y = [cos(theta), 0, -sin(theta); 0, 1, 0; sin(theta), 0, cos(theta)];
    R_z = [cos(psi), sin(psi), 0; -sin(psi), cos(psi), 0; 0, 0, 1];

    % Calculation of the corresponding rotation matrix by roll-pitch-yaw
    % order from Body to Global
    rot_mat = R_z*R_y*R_x;
end

% (1.c) Rotation matrix to Euler angles (in radians).
% returns an array of angles [phi, theta, psi]
function [theta, phi, psi] = rotmat2eul(R) 
    % From the rotation matrix obtained assuming roll-pitch-yaw order we
    % determined the Euler angles
    theta = asin(R(3,1));
    phi = -atan(R(3,2)/R(3,3));
    psi = -atan(R(2,1)/R(1,1));
end

% (2) Rotation Matrix and Translation Vector to Pose Transformation Matrix
function T = rotmatandtranslationvector2posetransformation(R,t)
    T=[R, t;
        zeros(1,3), 1];
end

% (3.a) Calculate ground vehicle Pose Transformation matrix from angle and 
% translation vector

function T = poseTransformationFromAngleAndTranslationVector(angleInDeg,t)
    R = calculateRotationMatrixFromDeg(angleInDeg);
    T = [R, t;
        zeros(1,2), 1];
end
% (3.a) Calculate ground vehicle Rotation Matrix R from angle
function R = calculateRotationMatrixFromDeg(angleInDeg)
    angleInRad = deg2rad(angleInDeg);
    R = [cos(angleInRad), sin(angleInRad);
         -sin(angleInRad), cos(angleInRad)];
end

% (3.b) create augmented vector(v)
function vBar = createAugmentedVector(v)
    vBar = [v;1];
end 
% (3.b) Calculate Translation in previous(could be origin) frame from pose transformation

function tInOriginFrame = getTranslationInOriginFrame(T)
    R_new = transpose(T(1:2,1:2));
    tInOriginFrame = R_new * T(1:2,3);
end
% (3.b) Function to Plot the data and make jpg
function plotPaths(data1, data2, fileName, resolution)
    % Check if the input matrices have the same number of columns
    assert(size(data1, 2) == size(data2, 2), 'Input matrices must have the same number of columns');

    % Create a new figure
    figure;

    % Plot the first set of data
    plot(data1(1, :), data1(2, :), 'o-', 'LineWidth', 2, 'DisplayName', 'Command');
    hold on;

    % Plot the second set of data
    plot(data2(1, :), data2(2, :), 's-', 'LineWidth', 2, 'DisplayName', 'Actual');

    % Add a grid
    grid on;

    % Set axis limits to make the plot tight
    axis tight;

    % Add title
    title('Plot of Ground Movement of the Autonomous Robot');

    % Add labels for x and y axes
    xlabel('X-axis [m]');
    ylabel('Y-axis [m]');

    % Add legend
    legend('Location', 'best');

    % Hold off to stop overlaying subsequent plots on the same figure
    hold off;

    % Save the plot as a high-resolution JPG file
    print(gcf, fileName, '-djpeg', ['-r' num2str(resolution)]);
    disp(['Plot saved as ', fileName, ' with resolution ', num2str(resolution), ' DPI']);
end


% (3.b) Function to calculate the error
function [errorX, errorY, totalError, errorAngle] = calculateErrorAndAngle(data1, data2)
    % Check if the input matrices have the same number of columns
    assert(size(data1, 2) == size(data2, 2), 'Input matrices must have the same number of columns');

    % Extract the last x and y values
    lastX1 = data1(1, end);
    lastY1 = data1(2, end);
    
    lastX2 = data2(1, end);
    lastY2 = data2(2, end);

    % Calculate the error in x and y
    errorX = abs(lastX1 - lastX2);
    errorY = abs(lastY1 - lastY2);

    % Calculate the total error
    totalError = sqrt(errorX^2 + errorY^2);

    % Calculate the angle (in degrees) between the vectors
    errorAngle = atan2(lastY2 - data2(2, end-1), lastX2 - data2(1, end-1));
end


    
