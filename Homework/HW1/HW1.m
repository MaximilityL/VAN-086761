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
psi = pi/7;
theta = pi/5;
phi = pi/4;
R_B2G = eul2rotmat(phi, theta, psi);
disp('The rotation matrix from Body to Global of the given input is:');
disp(R_B2G);

% (c) The required function is attached below in the "Functions" Section.
% The function's name is: rotmat2eul(R)

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


%% Functions
% (a) Euler angles (in radians) to rotation matrix.
function rot_mat = eul2rotmat(phi, theta, psi) 
    % Declaration of rotation matrices about the principal axes
    R_x = [1, 0, 0; 0, cos(phi), sin(phi); 0, -sin(phi), cos(phi)];
    R_y = [cos(theta), 0, 0; 0, 1, 0; sin(theta), 0, cos(theta)];
    R_z = [cos(psi), sin(psi), 0; -sin(psi), cos(psi), 0; 0, 0, 1];

    % Calculation of the corresponding rotation matrix by roll-pitch-yaw
    % order from Body to Global
    rot_mat = R_z*R_y*R_x;
end

% (a) Rotation matrix to Euler angles (in radians).
% returns an array of angles [phi, theta, psi]
function [theta, phi, psi] = rotmat2eul(R) 
    % From the rotation matrix obtained assuming roll-pitch-yaw order we
    % determined the Euler angles
    theta = asin(R(3,1));
    phi = -atan(R(3,2)/R(3,3));
    psi = -atan(R(2,1)/R(1,1));
end


    
