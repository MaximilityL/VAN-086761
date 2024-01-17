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

% (b) Calculation of the rotation matrix from Body to Global
psi_b = pi/7;
theta_b = pi/5;
phi_b = pi/4;
R_B2G = eul2rotmat(phi_b, theta_b, psi_b);
disp('The rotation matrix from Body to Global of the given input is:')
disp(R_B2G)

% (c) The required function is attached below in the "Functions" Section.

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



%% Functions
% (a) Euler angles (in radians) to rotation matrix.
function rot_mat = eul2rotmat(phi, theta, psi) 
    % Declaration of rotation matrices about the principal axes
    R_x = [1, 0, 0; 0, cos(phi), sin(phi); 0, -sin(phi), cos(phi)];
    R_y = [cos(theta), 0, -sin(theta); 0, 1, 0; sin(theta), 0, cos(theta)];
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


    
