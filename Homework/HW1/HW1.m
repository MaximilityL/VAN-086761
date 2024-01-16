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
psi = pi/7;
theta = pi/5;
phi = pi/4;
R_B2G = eul2rotmat(phi, theta, psi);
disp('The rotation matrix from Body to Global of the given input is:')
disp(R_B2G)

% (c) The required function is attached below in the "Functions" Section.






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
function eul_angles = rotmat2eul(R) 
    % Using Euler's rotation theorem, first we find the rotation angle and
    % axis of rotation of the rotation matrix input.
    theta = acos(0.5*(trace(R)-1));
    n_vec = 1/(2*sin(theta))*[R(3,2)-R(2,3); R(1,3)-R(3,1), R(2,1)-R(1,2)];
    eul_angles = n_vec;
end


    
