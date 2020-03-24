clear
clc
close all

% Useful Visualisations
% A = dlmread("build/A_matrix.txt");

% Reading velocities and Stream-Function
u = dlmread("build/u_velocity.txt");
v = dlmread("build/v_velocity.txt");
psi = dlmread("build/psi_final.txt");
omega = dlmread("build/omega_final.txt");

figure
psi_contour = contourf(flipud(psi))
figure
omega_contour = contourf(flipud(omega), "ShowText", "on")
figure
u_contour = contourf(flipud(u))
figure
v_contour = contourf(flipud(v))

figure
plot(u(ceil(size(u,1)/2),:))
figure
plot(v(ceil(size(v,2)/2),:))