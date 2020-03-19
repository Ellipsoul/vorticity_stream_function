clear
clc
close all

% Useful Visualisations
A = dlmread("build/A_matrix.txt");

% Reading velocities and Stream-Function
u = dlmread("build/u_velocity.txt");
v = dlmread("build/v_velocity.txt");
psi = dlmread("build/stream_matrix_new.txt");

contourf(flipud(psi));