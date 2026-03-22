%% ========================================================================
%  Script Name : generate_cart_pendulum_params.m
%  Description : Computes physical constants, system parameters and 
%                state-space matrices for the cart-pendulum system. 
%                Saves output to cart_pendulum_params.mat
%  Author      : Louis David
%  Date        : 2026-03-21
%  Version     : 1.1
%  Usage       : Run before any simulation
%  Dependencies: -
%  Output      : cart_pendulum_params.mat
%  Notes       : Linearised around upright equilibrium theta=0
% =========================================================================

clear ; clc ;

%% ── Physical parameters ─────────────────────────────────────────────────
g = 9.81; % gravitational acceleration (m/s^2)

%% ── Problem's parameters ────────────────────────────────────────────────
M = 0.5; % cart's mass (kg)
m = 0.2; % rod's mass (kg)
l = 0.3; % length from the pivot to the pendulum's COG (m)

I = (1/12) * m * l^2; % rod's moment of inertia about the COG (kg/m^2)

b = 0.1; % coefficient of friction for cart (N.s/m)
c = 0.1; % pendulum viscous friction coefficient (N.m.s/rad)

%% ── State-Space matrices ────────────────────────────────────────────────
den = (M+m)*I + m*M*l^2;

A = [ 0 1                0                0           ;
      0 -(I+m*l^2)*b/den -m^2*l^2*g/den   m*l*c/den   ;
      0 0                0                1           ;
      0 m*l*b/den        (M*m)*m*l*g/den  -(M+m)*c/den];
B= [0; (I+m*l^2)/den; 0; -m*l/den];
C = [1 0 0 0;
     0 0 1 0];
D = zeros(2,1);

clear('den');

system = ss(A,B,C,D);


%% ── Saving parameters ───────────────────────────────────────────────────
save('cart_pendulum_params','g','M','m','l','I', ...
     'b','c','A','B','C','D','system');

load('cart_pendulum_params.mat')
