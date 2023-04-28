% Análise de Malha Aberta para o modelo SISO Não-Linear de um Sistema de 2 Tanques

% Nome: STORMS Group
% Data: Abril/2023

% Pré-Comandos
  close all; clc; clear all;

% Parâmetros do Sistema Não-Linear   
  pho = 998.23; % kg/(m^3)
  Fnom = 30; % Hz
  g = 9.81; % m/(s^3)
  D1 = (5+3/8)*2.54/100; % diametro
  Area1 = pi*D1^2/4; % m^2
  H1max = 30*2.54/100;
  
  Dout = 0.05;
  Aout1 = pi*(Dout^2)/4; %0.0113; % m^2


  
% Eq de Bernoulli
  Qe = pho*Aout1*sqrt(2*g*H1max); % kg/s
  
  
  % Tanque 2
  Rmax = 0.176; % m
  Hmax = 0.7; % m
  Aout2 = Aout1*2/3; % m^2
  
% Condições Iniciais do Sistema
  h10 = 0.01; h20 = 0.01;
      
% Tempo de Simulação
  t0 = 0; dt = 0.01; tfinal = 2000;
  t = t0:dt:tfinal; st = size(t,2);
  
% Controle em Malha Aberta (F)
  F_rp = 30;
  
%   H1_rp = (F_rp.^6*Qe^2)/(2*Aout^2*Fnom^6*g*pho^2);
%   
%   H2_rp = (Aout^2*H1_rp)/Area4^2;
  
% Executando a simulação
  open('TwoTanks_NL_2018a'); 
  Saida = sim('TwoTanks_NL_2018a');
  
  
  plot(Saida.Nivel_h1{1}.Values.Time,Saida.Nivel_h1{1}.Values.Data);
  xlabel("Tempo[s]");
  ylabel("Altura H1[m]");
  figure();
  plot(Saida.Nivel_h2{1}.Values.Time,Saida.Nivel_h2{1}.Values.Data);
  xlabel("Tempo[s]");
  ylabel("Altura H2[m]");