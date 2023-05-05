% Análise de Malha Aberta para o modelo SISO Não-Linear de um Sistema de 2 Tanques

% Nome: STORMS Group
% Data: Maio/2023

% Pré-Comandos
  close all; clc; clear all;

% Parâmetros do Sistema Não-Linear   
 
  % Tanque 01
    D1 = (5+3/8)*2.54/100;          % m
    Area1 = pi*D1^2/4;              % m^2
    H1max = 1;                      % m          30*2.54/100
    Dout = 0.05;                    % m
    Aout1 = pi*(Dout^2)/4;          % m^2
  
  % Tanque 02
    R2max = 0.176;                  % m
    H2max = 1;                      % m
    Aout2 = Aout1*3/4;              % m^2
  
  % Outros
    g = 9.81;                       % m/(s^2)
    pho = 998.23;                   % kg/(m^3)
    Fnom = 60;                      % Hz
    Qnom = 16.67;                   % kg/s  -> 60 m3/h
    
  % Agregados
    num1 = Qnom; den1 = pho*Area1*(Fnom); a1 = num1/den1;
    num2 = Aout1*sqrt(2*g); den2 = Area1; a2 = num2/den2;
    num3 = Aout1*sqrt(2*g)*H2max^2; den3 = pi*R2max^2; a3 = num3/den3;
    num4 = Aout2*sqrt(2*g)*H2max^2; den4 = pi*R2max^2; a4 = num4/den4;
  
% Condições Iniciais do Sistema
  h10 = 0.25; h20 = 0.25;
      
% Tempo de Simulação
  t0 = 0; dt = 0.01; tfinal = 100;
  t = t0:dt:tfinal; st = size(t,2);
  
% Controle em Malha Aberta (F)
  F_rp = 20.2963 + 0*t;
  
  H1_rp = (F_rp.^2*a1^2)/a2^2;
  
  H2_rp = (H1_rp.*a3^2)/a4^2; 
    
% Executando a simulação
  sim('TwoTanks_NL_2018a');
  