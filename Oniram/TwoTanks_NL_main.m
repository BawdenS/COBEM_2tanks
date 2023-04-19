% Análise de Malha Aberta para o modelo SISO Não-Linear de um Sistema de 2 Tanques

% Nome: STORMS Group
% Data: Abril/2023

% Pré-Comandos
  close all; clc; clear all;

% Parâmetros do Sistema Não-Linear   
  pho = 998.23; % kg/(m^3)
  Qe = 2.2183/3; % kg/s
  Fnom = 25; % Hz
  g = 9.81; % m/(s^3)
  Area1 = 3.1416; % m^2
  Area2 = 0.0113; % m^2
  Area4 = 0.0314; % m^2
  Rmax = 1; % m
  Hmax = 3; % m
    
% Condições Iniciais do Sistema
  h10 = 1; h20 = 1;
      
% Tempo de Simulação
  t0 = 0; dt = 0.01; tfinal = 2000;
  t = t0:dt:tfinal; st = size(t,2);
  
% Controle em Malha Aberta (F)
  F_rp = 127.5377 + 0*t;
  
  H1_rp = (F_rp.^6*Qe^2)/(2*Area2^2*Fnom^6*g*pho^2);
  
  H2_rp = (Area2^2*H1_rp)/Area4^2;
  
% Executando a simulação
  open('TwoTanks_NL_2018a'); sim('TwoTanks_NL_2018a');
  