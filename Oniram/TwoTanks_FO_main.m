% Controle por Trajetória para o modelo SISO Não-Linear
% de um Sistema de 2 Tanques Diferencialmente Plano. 

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
  
  % Agrupamento dos parametros (d1H1)
    num1 = Qe;
    den1 = pho*Area1*(Fnom^3);
    a1 = num1/den1;
    
    num2 = Area2*sqrt(2*g);
    den2 = Area1;
    a2 = num2/den2;
    
  % Agrupamento dos parametros (d1H2)
    num3 = Area2*sqrt(2*g)*Hmax^2;
    den3 = pi*Rmax^2;
    a3 = num3/den3;
    
    num4 = Area4*sqrt(2*g)*Hmax^2;
    den4 = pi*Rmax^2; 
    a4 = num4/den4;
    
% Condições Iniciais do Sistema
  h10 = 1; h20 = 1;
      
% Tempo de Simulação
  t0 = 0; dt = 0.01; tfinal = 25;
  t = t0:dt:tfinal; st = size(t,2);
  
% Planejamento de Trajetória - Saída Plana (Z)

  % Trajetória Constante 
    Zd = 0.5 + 0*t; d1Zd = 0*t; d2Zd = 0*t;
  
  % Trajetória Cosseinodal
    % ro = 0.25; w = 0.50;    
    % Zd = 1.5 - ro*cos(w*t); d1Zd = ro*w*sin(w*t); d2Zd = ro*(w^2)*cos(w*t);
  
% Planejamento de Trajetória - Variáveis do Sistema
  H2d = Zd;
  
  H1d = (Zd.*(a4 + Zd.^(3/2).*d1Zd).^2)/a3^2;
     
  Ud = ((2*Zd.^(5/2)*a4 + 2*Zd.^4.*d1Zd)/(a1*a3^2)).*d2Zd + (a4^2*d1Zd + 4*Zd.^3.*d1Zd.^3 + 5*Zd.^(3/2)*a4.*d1Zd.^2 + a2*a3^2*((Zd.*(a4 + Zd.^(3/2).*d1Zd).^2)/a3^2).^(1/2))/(a1*a3^2);
  
  Fd = Ud.^(1/3);
  
% Parâmetros do Controlador
  p1 = 0.75; P1 = poly(-[p1 p1]); K1 = P1(1,2); K0 = P1(1,3);
  
% Executando a simulação e apresentando os resultados
  sim('TwoTanks_FO_R2018a'); plotResultados_NL;

