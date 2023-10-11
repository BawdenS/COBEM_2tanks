% Controle via PID para o modelo SISO Não-Linear
% de um Sistema de 2 Tanques Diferencialmente Plano. 

% Nome: STORMS Group
% Data: Junho/2023

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
  h10 = 0.10; h20 = 0.10;
      
% Tempo de Simulação
  t0 = 0; dt = 0.01; tfinal = 25*4;
  t = t0:dt:tfinal; st = size(t,2);
  
% Planejamento de Trajetória - Saída Plana (Z)

  % Trajetória Constante 
    % Zd = 0.80 + 0*t; d1Zd = 0*t; d2Zd = 0*t;
    
  % Trajetória Constante com 2 valores (STEP)
    fP1 = 0.30; fP2 = 0.60;
    Zd = 0.40 + 0*t; Zd((round(fP1*st)):end) = 0.60; Zd((round(fP2*st)):end) = 0.50;
    d1Zd = 0*t; d2Zd = 0*t;
      
% Planejamento de Trajetória - Variáveis do Sistema
  H2d = Zd;
  
  H1d = (Zd.*(a4 + Zd.^(3/2).*d1Zd).^2)/a3^2;
     
  Fd = ((2*Zd.^(5/2)*a4 + 2*Zd.^4.*d1Zd)/(a1*a3^2)).*d2Zd + (a4^2*d1Zd + 4*Zd.^3.*d1Zd.^3 + 5*Zd.^(3/2)*a4.*d1Zd.^2 + a2*a3^2*((Zd.*(a4 + Zd.^(3/2).*d1Zd).^2)/a3^2).^(1/2))/(a1*a3^2);

% Parâmetros do Controlador (PID)
  kP = 150;     % [50:250]
  kI = 30;      % [10:50]
  kD = 5;       % [1:10]
  
  % Saturação do Atuador
    satF_sup = 90;
    satF_inf = 10;
  
% Executando a simulação e apresentando os resultados
  sim('TwoTanks_PID_R2018a'); plotResultados_NL;

