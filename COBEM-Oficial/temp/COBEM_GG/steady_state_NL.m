% Pré-Comandos
  clear all; close all; clc

% Definindo variáveis simbólicas
  syms H1 H2 F;
  syms a1 a2 a3 a4
    
% Sistema Não-Linear
  d1H1 = a1*F - a2*sqrt(H1);
    
  d1H2 = a3*(sqrt(H1)/(H2^2)) - a4*(1/(H2^1.5));
  
% Dado F em regime permanente, determinando os valores de H1 e H2
  S1 = solve(d1H1 == 0, H1);
  
  S2 = solve(d1H2 == 0, H2);
 
         
 