% Pré-Comandos
  clear all; close all; clc

% Definindo variáveis simbólicas
  syms H1 H2 F U;
  syms a1 a2 a3 a4;
  syms d1H1 d1H2 d2H2
  
% Sistema Não-Linear
  Eq1 = d1H1 - (a1*U - a2*sqrt(H1)); % U = F^3
    
  Eq2 = d1H2 - (a3*(sqrt(H1)/(H2^2)) - a4/(H2^1.5));
  
% Parametrizacao Diferencial de H1 e F em funçao de H2 e suas derivadas
  h1 = simplifyFraction(solve(Eq2 == 0, H1))
    
  d1h1 = simplify(transpose(gradient(h1,[H2 d1H2]))*[d1H2; d2H2])
  
  u = simplifyFraction(solve(Eq1 == 0, U));
  u = collect(simplifyFraction(subs(u, [H1 d1H1], [h1 d1h1])), d2H2)
  
  
 

 
         
 