temp = 0;
k = 1; %constante de ganho da derivada
for i = 1:1:(size(ff,1)-1)
     temp = temp + k*abs(ff(i) -  ff(i+1));  
end
% disp(temp)

%calcular para h1f
temp_2 = (t*(1.5*abs(H1d'-h1f) + 3*abs(H2d'-h2f))); 
disp(t*(1.5*abs(H1d'-h1f) + 3*abs(H2d'-h2f)))
k = temp_2/(1*temp)  ; %constante de ganho da derivada

% disp(k)
% calcular para h2f


temp = 0;
for i = 1:1:(size(ff,1)-1)
     temp = temp + k*(ff(i) -  ff(i+1));  
end
% disp(temp)
sprintf('O valor do custo Du/Dt é %.15g',temp)
% disp(temp_2)
sprintf('O valor da função custo H1 e H2 é %.15g',temp_2)
temp = temp + temp_2;
% disp(temp)

sprintf('O valor da função custo total é %.15g',temp)

sprintf('O valor do K é %.15g',k) % 290