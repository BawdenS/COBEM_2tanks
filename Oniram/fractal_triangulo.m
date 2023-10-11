clear
close all
clc
k = 10^6;
base = 10;
triangulo_p1 = [0,0];
triangulo_p2 = [base,0];
triangulo_p3 = [base/2,base*cos(30*pi/180)];
ponto= zeros(k,2);
ponto(1,:) = triangulo_p1;
ponto(2,:) = triangulo_p2;
ponto(3,:) = triangulo_p3;
test = rand();
ponto(4,:) = [5,0];

for i = 5:1:k
   evento = randi(3);
   switch evento
       case 1
           ponto(i,1) = (ponto(i-1,1) +ponto(1,1))/2;
           ponto(i,2) = (ponto(i-1,2) +ponto(1,2))/2; 
       case 2
            ponto(i,1) = (ponto(i-1,1) +ponto(2,1))/2;
            ponto(i,2) = (ponto(i-1,2) +ponto(2,2))/2;
       case 3
           ponto(i,1) = (ponto(i-1,1) +ponto(3,1))/2;
            ponto(i,2) = (ponto(i-1,2) +ponto(3,2))/2;
   end
end

plot(ponto(:,1),ponto(:,2),'.k')