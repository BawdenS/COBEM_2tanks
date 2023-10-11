% Visualizaçao do Resultados

clear
clc
% Pré-Comandos
close all; LW =2;flagPrint=1;
temp_GWO = load('GWO_workspace.mat');         
                
temp_GWO.H1;
temp_GWO.H2;
temp_GWO.ff;



        
temp_PSO = load('resultado_semi_ideal.mat');
temp_PSO.h1f;
temp_PSO.h2f;
temp_PSO.ff;

H1d = temp_PSO.h1_G;
H2d = temp_PSO.h2_G;
Fd = temp_PSO.Fd;
tf = temp_PSO.tf;
%     str = '#F58C06';
    str = '#FF0000';
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
  % Vari?veis de Estado (H1 e H2) + Sinal de Control (U)
    a = figure(1);

    subplot(3,1,1); 
    plot(tf,H1d,'b','LineWidth', LW); 
%     hold on; 
%     plot(tf,temp_PSO.h1f,'Color',color,'LineStyle','--','LineWidth', LW);
%     hold on; 
%     plot(tf,temp_GWO.h1f,'black--','LineWidth', LW);
    
    
    
    hold on; 
    plot(tf,temp_PSO.h1f,'black--','LineWidth', LW);
    hold on; 
    plot(tf,temp_GWO.h1f,'Color',color,'LineStyle','--','LineWidth', LW);
    
    
    
    
%     legend({'Reference','Signal PSO','Signal GWO'},'fontsize', 12 )


    ylabel('H1 [m]', 'fontsize', 20);
    ylim([0 1])
    yticks([0 0.25 0.5 0.75 1])
    
    
    
    subplot(3,1,2); 
    plot(tf,H2d,'b','LineWidth', LW); 
    hold on; 
     plot(tf,temp_PSO.h2f,'black--','LineWidth', LW);
    hold on; 
    plot(tf,temp_GWO.h2f,'Color',color,'LineStyle','--','LineWidth', LW);
    legend({'Reference','Signal PSO','Signal GWO'},'fontsize', 12 )
 
    ylabel('H2 [m]', 'fontsize', 20);
    ylim([0 1])
    yticks([0 0.25 0.5 0.75 1])
    subplot(3,1,3); 
    plot(tf,Fd,'b','LineWidth', LW); 
    hold on; 
    plot(tf,temp_PSO.ff,'black--','LineWidth', LW);
    hold on; 
    plot(tf,temp_GWO.ff,'Color',color,'LineStyle','--','LineWidth', LW);
    
    
    
    
    
    
    
%     hold on; 
%     plot(tf,temp_PSO.ff,'Color',color,'LineStyle','--','LineWidth', LW);
%     hold on; 
%     plot(tf,temp_GWO.ff,'black--','LineWidth', LW);
%     legend({'Reference','Signal PSO','Signal GWO'},'fontsize', 12 )
    xlabel('t [s]', 'fontsize', 20); 
    ylabel('F [Hz]', 'fontsize', 20);
%     yticks([0 10 20 30 40 50]);
    ab = 40;
    a.Position = [150 150 1080 640];
    a.Children(1).YAxis.FontSize = 15;
    a.Children(3).YAxis.FontSize = 15;
    a.Children(4).YAxis.FontSize = 15;
    a.Children(1).XAxis.FontSize = 15;
    a.Children(3).XAxis.FontSize = 15;
    a.Children(4).XAxis.FontSize = 15;
    a.Children(1).YAxis.TickValues = [0 10 20 30 40 50];
%     a.Children(4).YAxis.YTick = [0 10 20 30 40 50];
%     ab = 40;
%     a.Position = [150 150 1080 640];
%     a.Children(1).XLabel.FontSize =  ab;
%     a.Children(1).YLabel.FontSize =  ab;
%     a.Children(2).XLabel.FontSize =  ab;
%     a.Children(2).YLabel.FontSize =  ab;
%     a.Children(3).XLabel.FontSize =  ab;
%     a.Children(3).YLabel.FontSize =  ab;
    
%     a.Children(1).FontSize = 10;
%     a.Children(2).FontSize = 10;
%     a.Children(3).FontSize = 10;

% 
%   % Pré-Comandos
    addpath('./Resultados/Export_fig/');
%     close all; LW = 2; flagPrint = 1;
% 
%   % Variáveis de Estado (H1 e H2) + Sinal de Control (U)
%     a = figure(1);
%     subplot(3,1,1); plot(t,H1d,'b','LineWidth', LW); hold on; plot(tf,h1f,'r--','LineWidth', LW);
%     % title('Evolução Temporal do Nível H1', 'fontsize', 14);
%     xlabel('t [s]', 'fontsize', 18); ylabel('H1 [m]', 'fontsize', 18);
%     
%     subplot(3,1,2); plot(t,H2d,'b','LineWidth', LW); hold on; plot(tf,h2f,'r--','LineWidth', LW);  
%     % title('Evolução Temporal do Nível H2', 'fontsize', 14);
%     xlabel('t [s]', 'fontsize', 18); ylabel('H2 [m]', 'fontsize', 18);
% 
%     subplot(3,1,3); plot(t,Fd,'b','LineWidth', LW); hold on; plot(tf,ff,'r--','LineWidth', LW);
%     % title('Evolução Temporal do Sinal de Controle', 'fontsize', 14);
%     xlabel('t [s]', 'fontsize', 18); ylabel('F [Hz]', 'fontsize', 18);
%     
%     ab = 40;
%     a.Position = [150 150 1080 640];
%     a.Children(1).XLabel.FontSize =  ab;
%     a.Children(1).YLabel.FontSize =  ab;
%     a.Children(2).XLabel.FontSize =  ab;
%     a.Children(2).YLabel.FontSize =  ab;
%     a.Children(3).XLabel.FontSize =  ab;
%     a.Children(3).YLabel.FontSize =  ab;
%     
%     a.Children(1).FontSize = 18;
%     a.Children(2).FontSize = 18;
%     a.Children(3).FontSize = 18;
    
    if(flagPrint == 1)
%       maximize(gcf); 
      export_fig('./Resultados/2Tanks', '-eps', '-transparent'); 
      eps2pdf('./Resultados/2Tanks.eps', './Resultados/correcao_pedida.pdf', 1);
    end