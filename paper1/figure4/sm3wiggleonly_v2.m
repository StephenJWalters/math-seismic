function figure3
%mkdir pics
clc
t=load('t.out');
nt=load('nimages.out');
numsensors=load('numsensors.out');
sensorsx=load('sensorsx.out');

totalsteps=load('totalsteps.out');

totaltime=nt*t;
time=linspace(0,totaltime,totalsteps);

figure(11);
    clf
    Rect    = [0.2, 0.1, 0.75, 0.75];
    rows=2;
    columns=1;
    fsize=24;
    AxisPos = myPlotPos(columns,rows,Rect);
    
for r = 1:rows
        axes('Position', AxisPos(r , :));
        


     if r==2
         cd data
         filename=strcat('vy.out');
     data = dlmread(filename);
     data=reshape(data,numsensors,totalsteps);   
%      plot(time,data(10,:))

    %     legend('Ux','Uy')
    cd ..
    wiggle(sensorsx,time,data,'2hI')
    ylim([sensorsx(1),sensorsx(numsensors)])
    set(gca,'YDir','Normal','FontSize',20)
%    title('Vy')
            xlabel('time (seconds)','FontSize',20)
            ylabel('position (km)','Rotation',90,'FontSize',20)
            set(gcf,'PaperPositionMode', 'auto')
            
text(0.612,1.1,'V_y','FontSize',20);
    %print('figure3Uy','-depsc','-opengl','-r300')
     end

if r==1
    cd data
    filename=strcat('vx.out');
     data = dlmread(filename);
     data=reshape(data,numsensors,totalsteps);   
     
    size(data);
    size(time);
    
cd ..
    wiggle(sensorsx,time,data,'2hI')
    ylim([sensorsx(1),sensorsx(numsensors)])
    set(gca,'YDir','Normal','FontSize',20,'XTickLabel','')
text(0.612,1.1,'V_x','FontSize',20);
    xlabel('')

            ylabel('position (km)','Rotation',90,'FontSize',20)
            set(gcf,'PaperPositionMode', 'auto')
    %print('figure3Ux','-depsc','-opengl','-r300')

end
end
            set(gcf,'PaperPositionMode', 'auto')
print('figure3VxVy','-depsc','-opengl','-r300')
end


function AxisPos = myPlotPos(nCol, nRow, defPos)
% Position of diagrams - a very light SUBPLOT
% This is the percent offset from the subplot grid of the plot box.
% Formula: Space = a + b * n
% Increase [b] to increase the space between the diagrams.
if nRow < 3
   BplusT = -0.00379 + 0.00102345 * nRow;
   %0.18;
else
   BplusT = -0.00379 + 0.00102345 * nRow;
end
if nCol < 3
   LplusR = 0.18;
else
   LplusR = -0.09 + 0.015 * nCol;
end
nPlot = nRow * nCol;
plots = 0:(nPlot - 1);
row   = (nRow - 1) - fix(plots(:) / nCol);
col   = rem(plots(:), nCol);
col_offset  = defPos(3) * LplusR / (nCol - LplusR);
row_offset  = defPos(4) * BplusT / (nRow - BplusT);
totalwidth  = defPos(3) + col_offset;
totalheight = defPos(4) + row_offset;
width       = totalwidth  / nCol - col_offset;
height      = totalheight / nRow - row_offset;
if width * 2 > totalwidth / nCol
   if height * 2 > totalheight / nRow
      AxisPos = [(defPos(1) + col * totalwidth / nCol), ...
            (defPos(2) + row * totalheight / nRow), ...
            width(ones(nPlot, 1), 1), ...
            height(ones(nPlot, 1), 1)];
   else
       AxisPos = [(defPos(1) + col * totalwidth / nCol), ...
            (defPos(2) + row * defPos(4) / nRow), ...
            width(ones(nPlot, 1), 1), ...
            (0.7 * defPos(ones(nPlot, 1), 4) / nRow)];
   end
else
   if height * 2 <= totalheight / nRow
      AxisPos = [(defPos(1) + col * defPos(3) / nCol), ...
            (defPos(2) + row * defPos(4) / nRow), ...
            (0.7 * defPos(ones(nPlot, 1), 3) / nCol), ...
            (0.7 * defPos(ones(nPlot, 1), 4) / nRow)];
   else
      AxisPos = [(defPos(1) + col * defPos(3) / nCol), ...
            (defPos(2) + row * totalheight / nRow), ...
            (0.7 * defPos(ones(nPlot, 1), 3) / nCol), ...
            height(ones(nPlot, 1), 1)];
    end
end
end