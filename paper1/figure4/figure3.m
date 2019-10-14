function figure3
%mkdir pics
clc
x1=load('x1.out')
x2=load('x2.out')
y1=load('y1.out')
y2=load('y2.out')
dw=load('dampwidth.out')
dwy=load('dampwidthy.out')
x=load('x.out');
[nx,~]=size(x)
y=load('y.out');
[ny,~]=size(y)
t=load('t.out');
nt=load('nimages.out');
sensory=load('sensory.out');
numsensors=load('numsensors.out');
sensorsx=load('sensorsx.out');
size(sensorsx);

totalsteps=load('totalsteps.out');

totaltime=nt*t;
time=linspace(0,totaltime,totalsteps);

cd data

filename=strcat('disp00.out');
disp2=reshape(dlmread(filename),ny,nx);
size(disp2)
disp=disp2(dwy+1:ny,dw+1:nx-dw);
size(disp)
x=x(dw+1:nx-dw);
y=y(dwy+1:ny);
x1=x(1);
x2=x(nx-2*dw);
y1=y(1);

figure(20)
clf
size(x);
size(y);
size(disp);

contour(x,y,disp,20)
title(filename)
axis equal

figure(11);
    clf
    Rect    = [0.2, 0.1, 0.75, 0.75];
    rows=4;
    columns=1;
    fsize=24;
    AxisPos = myPlotPos(columns,rows,Rect);
    
for r = 1:rows
        i=r
        axes('Position', AxisPos(i, :));
        filename=strcat('disp',int2str((r-1)*2),'.out');
        disp2=reshape(dlmread(filename),ny,nx);
        disp=disp2(dwy+1:ny,dw+1:nx-dw);
        contour(x,y,disp,20)
        %title(filename)
        hold on
        for s=1:numsensors
          plot(sensorsx(s),sensory,'r.')
        end
        axis equal
        xlim([x1 x2])
        ylim([y1 y2])
        set(gca,'ytick',[-0.3 -0.2 -0.1],'FontSize',16)
        if r<rows
            set(gca,'xtick',[])
        end
        if r==rows
            xlabel('x','FontSize',16)
        end
        ylabel('y','Rotation',0,'Position',[-0.056 -0.18],'FontSize',16)
end

cd ..        
        set(gcf,'PaperPositionMode', 'auto')
    print('figure3','-depsc','-opengl','-r300')
cd data
         filename=strcat('uy.out');
     data = dlmread(filename);
     data=reshape(data,numsensors,totalsteps);   
%      plot(time,data(10,:))

    %     legend('Ux','Uy')
size(sensorsx);
    size(data);
    size(time);
    figure(8)
    clf
    cd ..
    wiggle(sensorsx,time,data,'2hI')
    ylim([sensorsx(1),sensorsx(numsensors)])
    set(gca,'YDir','Normal')
    title('Uy')
            xlabel('time (seconds)','FontSize',20)
            ylabel('sensor location (km)','Rotation',90,'FontSize',20)
            set(gcf,'PaperPositionMode', 'auto')
    print('figure3Uy','-depsc','-opengl','-r300')

cd data
    filename=strcat('ux.out');
     data = dlmread(filename);
     data=reshape(data,numsensors,totalsteps);   
     
    size(data);
    size(time);
    figure(7)
    clf
cd ..
    wiggle(sensorsx,time,data,'2hI')
    ylim([sensorsx(1),sensorsx(numsensors)])
    set(gca,'YDir','Normal')
    title('Ux')
    xlabel('time (seconds)','FontSize',20)
            ylabel('sensor location (km)','Rotation',90,'FontSize',20)
            set(gcf,'PaperPositionMode', 'auto')
    print('figure3Ux','-depsc','-opengl','-r300')


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