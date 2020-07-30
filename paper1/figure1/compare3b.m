function compare3b
%mkdir pics
clc
dw=load('dampwidth.out')
x2=load('xout.out');
[nx,~]=size(x2)
y2=load('yout.out');
[ny,~]=size(y2)
t=load('tlin.out')
[nt,~]=size(t);
cd data
betax=load('betaxout.out');

betax2=reshape(betax,ny,nx);
betax2(1,:);
size(betax2);
figure(60)
surf(x2,y2,betax2)
betaxx=load('betaxxout.out');
size(betaxx)

betaxx2=reshape(betaxx,ny,nx);
size(betaxx2);
figure(61)
surf(x2,y2,betaxx2)
%dw=0;
nx2=nx-dw*2
ny2=ny-dw*2

x=x2(1+dw:nx2+dw);
y=y2(1+dw:ny2+dw);

filename2=strcat('disp00.out');
            disp2=reshape(dlmread(filename2),ny,nx);
            disp=disp2(1+dw:nx2+dw,1+dw:ny2+dw);

%zlim([-1 1])
%campos=[-23.5875  -30.7398   22.2436];
max(abs(disp(:)))
figure(11);
    clf
    Rect    = [0.1, 0.1, 0.75, 0.75];
    rows=3;
    columns=4;
    fsize=24;
    AxisPos = myPlotPos(columns,rows,Rect);
    
for c = 1:columns
    filename=strcat('umag',int2str(c+1),'.out');
        Umag2=reshape(dlmread(filename),ny,nx);
        Umag=Umag2(1+dw:nx2+dw,1+dw:ny2+dw);

    filename2=strcat('disp',int2str(c),'.out');
        disp2=reshape(dlmread(filename2),ny,nx);
        disp=disp2(1+dw:nx2+dw,1+dw:ny2+dw);
    %if c==1
        maxdisp=max(abs(disp(:)))
    %end
    for r = 1:rows
        i=(r-1)*columns+c;
        axes('Position', AxisPos(i, :));


        if r==1       
            surf(x,y,(Umag)','FaceColor','interp','EdgeColor','none',...
            'FaceLighting','gouraud')
            %axis equal    
            zlim([0 1.25])
            caxis([0 1.25])
            %campos=[-23.5875  -30.7398   22.2436];
            title(strcat('t = ',num2str(t(c+1))),'FontSize',24)
            if c==1
                zlabel('Analytic','FontSize',24,'Rotation',0,'Position',[-6 3 0]);  
            end
            if c==4
                colorbar('ticks',[0.25 0.5 0.75 1],'ticklabels',[0.25 0.5 0.75 1],'fontsize',24);  
            end 
        end
        if r==2
            surf(x,y,(disp)','FaceColor','interp','EdgeColor','none',...
            'FaceLighting','gouraud')
            %axis equal    
            zlim([0 1.25])
            caxis([0 1.25])
            %campos=[-23.5875  -30.7398   22.2436];
            if c==1
                zlabel('Numerical','FontSize',24,'Rotation',0,'Position',[-6 3 0]);  
            end
            if c==4
                colorbar('ticks',[0.25 0.5 0.75 1],'ticklabels',[0.25 0.5 0.75 1],'fontsize',24);  
            end 
        end
        if r==3
            diff=(-Umag+disp);
            difference=max(abs(diff(:)))
            surf(x,y,diff','FaceColor','interp','EdgeColor','none',...
            'FaceLighting','gouraud')
            %axis equal    
            zlim([-0.002 0.002])
            caxis([-0.002 0.002])
            %campos([-23.5875,  -30.7398,   1.2436]);
            if c==1
                zlabel('Difference','FontSize',24,'Rotation',0,'Position',[-6 3 0]);  
            end
            if c==4
                colorbar('ticks',[-0.001 0 0.001],'ticklabels',[-0.001 0 0.001],'fontsize',24);  
            end 
        end

    end
        im=2+4*(c-1)
        titlestring=strcat('t  =',{'  '},int2str(im));
    
end
% hp4 = get(subplot(1,4,4),'Position')
% colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.1  hp4(2)+hp4(3)*2.1])   
cd ..
        
        set(gcf,'PaperPositionMode', 'auto')
    print('fig1','-depsc','-opengl','-r300')
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