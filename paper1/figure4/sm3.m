%mkdir pics
function sm3
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
disp=disp2;
size(disp)
% figure(81+it)
%     surf(x,y,disp,'FaceColor','interp',...
%    'EdgeColor','none',...
%    'FaceLighting','gouraud')
%     axis equal
%     zlim([0 1])
%     campos=[-23.5875  -30.7398   22.2436];
%     %name=strcat('fig',int2str(it),'.png')
%     caxis([-0 0.5])
figure(20)
clf
size(x);
size(y);
size(disp);

contour(x,y,disp,20)
title(filename)
axis equal
for it=0:nt

filename=strcat('disp',int2str(it),'.out');
disp=reshape(dlmread(filename),ny,nx);

figure(21+it)
clf
contour(x,y,disp,20)
title(filename)
hold on
axis equal

dampx1=dw/nx*(x2-x1)+x1;
line([dampx1 dampx1],[y1 y2])

dampx1=(nx-dw)/nx*(x2-x1)+x1;
line([dampx1 dampx1],[y1 y2])

dampy1=dwy/ny*(y2-y1)+y1;
line([x1 x2],[dampy1 dampy1])

dampy1=(ny-dw)/ny*(y2-y1)+y1;
xlim([x1 x2])
ylim([y1 y2])
 set(gcf,'PaperPositionMode', 'auto')
 imagename=strcat('disp',int2str(it));
    print(imagename,'-depsc','-opengl','-r300')
for i=1:numsensors
    plot(sensorsx(i),sensory,'r.')
end

end

    filename=strcat('ux.out');
    filename
    data = dlmread(filename);
    size(data)
    data=reshape(data,numsensors,totalsteps);
    %data=data(:,1+dw:nx2+dw);
    size(data);
    size(time);
    size(sensorsx);
    size(data);
    size(time);
%      figure(43)
%      clf
%      plot(time,data(10,:))
%      hold on
     filename=strcat('vy.out');
     data = dlmread(filename);
     data=reshape(data,numsensors,totalsteps);   
%      plot(time,data(10,:))

    %     legend('Ux','Uy')
size(sensorsx);
    size(data);
    size(time);
    figure(1)
    clf
    cd ..
    wiggle(sensorsx,time,data,'2hI')
    ylim([sensorsx(1),sensorsx(numsensors)])
    set(gca,'YDir','Normal')
    title('Vy')
cd data
    filename=strcat('vx.out');
     data = dlmread(filename);
     data=reshape(data,numsensors,totalsteps);   
     
    size(data);
    size(time);
    figure(2)
    clf
cd ..
    wiggle(sensorsx,time,data,'2hI')
    ylim([sensorsx(1),sensorsx(numsensors)])
    set(gca,'YDir','Normal')
    title('Vx')
    end   
    function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp);
end