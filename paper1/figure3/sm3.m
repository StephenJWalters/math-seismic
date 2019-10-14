%mkdir pics
function sm3
clc
dw=load('dampwidth.out')
x2=load('x.out');
[nx,~]=size(x2)
y2=load('y.out');
[ny,~]=size(y2)
t=load('t.out');
nt=load('nimages.out');
sensory=load('sensory.out');
numsensors=load('numsensors.out');
sensorsx=load('sensorsx.out');
size(sensorsx);
%dw=0;
nx2=nx-dw*2;
ny2=ny-dw*2;
x=x2(1+dw:nx2+dw);
y=y2(1:ny2+dw);

totalsteps=load('totalsteps.out');

totaltime=nt*t;
time=linspace(0,totaltime,totalsteps);

cd data

for it=0:nt

filename=strcat('disp',int2str(it),'.out');
disp2=reshape(dlmread(filename),ny,nx);
size(disp2)
disp=disp2;
size(disp)
figure(81+it)
    surf(x2,y2,disp,'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','gouraud')
    axis equal
    zlim([0 1])
    campos=[-23.5875  -30.7398   22.2436];
    %name=strcat('fig',int2str(it),'.png')
    caxis([-0 0.5])
figure(21+it)
clf
size(x);
size(y);
size(disp);

contour(x2,y2,disp)
hold on
x1=dw/nx*2*pi-pi;
line([x1 x1],[-pi pi])

x1=(nx-dw)/nx*2*pi-pi;
line([x1 x1],[-pi pi])

y1=dw/ny*2*pi-pi;
line([-pi pi],[y1 y1])

y1=(ny-dw)/ny*2*pi-pi;
line([-pi pi],[y1 y1])

line([pi -pi],[-pi pi],'color','r')


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
     filename=strcat('uy.out');
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
    ylim([-2,2])
    set(gca,'YDir','Normal')
    title('Uy')
cd data
    filename=strcat('ux.out');
     data = dlmread(filename);
     data=reshape(data,numsensors,totalsteps);   
     
    size(data);
    size(time);
    figure(2)
    clf
cd ..
    wiggle(sensorsx,time,data,'2hI')
    ylim([-2,2])
    set(gca,'YDir','Normal')
    title('Ux')
    end   