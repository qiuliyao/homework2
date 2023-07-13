%% assignment 1 2.5.10
%%%%  2.5.10   U0=sin(2*pi*x);------------------------
clear all;
close all;
clc;
u=1;
Nx=100;
x=linspace(0,u,Nx);
dx=u/(Nx-1);
dt=dx*0.5; %
r=dt/(4*dx);

%  T(x,t=0)
U0=sin(2*pi*x);
% U0=0.4+sin(2*pi*x);

U1=zeros(size(U0)); %create an empty matrix for n+1 time layer
U2=zeros(size(U0));
xi=2:Nx-1;  % innner grid points

% FTCS get U1 for CTCS
U1(xi)=U0(xi)+U0(xi).*dt./dx.*(  U0(xi+1)-U0(xi-1)  );

drawnow;
set(gca,'xlim',[0 u],'ylim',[0 max(U0)+0.5])
count=0;
while count<60
    % while count<135
    count=count+1;

    %  central time central space sheme
    U2(xi)=U0(xi)+r*(  (U1(xi-1)+U1(xi)).^2-(U1(xi+1)+U1(xi)).^2  );
    U0=U1;
    U1=U2;
    % Neumann boundary condition: T_x(1)=0, T_x(end)=0
    U0(1)=U0(2);
    U0(end)=U0(end-1);
    U1(1)=U1(2);
    U1(end)=U1(end-1);
    U2(1)=U2(2);
    U2(end)=U2(end-1);

    mi=0.5*sum(U2.^2);
    m(count)=mi;
    if mod(count,5)==0
        plot(x,U2)
        xlim([0 u]);
        ylim([-(max(U0)+0.5) max(U0)+0.5]);
        title(['2.5.10a Time:' num2str(count*dt),'s']);
        xlabel('x')
        ylabel('T')
        pause(0.05);
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if count==5
            imwrite(imind,cm,'assignment1_Bugger’s_eq_2_5_10_a.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment1_Bugger’s_eq_2_5_10_a.gif','gif','WriteMode','append');
        end
    end 
end
figure
plot(m)
xlim([1 60])
% xlim([1 135])
ylim([20 60])
xlabel('Time')
ylabel('KE')
title(['2.5.10 initial condition：u_j^0=sin(2πx)']);

%%%%------2.5.10b   U0=0.5+sin(2*pi*x);-------------------------------------
clear all;
close all;
clc;
u=1;
Nx=100;
x=linspace(0,u,Nx);
dx=u/(Nx-1);
dt=dx*0.5; % 
r=dt/(4*dx);

% U0=sin(2*pi*x);
U0=0.5+sin(2*pi*x);

U1=zeros(size(U0)); %create an empty matrix for n+1 time layer
U2=zeros(size(U0));
xi=2:Nx-1;  % innner grid points

% FTCS get U1 for CTCS
U1(xi)=U0(xi)+U0(xi).*dt./dx.*(  U0(xi+1)-U0(xi-1)  );

drawnow;
set(gca,'xlim',[0 u],'ylim',[0 max(U0)+0.5])
count=0;
while count<135
    count=count+1;

    %  central time central space sheme
    U2(xi)=U0(xi)+r*(  (U1(xi-1)+U1(xi)).^2-(U1(xi+1)+U1(xi)).^2  );
    U0=U1;
    U1=U2;
    % Neumann boundary condition: T_x(1)=0, T_x(end)=0
    U0(1)=U0(2);
    U0(end)=U0(end-1);
    U1(1)=U1(2);
    U1(end)=U1(end-1);
    U2(1)=U2(2);
    U2(end)=U2(end-1);

    mi=0.5*sum(U2.^2);
    m(count)=mi;
    if mod(count,5)==0
        plot(x,U2)
        xlim([0 u]);
        ylim([-(max(U0)+0.5) max(U0)+0.5]);
        title(['2.5.10b Time:' num2str(count*dt),'s']);
        xlabel('x')
        ylabel('T')
        pause(0.05);
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if count==5
            imwrite(imind,cm,'assignment1_Bugger’s_eq_2_5_10_b.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment1_Bugger’s_eq_2_5_10_b.gif','gif','WriteMode','append');
        end
    end 
end
figure
plot(m)
xlim([1 60])
% xlim([1 135])
ylim([20 60])
xlabel('Time')
ylabel('KE')
title(['2.5.10 initial condition：u_j^0=0.5+sin(2πx)']);

%%---------------2.5.11 a-------------------------------------------
clear all;
close all;
clc;
u=1;
Nx=100;
x=linspace(0,u,Nx);
dx=u/(Nx-1);
dt=dx*0.5; % 
r=dt/(3*dx);

% 请尝试几种不同的初值温度分布 T(x,t=0)
U0=sin(2*pi*x);
% U0=0.5+sin(2*pi*x);

U1=zeros(size(U0)); %create an empty matrix for n+1 time layer
U2=zeros(size(U0));
xi=2:Nx-1;  % innner grid points

% FTCS get U1 for CTCS
U1(xi)=U0(xi)+U0(xi).*dt./dx.*(  U0(xi+1)-U0(xi-1)  );

drawnow;
set(gca,'xlim',[0 u],'ylim',[0 max(U0)+0.5])
count=0;
while count<90
    count=count+1;

    Um0=0.5.*(U0(xi-1)+U1(xi-1));
    Um1=0.5.*(U0(xi)+U1(xi));
    Um2=0.5.*(U0(xi+1)+U1(xi+1));
    U2(xi)=U0(xi)+r*(  Um1.*(Um2-Um0)+(Um2.^2-Um0.^2)  );

    U0=U1;
    U1=U2;
    % Neumann boundary condition: T_x(1)=0, T_x(end)=0
    U0(1)=U0(2);
    U0(end)=U0(end-1);
    U1(1)=U1(2);
    U1(end)=U1(end-1);
    U2(1)=U2(2);
    U2(end)=U2(end-1);


    mi=0.5*sum(U2.^2);%
    m(count)=mi;
    if mod(count,5)==0
        plot(x,U2)
        %         xlim([0 u]);
        %         ylim([-(max(U0)+0.5) max(U0)+0.5]);
        title(['2.5.11a Time:' num2str(count*dt),'s']);
        xlabel('x')
        ylabel('T')
        pause(0.05);
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if count==5
            imwrite(imind,cm,'assignment1_Bugger’s_eq_2_5_11_a.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment1_Bugger’s_eq_2_5_11_a.gif','gif','WriteMode','append');
        end
    end 
end

figure
plot(m)
% xlim([1 600])
% xlim([1 135])
% ylim([20 60])
xlabel('Time')
ylabel('总动能')
title(['初始条件：u_j^0=sin(2πx)']);

%%%%----------------2.5.11 b-------------------------------------
clear all;
close all;
clc;

u=1;
Nx=100;
x=linspace(0,u,Nx);
dx=u/(Nx-1);
dt=dx*0.5; % 1可改
r=dt/(3*dx);


% 请尝试几种不同的初值温度分布 T(x,t=0)
% U0=sin(2*pi*x);
U0=0.5+sin(2*pi*x);

U1=zeros(size(U0)); %create an empty matrix for n+1 time layer
U2=zeros(size(U0));
xi=2:Nx-1;  % innner grid points

% FTCS get U1 for CTCS
U1(xi)=U0(xi)+U0(xi).*dt./dx.*(  U0(xi+1)-U0(xi-1)  );


drawnow;
set(gca,'xlim',[0 u],'ylim',[0 max(U0)+0.5])
count=0;
while count<90
    count=count+1;

    Um0=0.5.*(U0(xi-1)+U1(xi-1));
    Um1=0.5.*(U0(xi)+U1(xi));
    Um2=0.5.*(U0(xi+1)+U1(xi+1));
    U2(xi)=U0(xi)+r*(  Um1.*(Um2-Um0)+(Um2.^2-Um0.^2)  );

    U0=U1;
    U1=U2;
    % Neumann boundary condition: T_x(1)=0, T_x(end)=0
    U0(1)=U0(2);
    U0(end)=U0(end-1);
    U1(1)=U1(2);
    U1(end)=U1(end-1);
    U2(1)=U2(2);
    U2(end)=U2(end-1);


    mi=0.5*sum(U2.^2);
    m(count)=mi;
    if mod(count,5)==0
        plot(x,U2)
        %         xlim([0 u]);
        %         ylim([-(max(U0)+0.5) max(U0)+0.5]);
        title(['1维Bugger’s非线性平流方程2.5.11b Time:' num2str(count*dt),'s']);
        xlabel('x')
        ylabel('T')
        pause(0.05);
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if count==5
            imwrite(imind,cm,'assignment2_1_Bugger’s_eq_2_5_11_b.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment2_1_Bugger’s_eq_2_5_11_b.gif','gif','WriteMode','append');
        end
    end %画图部分
end

figure
plot(m)
% xlim([1 600])
% xlim([1 135])
% ylim([20 60])
xlabel('Time')
ylabel('总动能')
title(['初始条件：u_j^0=0.5+sin(2πx)']);

