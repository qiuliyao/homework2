% the finite-difference model with Arakawa 1966 Jaccobian for 2D turbulence simulation
% with double periodic boundary condition, i.e. in both x and y directions,
% Adam -Bashforth 2step explicit scheme
clear all;
close all;
clc;


dx=1; %grid size
N=100; %number of grid in x or y direction
ix=2:N-1;%inner grid points
jy=2:N-1;

coeffi_Jacobian= -1/(12*dx^2);%Arakawa 1966 Eq(45,46)

r_d2x=1/(2*dx);
r_dx2=1/(dx^2);
Nts=8000; %total number of time step if no restart

Ah= 1.e-5;  %diffusion
bta=0; %strength of the beta effect
Amp=1;
% 3 time layers for vorticity field
zeta0= zeros(N,N);   % initial field, n-1 time layer
zeta1= zeros(N,N);  % n time layer
zeta2= zeros(N,N); % n+1 time layer
psi=zeros(N);%open a new matrix for psi, the solution to Poisson eq.
psi_guess=rand(N);

% initialize the inital relative vorticity field
for k=1:N
    for j=1:N
        r2=(j/N-0.35)^2*60+(k/N-0.25)^2*100;
        r3=(j/N-0.65)^2*30+(k/N-0.75)^2*100;
        zeta0(k,j) = exp(    -r2    )-exp(-r3)     ;

        %                 r1=(j/N-0.25)^2*200+(k/N-0.25)^2*300;
        %         r2=(j/N-0.35)^2*300+(k/N-0.75)^2*200;
        %         r3=(j/N-0.65)^2*400+(k/N-0.5)^2*400;
        %         zeta0(k,j) = exp(    -r2    )-exp(-r3)  +exp(-r1)  ;
    end
end
zeta0=Amp*zeta0 ;
[xx,yy]=meshgrid(1:N,1:N);
% zeta0=10*sin(2*pi/N*6*xx+2*pi*10/N*yy)-13*cos(2*pi/N*16*xx-2*pi*10/N*yy)-rand(N);
% zeta0=4*rand(N);

figure;imagesc(zeta0);title('initial vorticity'); colorbar
% zeta_test=(psi(ix+1,jy)+psi(ix-1,jy)+psi(ix,jy+1)+psi(ix,jy-1)-4*psi(ix,jy));
% figure;imagesc(zeta_test);title('test zeta')


% solve for the psi for the initial field using SOR with periodic boundary condition

SOR2d
figure
imagesc(psi);
title('SOR solution of \psi for the initial \zeta field')
colorbar

u=max(max(abs(diff(psi)/dx)));
dt=([dx/u])*0.2;


jac
zeta1=zeta0+J0*dt; %Euler forward in time to predict n+1 time layer
zeta0=zeta1;  SOR2d;  jac;
J1=J0;



tic;
count=0;
ts=0;
figure('position',[0,0,1200,400])
drawnow

T=0;

while( ts<Nts)
    ts=ts+1;
    T=T+dt;

    zeta2=zeta1+dt*(1.5*J1-0.5*J0) ;     % Adam-Bath 2nd order time scheme for the nonlinear term
    zeta2(jy,ix)=zeta2(jy,ix)+........
        (zeta2(jy+1,ix)+zeta2(jy-1,ix)+zeta2(jy,ix+1)+zeta2(jy,ix-1))*Ah*r_dx2*dt;
    %     -......
    %         dt*bta*(psi(jy,ix+1)-psi(jy,ix-1))*0.5;
    % diffusion_periobc



    zeta2=zeta2/(1+4*Ah*dt*r_dx2);




    zeta1=zeta2;

    zeta0=zeta2;   SOR2d;
    jac;        J2=J0;



    J0=J1;
    J1=J2;


    diag_ts(ts,1)=sum(sum(zeta2.^2));

    cputime_per_step=toc;
    if mod(ts,20)==0
        subplot(121);
        contourf(zeta2,20);colorbar;
        caxis([-1 1]*Amp*0.8);
        set(gca,'ydir','reverse')
        title([num2str(ts),' steps, cpu time ',......
            num2str(cputime_per_step,'%3.2f'),'s, model T=',num2str(T,'%3.2f'),'s']);

        subplot(122)
        plot(diag_ts,'bx-');
        xlabel('time step')
        ylabel('enstrophy')
        pause(0.1);
        u=max(max(abs(diff(psi)/dx)));
        dt=min([dx/u])*0.1;
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if ts==200
            imwrite(imind,cm,'assignment2_2_eddysimulate.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment2_2_eddysimulate.gif','gif','WriteMode','append');
        end
    end
end


