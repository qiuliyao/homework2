clear all;close all;clc

f=3.e-3; %Coriolis parameter s^-1
Cd=1.e-5; %Drag coefficiant, non-dimensional
g=9.8; %gravitational acceleration, m/s^2
L=1.e6;% dimension of the domain in meter
Nx=80; %no. of grid in each direction
Ny=Nx;
dx=L/Nx; %grid size, in m
xc=0.5*Nx; %center of the initial perturbation

% initialize the matrix
u=zeros(Nx,Ny); %u(n)
v=zeros(Nx,Ny); %v(n)
eta=zeros(Nx,Ny); %eta(n)
h=zeros(Nx,Ny);
u1=u; %u(n+1)
v1=v; %v(n+1)
eta1=eta; %eta(n+1)
d=u;
e=u;
beta1=u; %coe used in the frictin term
beta2=u; %coe used in the friction term
us=u; %imnmediate step u*
vs=v; %immediate step v*


[X,Y]=meshgrid(1:Nx,1:Ny);
H0=2000; %mean water depth,
h=H0*(1-0.8*exp( (-(X-xc*1.6).^2-(Y-xc*1.6).^2))); %add bathmetry
figure;mesh((1:Nx)*dx*1.e-3,(1:Nx)*dx*1.e-3,-h);colorbar;title('bathmetry, m')
xlabel('x, km')
ylabel('y, km')


dt=min([dx/sqrt(g*max(h(:))),Cd/f^2]); %time step, s
% dt=10; %time step in seconds
ns=max([10/dt, 20]);

%initial surf ssh
eta=0.03*exp( -(X-xc).^2/8^2-(Y-xc).^2/8^2);
H=eta+h;

ju=2:Nx-1; %x-grid (eastward) indices for u and d
ku=1:Ny-1; %y-grid (northward)indicces for u and d
d(ju,ku)=0.5*( H(ju,ku)+H(ju-1,ku) );

jv=1:Nx-1; %x-grid (eastward) indices for v and e
kv=2:Ny-1; %y-grid (northward)indicces for v and e
e(jv,kv)=0.5*( H(jv,kv)+H(jv,kv-1) );

je=jv;   % x-grid (eastward) indices for eta
ke=ku;  % y-grid (northward) indices for eta

ts=sum(eta(:)); %used for check mass conservation;
disp(['initial mass is ',num2str(ts),'!']);



niter=0;T=0;
while niter*dt<(1801.1*10)/3*2
    % ======================================第一步=============================
    % 3.5.1第一式先算水位
    eta1(je,ke)=eta(je,ke)-(d(je+1,ke).*u(je+1,ke)-d(je,ke).*u(je,ke)......
        +e(je,ke+1).*v(je,ke+1)-e(je,ke).*v(je,ke))*dt/dx;

    % update H， d， e， and the coeff beta1 beta2 used in eq. 3.5.1'
    H=eta1+h;
    d(ju,ku)=0.5*( H(ju,ku)+H(ju-1,ku) );
    e(jv,kv)=0.5*( H(jv,kv)+H(jv,kv-1) );
    beta1(ju,ku)=Cd*sqrt( u(ju,ku).^2 +......
        0.25^2*(v(ju,ku+1)+v(ju,ku)+v(ju-1,ku)+v(ju-1,ku+1)).^2 ......
        )./d(ju,ku);
    beta2(jv,kv)=Cd*sqrt( v(jv,kv).^2 +......
        0.25^2*(u(jv,kv)+u(jv+1,kv)+u(jv,kv-1)+u(jv+1,kv-1)).^2 ......
        )./e(jv,kv);

    % % % eq. 3.5.1': 中间时刻u*=us v*=vs,
    us(ju,ku) = u(ju,ku) + f*dt*0.25*( v(ju,ku+1)+v(ju,ku)+v(ju-1,ku)+v(ju-1,ku+1) ) ......
        - g*dt/dx*( eta1(ju,ku)-eta1(ju-1,ku) );
    us(ju,ku) = us(ju,ku)./(1+beta1(ju,ku)*dt);
    vs(jv,kv) = v(jv,kv) - f*dt*0.25*( u(jv,kv)+u(jv+1,kv)+u(jv,kv-1)+u(jv+1,kv-1) ) ......
        - g*dt/dx*( eta1(jv,kv)-eta1(jv,kv-1) );
    vs(jv,kv) = vs(jv,kv)./(1+beta2(jv,kv)*dt);


    %================================= 第二步  ================================
    % Coeffs in the 2nd step, eq.3.5.2'
    lambd=zeros(Nx,Ny);
    lambdr=zeros(Nx,Ny);
    lambdl=zeros(Nx,Ny);
    nu=zeros(Nx,Ny);
    nuu=zeros(Nx,Ny);
    nud=zeros(Nx,Ny);

    lambd(ju,ku)=0.25*(us(ju-1,ku)+2*us(ju,ku)+us(ju+1,ku))*dt/(2*dx);
    nuu(ju,ku)=(vs(ju-1,ku+1)+vs(ju,ku+1))*dt/(4*dx);
    nud(ju,ku)=(vs(ju-1,ku)+vs(ju,ku))*dt/(4*dx);

    nu(jv,kv)=0.25*(vs(jv,kv-1)+2*vs(jv,kv)+vs(jv,kv+1))*dt/(2*dx);
    lambdr(jv,kv)=(us(jv+1,kv)+us(jv+1,kv-1))*dt/(4*dx);
    lambdl(jv,kv)=(us(jv,kv)+us(jv,kv-1))*dt/(4*dx);

    for jj=2:Nx-1
        for kk=1:Ny-1
            if kk>1
                u1(jj,kk)=(1+lambd(jj,kk)+nuu(jj,kk))*us(jj,kk);
            else
                u1(jj,kk)=(1+lambd(jj,kk)+nuu(jj,kk))*us(jj,kk);
            end
            u1(jj,kk)=u1(jj,kk)/(1+lambd(jj,kk)+nud(jj,kk));
        end
    end


    for kk=2:Ny-1
        for jj=1:Nx-1
            if jj>1
                v1(jj,kk)=(1+lambdr(jj,kk)+nu(jj,kk))*vs(jj,kk);
            else
                v1(jj,kk)=(1+lambdr(jj,kk)+nu(jj,kk))*vs(jj,kk);
            end
            v1(jj,kk)=v1(jj,kk)/(1+lambdl(jj,kk)+nu(jj,kk));
        end
    end

    %======================= 第三步 n+1时刻(eta, d, e) -> n+2时刻的eta=============
    eta(je,ke)=eta1(je,ke)-(d(je+1,ke).*u1(je+1,ke)-d(je,ke).*u1(je,ke)+......
        e(je,ke+1).*v1(je,ke+1)-e(je,ke).*v1(je,ke))*dt/dx;

    % update beta1 beta2 again for eq.3.5.1'
    H=eta+h;
    d(ju,ku)=0.5*( H(ju,ku)+H(ju-1,ku) );
    e(jv,kv)=0.5*( H(jv,kv)+H(jv,kv-1) );
    beta1(ju,ku)=Cd*sqrt( u1(ju,ku).^2 +......
        0.25^2*(v1(ju,ku+1)+v1(ju,ku)+v1(ju-1,ku)+v1(ju-1,ku+1)).^2 ......
        )./d(ju,ku);
    beta2(jv,kv)=Cd*sqrt( v1(jv,kv).^2 +......
        0.25^2*(u1(jv,kv)+u1(jv+1,kv)+u1(jv,kv-1)+u1(jv+1,kv-1)).^2 ......
        )./e(jv,kv);

    % % % eq. 3.5.1': 注意是n+1和n+2时刻的中间时刻
    us(ju,ku)=u1(ju,ku)+f*dt*0.25*(v1(ju,ku+1)+v1(ju,ku)+v1(ju-1,ku)+v1(ju-1,ku+1))......
        -g*dt/dx*( eta(ju,ku)-eta(ju-1,ku) );
    us(ju,ku)=us(ju,ku)./(1+beta1(ju,ku)*dt);
    vs(jv,kv)=v1(jv,kv)-f*dt*0.25*(u1(jv,kv)+u1(jv+1,kv)+u1(jv,kv-1)+u1(jv+1,kv-1))......
        -g*dt/dx*( eta(jv,kv)-eta(jv,kv-1) );
    vs(jv,kv)=vs(jv,kv)./(1+beta2(jv,kv)*dt);


    %================================= 第四步  ================================
    % Coeffs in the 2nd step, eq.3.5.2'
    lambd=zeros(Nx,Ny);
    lambdr=zeros(Nx,Ny);
    lambdl=zeros(Nx,Ny);
    nu=zeros(Nx,Ny);
    nuu=zeros(Nx,Ny);
    nud=zeros(Nx,Ny);

    lambd(ju,ku)=0.25*(us(ju-1,ku)+2*us(ju,ku)+us(ju+1,ku))*dt/(2*dx);
    nuu(ju,ku)=(vs(ju-1,ku+1)+vs(ju,ku+1))*dt/(4*dx);
    nud(ju,ku)=(vs(ju-1,ku)+vs(ju,ku))*dt/(4*dx);

    nu(jv,kv)=0.25*(vs(jv,kv-1)+2*vs(jv,kv)+vs(jv,kv+1))*dt/(2*dx);
    lambdr(jv,kv)=(us(jv+1,kv)+us(jv+1,kv-1))*dt/(4*dx);
    lambdl(jv,kv)=(us(jv,kv)+us(jv,kv-1))*dt/(4*dx);

    for jj=Nx-1:-1:2
        for kk=Ny-1:-1:1
            if kk>1
                u(jj,kk)=(1+lambd(jj,kk)+nud(jj,kk))*us(jj,kk);
            else
                u(jj,kk)=(1+lambd(jj,kk)+nud(jj,kk))*us(jj,kk);
            end
            u(jj,kk)=u(jj,kk)/(1+lambd(jj,kk)+nuu(jj,kk));
        end
    end

    for kk=Ny-1:-1:2
        for jj=Nx-1:-1:1
            if jj>1
                v(jj,kk)=(1+lambdl(jj,kk)+nu(jj,kk))*vs(jj,kk);
            else
                v(jj,kk)=(1+lambdl(jj,kk)+nu(jj,kk))*vs(jj,kk);
            end
            v(jj,kk)=v(jj,kk)/(1+lambdr(jj,kk)+nu(jj,kk));
        end
    end

    drawnow
    if mod(niter,ns)==0
        clf
        mesh((1:Nx)*dx*1.e-3,(1:Nx)*dx*1.e-3,eta);
        set(gca,'zlim',[-0.02  0.05])
        zlabel('\eta, m')
%                 imagesc((1:Nx)*dx*1.e-3,(1:Nx)*dx*1.e-3,eta');

        title([num2str(round(T/60)),' min'])
        xlabel('x, km')
        ylabel('y, km')
        zlabel('\eta, m')


        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if niter == 0
            imwrite(imind,cm,'assignment2_4_nonlinear-HN-eq3-4-2.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment2_4_nonlinear-HN-eq3-4-2.gif','gif','WriteMode','append');
        end

        dt=min([dx/sqrt(g*max(h(:))),Cd/f^2]); %time step, s

    end

    niter=niter+1;
    T=T+dt;
    ts=[ts sum(eta(:))];
end







