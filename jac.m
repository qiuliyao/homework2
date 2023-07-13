% function J1=jac(psi,zeta)

% global coeffi_Jacobian
J=zeros(size(psi));
% [ny,nx]=size(psi);
% Arakawa 1966 Eq(46), but note that i, j in the equation are to the East and North respectively,
%  but here jy is to the south therefore j+1£¬in the original eq 46 is jy-1 in our martrix, j - 1 is jy+1 to the
% South,   % J1=Jacobian(zeta,psi)

% interior points
ixx=2:N-1;jyy=2:N-1;
ixxm=ixx-1; ixxp=ixx+1;
jyym=jyy-1; jyyp=jyy+1;
Arakawa66_eq46;


%     if double periodic b.c is implied, then 4 corner points and 4
%     boundaries  need spectial attention 



ixx=1;   jyy=1; % 1 of 4 corner points
ixxm=N; ixxp=ixx+1;
jyym=N; jyyp=jyy+1;
Arakawa66_eq46;


ixx=1;jyy=N;%lower left corner
ixxm=N; ixxp=ixx+1;
jyym=jyy-1; jyyp=1;
Arakawa66_eq46;

ixx=N;  jyy=1; % upper right
ixxm=ixx-1; ixxp=1;
jyym=N; jyyp=jyy+1;
Arakawa66_eq46;

ixx=N;  jyy=N; % lower right
ixxm=ixx-1; ixxp=1;
jyym=jyy-1; jyyp=1;
Arakawa66_eq46;



% -------------------4 boundaries with corner points excluded-----------
ixx=2:N-1; jyy=1;  % south line
ixxm=ixx-1; ixxp=ixx+1;
jyym=N; jyyp=jyy+1;
Arakawa66_eq46;


ixx=2:N-1; jyy=N;  % north line;
ixxm=ixx-1; ixxp=ixx+1;
jyym=jyy-1; jyyp=1;
Arakawa66_eq46;


ixx=1; jyy=2:N-1;  % west line
ixxm=N; ixxp=ixx+1;
jyym=jyy-1; jyyp=jyy+1;
Arakawa66_eq46;


ixx=N; jyy=2:N-1;  % east line
ixxm=ixx-1; ixxp=1;
jyym=jyy-1; jyyp=jyy+1;
Arakawa66_eq46;

% end if 

J0=coeffi_Jacobian*J0;
% end


