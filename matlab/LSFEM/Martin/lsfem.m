function lsfem(u,v,nx,nt,dx,dt)

% diffusion constant
sigma = 1*ones(nx-1, nt-1);
% sigma = 0.1*ones(nx-1,nt-1);

% J = |u_t - div(v)|^2 + |v-sigma*u_x|^2 with div(v)=v_x
% and traditional Neumann b.c., i.e. v=0 on \partial\Omega
% Ansatz space: bilinear finite elements on cartesian grid u(t,x)
% vectorization of u: time-major (i.e. ut0x0, ut0x1, ...)
[u,v] = newtonStep(u,v,sigma,nx,nt,dx,dt);

%figure(1); surf(reshape(u,nx,nt));
%figure(2); surf(reshape(v,nx,nt));

[u,v,r] = newtonStep(u,v,sigma,nx,nt,dx,dt);
[u,v,r] = newtonStep(u,v,sigma,nx,nt,dx,dt);

%figure(3); surf(reshape(u,nx,nt));
%figure(4); surf(reshape(v,nx,nt));
%figure(5); surf(r/dx/dt);

%# keyboard

function [u,v,r] = newtonStep(u,v,sigma,nx,nt,dx,dt)
  J = 0;
  Ju = zeros(nx*nt,1);
  Jv = Ju;
  Juu = sparse(nx*nt,nx*nt);
  Juv = Juu;
  Jvu = Juu;
  Jvv = Juu;

  r = sigma;

  for it = 1:nt-1
    for ix = 1:nx-1
      idx = [(it-1)*nx+ix; (it-1)*nx+ix+1; it*nx+ix; it*nx+ix+1];
      [j,ju,jv,juu,juv,jvu,jvv] = Jloc(u(idx),v(idx),dt,dx,sigma(ix,it));
      J = J + j;
      r(ix,it) = j;
      Ju(idx) = Ju(idx)+ju;
      Jv(idx) = Jv(idx)+jv;
      Juu(idx,idx) = Juu(idx,idx)+juu;
      Juv(idx,idx) = Juv(idx,idx)+juv;
      Jvu(idx,idx) = Jvu(idx,idx)+jvu;
      Jvv(idx,idx) = Jvv(idx,idx)+jvv;
    end
  end

  % take out degrees of freedom with Dirichlet b.c.. This is u for t=0 and v on
  % the spatial domain boundaries
  idxu = nx+1:nx*nt;
  idxv = [];
  for it=1:nt
    idxv = [idxv (it-1)*nx+(2:nx-1)];
  end
  ddJ = [Juu(idxu,idxu) Juv(idxu,idxv);
         Jvu(idxv,idxu) Jvv(idxv,idxv)];
  dudv = - (ddJ \ [Ju(idxu); Jv(idxv)]);
  u(idxu) = u(idxu) + dudv(1:length(idxu));
  v(idxv) = v(idxv) + dudv(length(idxu)+1:end);
  
%   Juu = 0.5*full(Juu);
%   Juv = 0.5*full(Juv);
%   Jvu = 0.5*full(Jvu);
%   Jvv = 0.5*full(Jvv);

  
  save('matrices_Mart.mat', 'Juu', 'Juv', 'Jvu', 'Jvv');
  fprintf('residual J = %g\n',J);

% Jloc computes J on a single cell for given u and v
% in vector form as [u00, u01, u10 u11]' (entries utx)
function [J,Ju,Jv,Juu,Juv,Jvu,Jvv] = Jloc(u,v,dt,dx,sigma)
  % mass matrix of bilinear elements
  M = [ 4 2 2 1
        2 4 1 2
        2 1 4 2
        1 2 2 4 ] * dx*dt/36;
  % Differentiation matrix mapping space gradient to bilinear basis
  DxB = [ -1 1  0 0
          -1 1  0 0
           0 0 -1 1
           0 0 -1 1 ] / dx;
  % Differentiation matrix mapping time gradient to bilinear basis
  DtB = [ -1  0 1 0
           0 -1 0 1
          -1  0 1 0
           0 -1 0 1 ] / dt;


  J = (DtB*u-DxB*v)'*M*(DtB*u-DxB*v) + (v-sigma*DxB*u)' * M * (v-sigma*DxB*u);
  Ju = 2*DtB'*M*(DtB*u-DxB*v) - 2*sigma*DxB'*M*(v-sigma*DxB*u);
  Jv = -2*DxB'*M*(DtB*u-DxB*v) + 2*M*(v-sigma*DxB*u);

  Juu = 2*DtB'*M*DtB + 2*sigma*DxB'*M*sigma*DxB;
  Juv = -2*DtB'*M*DxB - 2*DxB'*sigma*M; 
  Jvu = -2*DxB'*M*DtB - 2*M*sigma*DxB;
  Jvv = 2*DxB'*M*DxB + 2*M;

