function lsfem(u,v,nx,nt,dx,dt,diffusion)
eps = 10^(-7);
% diffusion constant
sigma = diffusion*ones(nx-1,nt-1);

% J = |u_t - div(v)-f(u)|^2 + |v-sigma*u_x|^2 with div(v)=v_x
% and traditional Neumann b.c., i.e. v=0 on \partial\Omega
% Ansatz space: bilinear finite elements on cartesian grid u(t,x)
% vectorization of u: time-major (i.e. ut0x0, ut0x1, ...)
damp = 0.0;
for k=1:120 %41
    
  damp = (1+2*damp)/3;
  [u,v,~,damp, norm_res_J] = newtonStep(u,v,sigma,nx,nt,dx,dt,damp);

  %figure(1); surf(reshape(u,nx,nt));
  %figure(2); surf(reshape(v,nx,nt));
   
  if(norm_res_J < eps)
      fprintf('convergence criterion reached\n');
      break;
  end
  
  %fprintf('press return');
  %pause
end

figure(1); surf(reshape(u,nx,nt));

%# [u,v,r] = newtonStep(u,v,sigma,f,nx,nt,dx,dt);
%# [u,v,r] = newtonStep(u,v,sigma,f,nx,nt,dx,dt);
% 
%# figure(3); surf(reshape(u,nx,nt));
%# figure(4); surf(reshape(v,nx,nt));
%# figure(5); surf(r/dx/dt);

%# keyboard
end

function [u,v,r,damp, norm_res_J] = newtonStep(u,v,sigma,nx,nt,dx,dt,damp)
  J = 0;
  Ju = zeros(nx*nt,1);
  Jv = Ju;
  Juu = sparse(nx*nt,nx*nt);
  Juv = Juu;
  Jvu = Juu;
  Jvv = Juu;
  
  Jvu_lin = zeros(nx*nt);
  Jvu_nonlin = zeros(nx*nt);
  
  Juu_lin  = zeros(nx*nt);
  Juu_d2f_ut = zeros(nx*nt);
  Juu_d2f_sx = zeros(nx*nt); 
  Juu_d2f_f  = zeros(nx*nt);
  Juu_df_df = zeros(nx*nt);
  Juu_df_vt = zeros(nx*nt);
  Juu_vt_df = zeros(nx*nt);

  r = sigma;
  
  for it = 1:nt-1
    for ix = 1:nx-1
      idx = [(it-1)*nx+ix; (it-1)*nx+ix+1; it*nx+ix; it*nx+ix+1];
      [j,ju,jv,juu,juv,jvu,jvv, jvu_lin, jvu_nonlin, juu_lin, juu_d2f_ut, juu_d2f_sx, juu_d2f_f, juu_df_df, juu_df_vt, juu_vt_df] = Jloc(u(idx),v(idx),dt,dx,sigma(ix,it));
      J = J + j;
      r(ix,it) = j;
      Ju(idx) = Ju(idx)+ju;
      Jv(idx) = Jv(idx)+jv;
      
      
      Juu(idx,idx) = Juu(idx,idx)+juu;
      Juv(idx,idx) = Juv(idx,idx)+juv;
      Jvu(idx,idx) = Jvu(idx,idx)+jvu;
      Jvv(idx,idx) = Jvv(idx,idx)+jvv;
      
      Juu_lin(idx,idx) = Juu_lin(idx,idx) + juu_lin;
      Juu_d2f_ut(idx,idx) =  Juu_d2f_ut(idx,idx) + juu_d2f_ut;
      Juu_d2f_sx(idx,idx) = Juu_d2f_sx(idx,idx) + juu_d2f_sx;
      Juu_d2f_f(idx,idx) = Juu_d2f_f(idx,idx) + juu_d2f_f;
      Juu_df_df(idx,idx) = Juu_df_df(idx,idx) + juu_df_df;
      Juu_df_vt(idx,idx) = Juu_df_vt(idx,idx) + juu_df_vt;
      Juu_vt_df(idx,idx) = Juu_vt_df(idx,idx) + juu_vt_df;
      
      Jvu_lin(idx, idx) = Jvu_lin(idx, idx) + jvu_lin;
      Jvu_nonlin(idx, idx) = Jvu_nonlin(idx, idx) + jvu_nonlin;
      
    end
  end
  
  % take out degrees of freedom with Dirichlet b.c.. This is u for t=0 and v on
  % the spatial domain boundaries
  idxu = nx+1:nx*nt;
  idxv = [];
  for it=1:nt
    idxv = [idxv (it-1)*nx+(2:nx-1)];
  end
  fprintf('norm(dJ) = %g   ',norm([Ju(idxu); Jv(idxv)]));
  fprintf('norm(dJ_u) = %g   ',norm(Ju(idxu)));
  fprintf('norm(dJ_v) = %g   ',norm(Jv(idxv)));

  norm_res_J = norm([Ju(idxu); Jv(idxv)]);
  
  ddJ = [Juu(idxu,idxu) Juv(idxu,idxv);
         Jvu(idxv,idxu) Jvv(idxv,idxv)];
  dudv = - (ddJ \ [Ju(idxu); Jv(idxv)]);
  
    
    fprintf('norm Jv = %g\n',norm(Jv(idxv)));  
    fprintf('norm Ju = %g\n',norm(Ju(idxu)));
    
    fprintf('norm Jvu = %g\n',norm(full(Jvu(idxv, idxu))));
    fprintf('norm Jvu_lin = %g\n',norm(full(Jvu_lin(idxv, idxu))));
    fprintf('norm Jvu_nonlin = %g\n',norm(full(Jvu_nonlin(idxv, idxu))));
    
    fprintf('norm Juu_lin = %g\n',norm(full(Juu_lin(idxu, idxu))));
    fprintf('norm Juu_d2f_ut = %g\n',norm(full(Juu_d2f_ut(idxu, idxu))));
    fprintf('norm Juu_d2f_sx = %g\n',norm(full(Juu_d2f_sx(idxu, idxu))));
    fprintf('norm Juu_d2f_f = %g\n',norm(full(Juu_d2f_f(idxu, idxu))));
    fprintf('norm Juu_df_df = %g\n',norm(full(Juu_df_df(idxu, idxu))));
    fprintf('norm Juu_df_vt = %g\n',norm(full(Juu_df_vt(idxu, idxu))));
    fprintf('norm Juu_vt_df = %g\n',norm(full(Juu_vt_df(idxu, idxu))));
  
%   norm(dudv)
%    norm_ddJ = norm(full(ddJ))
%    norm_Juu = norm(full(Juu(idxu,idxu)))
%    norm_Jvu = norm(full(Jvu(idxv,idxu)))
%    norm_Juv = norm(full(Juv(idxu,idxv)))
%    norm_Jvv = norm(full(Jvv(idxv,idxv)))

%   norm([Ju(idxu); Jv(idxv)])
%   

  
  % quick hack for nonconvexity:
  if dudv'*[Ju(idxu); Jv(idxv)] > 0
    dudv = -dudv;
  end
  
  utrial = u;
  vtrial = v;
  while true
    utrial(idxu) = u(idxu) + damp*dudv(1:length(idxu));
    vtrial(idxv) = v(idxv) + damp*dudv(length(idxu)+1:end);
   
    Jtrial = 0;
    for it = 1:nt-1
      for ix = 1:nx-1
        idx = [(it-1)*nx+ix; (it-1)*nx+ix+1; it*nx+ix; it*nx+ix+1];
        j = Jloc(utrial(idx),vtrial(idx),dt,dx,sigma(ix,it));
        Jtrial = Jtrial + j;
      end
    end
   
    fprintf('damping = %g, Jnew = %g, Jold = %g\n',damp,Jtrial,J);


    if Jtrial < J
      break;
    end
    
    damp = damp/2;
  end
 
  
  u = utrial;
  v = vtrial;
  
  fprintf('residual J = %g\n',J);
end

% Jloc computes J on a single cell for given u and v
% in vector form as [u00, u01, u10 u11]' (entries utx)
function [J,Ju,Jv,Juu,Juv,Jvu,Jvv, Jvu_lin, Jvu_nonlin, Juu_lin, Juu_d2f_ut, Juu_d2f_sx, Juu_d2f_f, Juu_df_df, Juu_df_vt, Juu_vt_df] = Jloc(u,v,dt,dx,sigma)
  % mass matrix of bilinear elements
  M = [ 4 2 2 1
        2 4 1 2
        2 1 4 2
        1 2 2 4 ] * dx*dt / 36;
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
  % right hand side evaluated at node points
  
  %u = [1;0;1;0];
  %v = [6;2;4;-5];
  
  [fu,dfu,ddfu] = F(u);

  J = (DtB*u-DxB*v-fu)'*M*(DtB*u-DxB*v-fu) ...
      + (v-sigma*DxB*u)' * M * (v-sigma*DxB*u);
  
  J_ut = (DtB*u)'*M*(DtB*u);
  J_ux = (-sigma*DxB*u)' * M * (-sigma*DxB*u);
 
  
  Ju = 2*(DtB-diag(dfu))'*M*(DtB*u-DxB*v-fu) ...
       - 2*sigma*DxB'*M*(v-sigma*DxB*u);
  Jv = -2*DxB'*M*(DtB*u-DxB*v-fu) + 2*M*(v-sigma*DxB*u);

  Juu = -2*diag(ddfu.*M*(DtB*u-DxB*v-fu)) + 2*(DtB-diag(dfu))'*M*(DtB-diag(dfu)) ...
        + 2*sigma*DxB'*M*sigma*DxB;
    
  Juu_lin = 2*DtB'*M*DtB + 2*sigma*DxB'*M*sigma*DxB;
  Juu_d2f_ut = -2*diag(ddfu.*M*(DtB*u));
  Juu_d2f_sx = -2*diag(ddfu.*M*(-DxB*v));
  Juu_d2f_f = -2*diag(ddfu.*M*(-fu));
  Juu_df_vt = 2*(-diag(dfu))'*M*DtB;
  Juu_df_df = 2*(-diag(dfu))'*M*(-diag(dfu));
  Juu_vt_df = 2*(DtB)'*M*(-diag(dfu)); 
  
  Juv = -2*(DtB-diag(dfu))'*M*DxB - 2*DxB'*sigma*M;
  Jvu = -2*DxB'*M*(DtB-diag(dfu)) - 2*M*sigma*DxB;
  
  Jvu_lin = -2*DxB'*M*(DtB) - 2*M*sigma*DxB;
  Jvu_nonlin = -2*DxB'*M*diag(dfu);
  
  Jvv = 2*DxB'*M*DxB + 2*M;
end

function [f,df,ddf] = F(u)
%#   f = -4*u;
%#   df = -4*ones(length(u),1);
%#   ddf = 0*u;
  alpha = 1; 
  th = 0.1;
  f = alpha*(-u.*(u-th).*(u-1));
  df = alpha*(-((u-th).*(u-1) + u.*(u-1) + u.*(u-th)));
  ddf = alpha*(-((u-th)+(u-1) + u+(u-1) + u+(u-th)));
end

