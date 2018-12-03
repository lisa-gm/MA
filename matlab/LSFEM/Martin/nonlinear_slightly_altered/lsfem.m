function lsfem(u,v,nx,nt,dx,dt,diffusion)

% diffusion constant
sigma = diffusion*ones(nx-1,nt-1);

% J = |u_t - div(v)-f(u)|^2 + |v-sigma*u_x|^2 with div(v)=v_x
% and traditional Neumann b.c., i.e. v=0 on \partial\Omega
% Ansatz space: bilinear finite elements on cartesian grid u(t,x)
% vectorization of u: time-major (i.e. ut0x0, ut0x1, ...)
damp = 0.0;
for k=1:41
  damp = (1+2*damp)/3;
  [u,v,~,damp] = newtonStep(u,v,sigma,nx,nt,dx,dt,damp);

  figure(1); surf(reshape(u,nx,nt));
  xlabel('time');
  ylabel('space');
  % figure(2); surf(reshape(v,nx,nt));

  fprintf('press return');
  pause
end

% # [u,v,r] = newtonStep(u,v,sigma,f,nx,nt,dx,dt);
% # [u,v,r] = newtonStep(u,v,sigma,f,nx,nt,dx,dt);
% # 
% # figure(3); surf(reshape(u,nx,nt));
% # figure(4); surf(reshape(v,nx,nt));
% # figure(5); surf(r/dx/dt);
% 
% # keyboard
end

function [u,v,r,damp] = newtonStep(u,v,sigma,nx,nt,dx,dt,damp)
  J = 0;
  Ju = zeros(nx*nt,1);
  Jv = Ju;
  Juu = sparse(nx*nt,nx*nt);
  Juv = Juu;
  Jvu = Juu;
  Jvv = Juu;
  
  J_u_s = zeros(nx*nt,1);
  J_u_u = zeros(nx*nt,1);
  
  J_v_s = zeros(nx*nt,1);
  J_v_u = zeros(nx*nt,1);
 
  J_vt_fu = zeros(nx*nt,1);
  J_df_ut = zeros(nx*nt,1);
  J_df_sx = zeros(nx*nt,1);
  J_df_f = zeros(nx*nt,1);
  
  J_f_tx = zeros(nx*nt,1);

  r = sigma;

  for it = 1:nt-1
    for ix = 1:nx-1
      idx = [(it-1)*nx+ix; (it-1)*nx+ix+1; it*nx+ix; it*nx+ix+1];
      [j,ju,jv,juu,juv,jvu,jvv, j_vt_fu, j_df_ut, j_df_sx, j_df_f, j_f_tx, j_u_s, j_u_u, j_v_u, j_v_s] = Jloc(u(idx),v(idx),dt,dx,sigma(ix,it));
      J = J + j;
      r(ix,it) = j;
      Ju(idx) = Ju(idx)+ju;
      Jv(idx) = Jv(idx)+jv;
      
      J_u_u(idx) = J_u_u(idx) + j_u_u;
      J_u_s(idx) = J_u_s(idx) + j_u_s;
      
      J_v_u(idx) = J_v_u(idx) + j_v_u;
      J_v_s(idx) = J_v_s(idx) + j_v_s;
      
      J_vt_fu(idx) = J_vt_fu(idx) + j_vt_fu;
      J_df_ut(idx) = J_df_ut(idx) + j_df_ut;
      J_df_sx(idx) = J_df_sx(idx) + j_df_sx;
      J_df_f(idx) = J_df_f(idx) + j_df_f;
      
      J_f_tx(idx) = J_f_tx(idx) + j_f_tx;
      
      Ju_comp = J_u_u +J_u_s + J_vt_fu + J_df_ut + J_df_sx + J_df_f;
      Jv_comp = J_v_u + J_v_s + J_f_tx;
      
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
  fprintf('norm(dJ) = %g   ',norm([Ju(idxu); Jv(idxv)]));
  ddJ = [Juu(idxu,idxu) Juv(idxu,idxv);
         Jvu(idxv,idxu) Jvv(idxv,idxv)];
  dudv = - (ddJ \ [Ju(idxu); Jv(idxv)]);
  
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
    
      fprintf('norm Ju: %d\n', norm(Ju));
      fprintf('norm Jv: %d\n', norm(Jv));

      fprintf('norm Ju_comp: %d\n', norm(Ju_comp));
      fprintf('norm Jv_comp: %d\n', norm(Jv_comp));
      
     
    
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
function [J,Ju,Jv,Juu,Juv,Jvu,Jvv, J_vt_fu, J_df_ut, J_df_sx, J_df_f, J_f_tx, J_u_s, J_u_u, J_v_u, J_v_s] = Jloc(u,v,dt,dx,sigma)
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
  [fu,dfu,ddfu] = F(u);

  J = (DtB*u-DxB*v-fu)'*M*(DtB*u-DxB*v-fu) ...
      + (v-sigma*DxB*u)' * M * (v-sigma*DxB*u);
  Ju = 2*(DtB-diag(dfu)*eye(4))'*M*(DtB*u-DxB*v-fu) ...
       - 2*sigma*DxB'*M*(v-sigma*DxB*u);
  Jv = -2*DxB'*M*(DtB*u-DxB*v-fu) + 2*M*(v-sigma*DxB*u);
  
  
  J_vt_fu = -2*DtB'*M*fu;
  J_df_ut = -2*diag(dfu)*M*DtB*u;
  J_df_sx = 2*diag(dfu)*M*DxB*v;
  J_df_f = 2*diag(dfu)*M*fu;
   
  J_u_u = 2*sigma*DxB'*M*sigma*DxB*u + 2*DtB'*M*DtB*u;
  J_u_s = -2*DtB'*M*DxB*v - 2*sigma*DxB'*M*v;

  Jv = -2*DxB'*M*(DtB*u-DxB*v-fu) + 2*M*(v-sigma*DxB*u);
  
  J_v_u = -2*DxB'*M*DtB*u - 2*M*sigma*DxB*u;
  J_v_s = 2*DxB'*M*DxB*v + 2*M*v;
  J_f_tx = 2*DxB'*M*fu;
   

  Juu = -2*diag(ddfu)'*M*(DtB*u-DxB*v-fu) + 2*(DtB-diag(dfu))'*M*(DtB-diag(dfu)) ...
        + 2*sigma*DxB'*M*sigma*DxB;
  Juv = -2*(DtB-diag(dfu))'*M*DxB - 2*DxB'*sigma*M;
  Jvu = -2*DxB'*M*(DtB-diag(dfu)) - 2*M*sigma*DxB;
  Jvv = 2*DxB'*M*DxB + 2*M;
end

function [f,df,ddf] = F(u)
% #   f = -4*u;
% #   df = -4*ones(length(u),1);
% #   ddf = 0*u;
  alpha = 10;
  th = 0.1;
  f = -alpha*u.*(u-th).*(u-1);
  df = -alpha*((u-th).*(u-1) + alpha*u.*(u-1) + alpha*u.*(u-th));
  ddf = -alpha*((u-th)+alpha*(u-1) + alpha*u+alpha*(u-1) + alpha*u+alpha*(u-th));
end

