%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% u=displacement unknown   %%%%%%%%%%%%%%%
%%%% v=displacement test      %%%%%%%%%%%%%%%
%%%% sigma=stress unknown     %%%%%%%%%%%%%%%
%%%% tau=stress test          %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms uxx uxy uxz 
syms uyx uyy uyz
syms uzx uzy uzz

syms vxx vxy vxz
syms vyx vyy vyz
syms vzx vzy vzz

syms sxx sxy sxz
syms syx syy syz
syms szx szy szz

syms txx txy txz
syms tyx tyy tyz
syms tzx tzy tzz

syms dx_sxx dy_syy dz_szz
syms dx_txx dy_tyy dz_tzz

syms phi1x phi1y phi1z
syms phi2x phi2y phi2z
syms phi3x phi3y phi3z 

syms psi1x psi1y psi1z
syms psi2x psi2y psi2z
syms psi3x psi3y psi3z



Casym=0;
Cconst=1
Ceq=1

Lambda = 1
Mu = 1
Dim=3

beta=(1.0/(2*Mu))
alpha = (-beta * Lambda /(Dim*Lambda +2*Mu))


% beta=0.5
% alpha=-1.000000000000000e-01

% syms beta alpha
% syms Ceq Cconst Casym

syms fx fy fz

sigma=[ sxx, sxy, sxz ;...
        syx, syy, syz ;...
        szx, szy, szz ;];
    
tau=[txx, txy, txz ;...
     tyx, tyy, tyz ;...
     tzx, tzy, tzz ; ];
 
divsigma= [fx+dx_sxx; fy+dy_syy; fz+dz_szz];

divtau= [dx_txx;dy_tyy; dz_tzz];

Asigma=beta *sigma + alpha * (sxx+syy+szz) * [1 0 0; 0 1 0; 0 0 1];
Atau=  beta * tau  + alpha * (txx+tyy+tzz) * [1 0 0; 0 1 0; 0 0 1];



 epsilonu= [uxx                  0.5*(uxy+uyx)    0.5*(uxz+uzx) ;  ...
            0.5*(uxy+uyx)        uyy              0.5*(uyz+uzy) ;  ...
            0.5*(uxz+uzx)        0.5*(uyz+uzy)    uzz           ;]; 
        
 epsilonv= [vxx                  0.5*(vxy+vyx)    0.5*(vxz+vzx) ;  ...
            0.5*(vxy+vyx)        vyy              0.5*(vyz+vzy) ;  ...
            0.5*(vxz+vzx)        0.5*(vyz+vzy)    vzz           ;]; 

        
Constitutive1=Asigma -epsilonu;
Constitutive2=Atau - epsilonv;
Asymmetry1=sigma-transpose(sigma);
Asymmetry2=tau-transpose(tau);



F=  Cconst * sum(sum(Constitutive1 .*Constitutive2))
F= F + Ceq * sum(divsigma' * divtau);
F= F + 0.5*Casym * sum(sum(Asymmetry1.*Asymmetry2 ));
F1 = diff(F,txx)* txx + diff(F,txy)* txy + diff(F,txz)* txz + diff(F,dx_txx)* dx_txx;
F2 = diff(F,tyx)* tyx + diff(F,tyy)* tyy + diff(F,tyz)* tyz + diff(F,dy_tyy)* dy_tyy;
F3 = diff(F,tzx)* tzx + diff(F,tzy)* tzy + diff(F,tzz)* tzz + diff(F,dz_tzz)* dz_tzz;
F4 = diff(F,vxx)*vxx + diff(F,vxy)*vxy + diff(F,vxz)*vxz;
F5 = diff(F,vyx)*vyx + diff(F,vyy)*vyy + diff(F,vyz)*vyz;
F6 = diff(F,vzx)*vzx + diff(F,vzy)*vzy + diff(F,vzz)*vzz;

%%%% first row %%%%
F11 = diff(F1,sxx)*sxx+diff(F1,sxy)*sxy+diff(F1,sxz)*sxz+diff(F1,dx_sxx)*dx_sxx;
F12 = diff(F1,syx)*syx+diff(F1,syy)*syy+diff(F1,syz)*syz+diff(F1,dy_syy)*dy_syy;
F13 = diff(F1,szx)*szx+diff(F1,szy)*szy+diff(F1,szz)*szz+diff(F1,dz_szz)*dz_szz;
F14 = diff(F1,uxx)*uxx + diff(F1,uxy)*uxy + diff(F1,uxz)*uxz;
F15 = diff(F1,uyx)*uyx + diff(F1,uyy)*uyy + diff(F1,uyz)*uyz;
F16 = diff(F1,uzx)*uzx + diff(F1,uzy)*uzy + diff(F1,uzz)*uzz;
%%%% second row %%%%
F21 = diff(F2,sxx)*sxx+diff(F2,sxy)*sxy+diff(F2,sxz)*sxz+diff(F2,dx_sxx)*dx_sxx;
F22 = diff(F2,syx)*syx+diff(F2,syy)*syy+diff(F2,syz)*syz+diff(F2,dy_syy)*dy_syy;
F23 = diff(F2,szx)*szx+diff(F2,szy)*szy+diff(F2,szz)*szz+diff(F2,dz_szz)*dz_szz;
F24 = diff(F2,uxx)*uxx + diff(F2,uxy)*uxy + diff(F2,uxz)*uxz;
F25 = diff(F2,uyx)*uyx + diff(F2,uyy)*uyy + diff(F2,uyz)*uyz;
F26 = diff(F2,uzx)*uzx + diff(F2,uzy)*uzy + diff(F2,uzz)*uzz;
%%%% third row %%%%
F31 = diff(F3,sxx)*sxx+diff(F3,sxy)*sxy+diff(F3,sxz)*sxz+diff(F3,dx_sxx)*dx_sxx;
F32 = diff(F3,syx)*syx+diff(F3,syy)*syy+diff(F3,syz)*syz+diff(F3,dy_syy)*dy_syy;
F33 = diff(F3,szx)*szx+diff(F3,szy)*szy+diff(F3,szz)*szz+diff(F3,dz_szz)*dz_szz;
F34 = diff(F3,uxx)*uxx + diff(F3,uxy)*uxy + diff(F3,uxz)*uxz;
F35 = diff(F3,uyx)*uyx + diff(F3,uyy)*uyy + diff(F3,uyz)*uyz;
F36 = diff(F3,uzx)*uzx + diff(F3,uzy)*uzy + diff(F3,uzz)*uzz;
%%%% fourth row %%%%
F41 = diff(F4,sxx)*sxx+diff(F4,sxy)*sxy+diff(F4,sxz)*sxz+diff(F4,dx_sxx)*dx_sxx;
F42 = diff(F4,syx)*syx+diff(F4,syy)*syy+diff(F4,syz)*syz+diff(F4,dy_syy)*dy_syy;
F43 = diff(F4,szx)*szx+diff(F4,szy)*szy+diff(F4,szz)*szz+diff(F4,dz_szz)*dz_szz;
F44 = diff(F4,uxx)*uxx + diff(F4,uxy)*uxy + diff(F4,uxz)*uxz;
F45 = diff(F4,uyx)*uyx + diff(F4,uyy)*uyy + diff(F4,uyz)*uyz;
F46 = diff(F4,uzx)*uzx + diff(F4,uzy)*uzy + diff(F4,uzz)*uzz;
%%%% fifth row %%%%
F51 = diff(F5,sxx)*sxx+diff(F5,sxy)*sxy+diff(F5,sxz)*sxz+diff(F5,dx_sxx)*dx_sxx;
F52 = diff(F5,syx)*syx+diff(F5,syy)*syy+diff(F5,syz)*syz+diff(F5,dy_syy)*dy_syy;
F53 = diff(F5,szx)*szx+diff(F5,szy)*szy+diff(F5,szz)*szz+diff(F5,dz_szz)*dz_szz;
F54 = diff(F5,uxx)*uxx + diff(F5,uxy)*uxy + diff(F5,uxz)*uxz;
F55 = diff(F5,uyx)*uyx + diff(F5,uyy)*uyy + diff(F5,uyz)*uyz;
F56 = diff(F5,uzx)*uzx + diff(F5,uzy)*uzy + diff(F5,uzz)*uzz;
%%%% sixth row %%%%
F61 = diff(F6,sxx)*sxx+diff(F6,sxy)*sxy+diff(F6,sxz)*sxz+diff(F6,dx_sxx)*dx_sxx;
F62 = diff(F6,syx)*syx+diff(F6,syy)*syy+diff(F6,syz)*syz+diff(F6,dy_syy)*dy_syy;
F63 = diff(F6,szx)*szx+diff(F6,szy)*szy+diff(F6,szz)*szz+diff(F6,dz_szz)*dz_szz;
F64 = diff(F6,uxx)*uxx + diff(F6,uxy)*uxy + diff(F6,uxz)*uxz;
F65 = diff(F6,uyx)*uyx + diff(F6,uyy)*uyy + diff(F6,uyz)*uyz;
F66 = diff(F6,uzx)*uzx + diff(F6,uzy)*uzy + diff(F6,uzz)*uzz;

b1= diff(F,fx)* fx ;
b2= diff(F,fy)* fy ;
b3= diff(F,fz)* fz ;

Mss=[F11 F12 F13;
     F21 F22 F23;
     F31 F32 F33;];
Msu=[F14 F15 F16;
     F24 F25 F26;
     F34 F35 F36;];
Muu=[F44 F45 F46;
     F54 F55 F56;
     F64 F32 F66;];
syms x y z
syms ux uy uz
syms mu lambda



u=[sin(pi*x).*sin(pi*y).*sin(pi*z);sin(pi*x).*sin(pi*y).*sin(pi*z);sin(pi*x).*sin(pi*y).*sin(pi*z);];
gradu=[diff(u,x),diff(u,y),diff(u,z)];
epsu=0.5*(gradu+transpose(gradu));
divepsu=diff(u(1),x)+diff(u(2),y)+diff(u(3),z);

sigma=2*mu * epsu + lambda * divepsu * eye(3);
divsigma1 = diff(sigma(1,1),x) + diff(sigma(1,2),y) + diff(sigma(1,3),z);
divsigma2 = diff(sigma(2,1),x) + diff(sigma(2,2),y) + diff(sigma(2,3),z);
divsigma3 = diff(sigma(3,1),x) + diff(sigma(3,2),y) + diff(sigma(3,3),z);
x=0.5
y=0.5
z=0.5
lambda=1
mu=1
f(1)=-divsigma1;
f(2)=-divsigma2;
f(3)=-divsigma3;