%2 dimension 4 variables (Fig1g)
%EGF (E), Notch signal (N), Delta (D), AS-C (A)
%Hill function
dn=0.25;dc=0.25;
nt=1;kt=1;ktnt=kt^nt;%nt: Hill's co-efficient@kt:activation co-efficient  for trans-activation
nc=5;kc=0.09;kcnc=kc^nc;%nc: Hill's co-efficient@kc:activation co-efficient for cis-inhibition

kd=1;ae=1;ad=1;
ea=100;de=1;ke=1;kn=1;Sigma=0;
dx = 2;dx2=dx * dx;dt=0.01;Tmax = 2000;
Xmax=25;Xsteps=1:Xmax;
Emax=0.1; Nmax=0.05;Dmax=0.1;
colorlimit_u=[0 1]; colorlimit_e=[0 Emax]; colorlimit_n=[0 Nmax];colorlimit_d=[0 Dmax];

E= zeros(Xmax,Xmax,Tmax);
N= zeros(Xmax,Xmax,Tmax);
D= zeros(Xmax,Xmax,Tmax);
A= zeros(Xmax,Xmax,Tmax);
Temp= zeros(Xmax,Xmax,Tmax);
  
A(:,1,1)=0.90*ones(Xmax,1,1);
A(:,2,1)=0.31*ones(Xmax,1,1);
A(:,3,1)=0.02*ones(Xmax,1,1);
E(:,1,1)=0.054*ones(Xmax,1,1);
E(:,2,1)=0.021*ones(Xmax,1,1);
E(:,3,1)=0.0016*ones(Xmax,1,1); 
D(:,1,1)=0.062*ones(Xmax,1,1);
D(:,2,1)=0.021*ones(Xmax,1,1);
D(:,3,1)=0.0013*ones(Xmax,1,1);

for T=1: Tmax-1
    if A(13,13,T)>0.5
        break
    end  
    Etemp=E(:,:,T);
    Eright(:,Xmax)=Etemp(:,Xmax);Eright(:,1:Xmax-1)=Etemp(:,2:Xmax);
    Eleft(:,1)=Etemp(:,1);Eleft(:,2:Xmax)=Etemp(:,1:Xmax-1);
    Etop(Xmax,:)=Etemp(Xmax,:);Etop(1:Xmax-1,:)=Etemp(2:Xmax,:);
    Ebottom(1,:)=Etemp(1,:);Ebottom(2:Xmax,:)=Etemp(1:Xmax-1,:);
    DeltaE=(Eright+Eleft+Etop+Ebottom-4*Etemp)/dx2;

    E(:,:,T+1)= E(:,:,T) + dt*( ae.*A(:,:,T).*(1-A(:,:,T)) - ke.*E(:,:,T) +de.*DeltaE);

    D(:,:,T+1)= D(:,:,T) + dt*( ad.*A(:,:,T).*(1-A(:,:,T)) - kd.*D(:,:,T) );     

    for X=2: Xmax-1
        for Y=2:Xmax-1
            SigmaD=D(Y,X+1,T)+D(Y,X-1,T)+D(Y+1,X,T)+D(Y-1,X,T);SigmaDNt=SigmaD.^nt;DNc=D(Y,X,T)^nc;
            N(Y,X,T+1)= N(Y,X,T) + dt*( dn* SigmaDNt./(ktnt+SigmaDNt) -dc*DNc/(kcnc+DNc) -kn*N(Y,X,T) );                     
        end
    end
    for X=2: Xmax-1 
        SigmaD=D(1,X+1,T)+D(1,X-1,T)+D(2,X,T);SigmaDNt=SigmaD.^nt;DNc=D(1,X,T)^nc;
        N(1,X,T+1)= N(1,X,T) + dt*( dn* SigmaDNt./(ktnt+SigmaDNt) -dc*DNc/(kcnc+DNc) -kn*N(1,X,T) );  
        SigmaD=D(Xmax,X+1,T)+D(Xmax,X-1,T)+D(Xmax-1,X,T);SigmaDNt=SigmaD.^nt;DNc=D(Xmax,X,T)^nc;
        N(Xmax,X,T+1)= N(Xmax,X,T) + dt*( dn* SigmaDNt./(ktnt+SigmaDNt) -dc*DNc/(kcnc+DNc) -kn*N(Xmax,X,T)  );         
    end
    for Y=2: Xmax-1        
        SigmaD=D(Y+1,1,T)+D(Y-1,1,T)+D(Y,2,T);SigmaDNt=SigmaD.^nt;DNc=D(Y,1,T)^nc;
        N(Y,1,T+1)= N(Y,1,T) + dt*( dn* SigmaDNt./(ktnt+SigmaDNt) -dc*DNc/(kcnc+DNc) -kn*N(Y,1,T)); 
        SigmaD=D(Y+1,Xmax,T)+D(Y-1,Xmax,T)+D(Y,Xmax-1,T);SigmaDNt=SigmaD.^nt;DNc=D(Y,Xmax,T)^nc;
        N(Y,Xmax,T+1)= N(Y,Xmax,T) + dt*( dn* SigmaDNt./(ktnt+SigmaDNt) -dc*DNc/(kcnc+DNc) -kn*N(Y,Xmax,T) );         
    end     
    SigmaD=D(1,2,T)+D(2,1,T);SigmaDNt=SigmaD.^nt;DNc=D(1,1,T)^nc;
    N(1,1,T+1)= N(1,1,T) + dt*( dn* SigmaDNt./(ktnt+SigmaDNt) -dc*DNc/(kcnc+DNc) -kn*N(1,1,T));  
    SigmaD=D(Xmax,2,T)+D(Xmax-1,1,T);SigmaDNt=SigmaD.^nt;DNc=D(Xmax,1,T)^nc;
    N(Xmax,1,T+1)= N(Xmax,1,T) + dt*( dn* SigmaDNt./(ktnt+SigmaDNt) -dc*DNc/(kcnc+DNc) -kn*N(Xmax,1,T) );  
    SigmaD=D(1,Xmax-1,T)+D(2,Xmax,T);SigmaDNt=SigmaD.^nt;DNc=D(1,Xmax,T)^nc;
    N(1,Xmax,T+1)= N(1,Xmax,T) + dt*( dn* SigmaDNt./(ktnt+SigmaDNt) -dc*DNc/(kcnc+DNc) -kn*N(1,Xmax,T) );
    SigmaD=D(Xmax,Xmax-1,T)+D(Xmax-1,Xmax,T);SigmaDNt=SigmaD.^nt;DNc=D(Xmax,Xmax,T)^nc;
    N(Xmax,Xmax,T+1)= N(Xmax,Xmax,T) + dt*( dn* SigmaDNt./(ktnt+SigmaDNt) -dc*DNc/(kcnc+DNc) -kn*N(Xmax,Xmax,T) ); 

    Temp(:,:,T)=(E(:,:,T) > N(:,:,T)).*(E(:,:,T)-N(:,:,T)-Sigma);
    A(:,:,T+1)=A(:,:,T)+dt*ea.*Temp(:,:,T).*(1-A(:,:,T));
end

Tmax=T;
% figure('Position',[2000 500 800 600]);colormap jet;
% subplot(2,2,1);imagesc(Xsteps, Xsteps, A(:,:,Tmax),colorlimit_u);set(gca,'YDir','normal');title('AS-C');colorbar;
% subplot(2,2,2);imagesc(Xsteps, Xsteps, E(:,:,Tmax),colorlimit_e);set(gca,'YDir','normal');title('EGF');colorbar;
% subplot(2,2,3);imagesc(Xsteps, Xsteps, N(:,:,Tmax),colorlimit_n);set(gca,'YDir','normal');title('Notch');colorbar;
% subplot(2,2,4);imagesc(Xsteps, Xsteps, D(:,:,Tmax),colorlimit_d);set(gca,'YDir','normal');title('Delta');colorbar;

figure('Position',[1800 800 800 400]);
Xmax=0.1;X=0:Xmax/20:Xmax;
SigmaXNt=(4*X).^nt;trans=dn*SigmaXNt./(ktnt+SigmaXNt);
XNc=X.^nc;cis=dc*XNc./(kcnc+XNc);
subplot(1,2,1); 
plot(X,cis,'r',X,trans,'k','LineWidth',2);xlim([0 Xmax]);ylim([0 0.16]);
for T=1:2:Tmax
subplot(1,2,2); 
plot(Xsteps,A(13,:,T),'b',Xsteps, 10*E(13,:,T),'k',Xsteps,10*D(13,:,T),'m',Xsteps,30*N(13,:,T),'g','LineWidth',2);ylim([-0.1 1]);
pause(0.01);
end

