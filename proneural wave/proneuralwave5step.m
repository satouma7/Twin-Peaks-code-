%2 dimension 5 variables 
%EGF (E), full length Notch (F), Notch signal (S), Delta (D), AS-C (A)
%Step function
%full length Notch degrades by cis-inhibition and trans-activation
dn=0.25;%trans activation
dc=10;%cis inhibition 5<dc<10
ec=0.1;
n=1;%full length Notch production
kf=0.1;%full length Notch degradation (spontaneous)
kt=0.1;%full length Notch degradation by trans-activation

kd=1;ae=1;ad=1;
ea=10;de=5;ke=1;ks=1;
dx = 2;dx2=dx * dx;dt=0.01;Tmax = 1000;
Xmax=25;Xsteps=1:Xmax;
Emax=0.2; Nmax=0.1;Dmax=0.2;
colorlimit_u=[0 1]; colorlimit_e=[0 Emax]; colorlimit_n=[0 Nmax];colorlimit_d=[0 Dmax];

E= zeros(Xmax,Xmax,Tmax);
F= ones(Xmax,Xmax,Tmax);%full length Notch
S= zeros(Xmax,Xmax,Tmax);%SuH activity
%H= zeros(Xmax,Xmax,Tmax);%H binding
D= zeros(Xmax,Xmax,Tmax);
A= zeros(Xmax,Xmax,Tmax);
ESH= zeros(Xmax,Xmax);
SigmaD= zeros(Xmax,Xmax);
  
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
    Etemp=E(:,:,T);
    Eright(:,Xmax)=Etemp(:,Xmax);Eright(:,1:Xmax-1)=Etemp(:,2:Xmax);
    Eleft(:,1)=Etemp(:,1);Eleft(:,2:Xmax)=Etemp(:,1:Xmax-1);
    Etop(Xmax,:)=Etemp(Xmax,:);Etop(1:Xmax-1,:)=Etemp(2:Xmax,:);
    Ebottom(1,:)=Etemp(1,:);Ebottom(2:Xmax,:)=Etemp(1:Xmax-1,:);
    DeltaE=(Eright+Eleft+Etop+Ebottom-4*Etemp)/dx2;
    E(:,:,T+1)= E(:,:,T) + dt*( ae*A(:,:,T).*(1-A(:,:,T)) - ke*E(:,:,T) +de*DeltaE);
    D(:,:,T+1)= D(:,:,T) + dt*( ad*A(:,:,T).*(1-A(:,:,T)) - kd*D(:,:,T) );     
        
    Dtemp=D(:,:,T);
    Dright=zeros(Xmax,Xmax);Dright(:,1:Xmax-1)=Dtemp(:,2:Xmax);
    Dleft=zeros(Xmax,Xmax);Dleft(:,2:Xmax)=Dtemp(:,1:Xmax-1);
    Dtop=zeros(Xmax,Xmax);Dtop(1:Xmax-1,:)=Dtemp(2:Xmax,:);
    Dbottom=zeros(Xmax,Xmax);Dbottom(2:Xmax,:)=Dtemp(1:Xmax-1,:);
    SigmaD=Dright+Dleft+Dtop+Dbottom;
    S(:,:,T+1)= S(:,:,T) + dt*( dn* SigmaD.*F(:,:,T)  -ks*S(:,:,T) ); 
    F(:,:,T+1)= F(:,:,T) + dt*( n*(1-F(:,:,T)) - dc*(D(:,:,T)>=ec).*F(:,:,T)-kf*F(:,:,T) -kt*SigmaD.*F(:,:,T) );

    ESH=E(:,:,T+1)-S(:,:,T+1);ESH=ESH.*(ESH>0);
    A(:,:,T+1)=A(:,:,T)+dt*ea*ESH.*(1-A(:,:,T));
end

figure('Position',[2000 500 800 900]);colormap jet;
for T=1:20:800
    subplot(3,2,1);imagesc(A(:,:,T),colorlimit_u);set(gca,'YDir','normal');title('AS-C');colorbar;
    subplot(3,2,2);imagesc(E(:,:,T),colorlimit_e);set(gca,'YDir','normal');title('EGF');colorbar;
    subplot(3,2,3);imagesc(S(:,:,T),colorlimit_n);set(gca,'YDir','normal');title('Su(H)');colorbar;
    subplot(3,2,4);imagesc(D(:,:,T),colorlimit_d);set(gca,'YDir','normal');title('Delta');colorbar;
    subplot(3,2,5);imagesc(F(:,:,T),colorlimit_u);set(gca,'YDir','normal');title('full length N');colorbar;
    subplot(3,2,6);
    plot(Xsteps,A(13,:,T),'b',Xsteps, 5*E(13,:,T),'k',Xsteps,5*D(13,:,T),'g',Xsteps,20*S(13,:,T),'r',Xsteps,F(13,:,T),'m');ylim([-0.1 1.2]);
    pause(0.001);
end