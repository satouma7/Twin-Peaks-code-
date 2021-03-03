%2 dimension 5 variables Phase Diagram (Supplementary Fig12)
%EGF (E), full length Notch (F), Notch signal (S), Delta (D), AS-C (A)
%Hill function
%full length Notch degrades by cis-inhibition and trans-activation
dn=0.25;%trans activation:dt
dc=5;%cis inhibition
n=1;%full length Notch production
kf=0.1;%full length Notch degradation (spontaneous)
kt=0.1;%full length Notch degradation by trans-activation
nt=1;K1=0.3;Knt=K1^nt;%nt: Hill's co-efficient@K1:activation co-efficient  for trans-activation

kd=1;ae=1;ad=1;
ea=10;de=5;ke=1;ks=2;Sigma=0;
dx = 2;dx2=dx * dx;dt=0.01;Tmax = 800;
Xmax=25;Xsteps=1:Xmax;
Emax=0.2; Nmax=0.1;Dmax=0.2;
colorlimit_u=[0 1]; colorlimit_e=[0 Emax]; colorlimit_n=[0 Nmax];colorlimit_d=[0 Dmax];

nci=[1 3 5 7];%nci: Hill's co-efficient@for cis-inhibition
K2j=[0.05 0.1 0.2 0.3];%K2j:activation co-efficient for cis-inhibition
Imax=size(nci,2);Jmax=size(K2j,2);Count=1;
figure('Position',[2000 500 1000 1000]);colormap jet;

for I=1:Imax%dc
    nc=nci(I);
    E= zeros(Xmax,Xmax,Imax,Jmax);
    F= ones(Xmax,Xmax,Imax,Jmax);%full length Notch
    S= zeros(Xmax,Xmax,Imax,Jmax);%SuH activity
    D= zeros(Xmax,Xmax,Imax,Jmax);
    A= zeros(Xmax,Xmax,Imax,Jmax);
    A(:,1,:,:)=0.90*ones(Xmax,1,Imax,Jmax);
    A(:,2,:,:)=0.31*ones(Xmax,1,Imax,Jmax);
    A(:,3,:,:)=0.02*ones(Xmax,1,Imax,Jmax);
    E(:,1,:,:)=0.054*ones(Xmax,1,Imax,Jmax);
    E(:,2,:,:)=0.021*ones(Xmax,1,Imax,Jmax);
    E(:,3,:,:)=0.0016*ones(Xmax,1,Imax,Jmax); 
    D(:,1,:,:)=0.062*ones(Xmax,1,Imax,Jmax);
    D(:,2,:,:)=0.021*ones(Xmax,1,Imax,Jmax);
    D(:,3,:,:)=0.0013*ones(Xmax,1,Imax,Jmax);
    
    for J=1:Jmax%ec Parfor
        K2=K2j(J);Knc=K2^nc;
        for T=1: Tmax-1
            Etemp=E(:,:,I,J);Ftemp=F(:,:,I,J);Stemp=S(:,:,I,J);Dtemp=D(:,:,I,J);Atemp=A(:,:,I,J);
            ESH= zeros(Xmax,Xmax);SigmaD= zeros(Xmax,Xmax);DeltaE= zeros(Xmax,Xmax);Hill= zeros(Xmax,Xmax);
            Eright=zeros(Xmax,Xmax);Eleft=zeros(Xmax,Xmax);Etop=zeros(Xmax,Xmax);Ebottom=zeros(Xmax,Xmax);
            Eright(:,Xmax)=Etemp(:,Xmax);Eright(:,1:Xmax-1)=Etemp(:,2:Xmax);
            Eleft(:,1)=Etemp(:,1);Eleft(:,2:Xmax)=Etemp(:,1:Xmax-1);
            Etop(Xmax,:)=Etemp(Xmax,:);Etop(1:Xmax-1,:)=Etemp(2:Xmax,:);
            Ebottom(1,:)=Etemp(1,:);Ebottom(2:Xmax,:)=Etemp(1:Xmax-1,:);
            DeltaE=(Eright+Eleft+Etop+Ebottom-4*Etemp)/dx2;
            E(:,:,I,J)= Etemp + dt*( ae*Atemp.*(1-Atemp) - ke*Etemp +de*DeltaE);
            D(:,:,I,J)= Dtemp + dt*( ad*Atemp.*(1-Atemp) - kd*Dtemp );     
            
            Dright=zeros(Xmax,Xmax);Dright(:,1:Xmax-1)=Dtemp(:,2:Xmax);
            Dleft=zeros(Xmax,Xmax);Dleft(:,2:Xmax)=Dtemp(:,1:Xmax-1);
            Dtop=zeros(Xmax,Xmax);Dtop(1:Xmax-1,:)=Dtemp(2:Xmax,:);
            Dbottom=zeros(Xmax,Xmax);Dbottom(2:Xmax,:)=Dtemp(1:Xmax-1,:);
            SigmaD=Dright+Dleft+Dtop+Dbottom;
            SigmaDNt=SigmaD.^nt;Hill=SigmaDNt./(Knt+SigmaDNt);
            S(:,:,I,J)= Stemp + dt*( dn* Hill.*Ftemp  -ks*Stemp); 
            
            DNc=Dtemp.^nc;
            F(:,:,I,J)= Ftemp + dt*( n*(1-Ftemp) - dc*DNc./(Knc+DNc).*Ftemp-kf*Ftemp -kt*Hill.*Ftemp );
            ESH=E(:,:,I,J)-S(:,:,I,J);ESH=ESH.*(ESH>0);
            A(:,:,I,J)=Atemp +dt*ea*ESH.*(1-Atemp);
        end    
    end
    for J=1:Jmax
        subplot(Imax,Jmax,Count);plot(1:Xmax,S(13,:,I,J)*10,'r',1:Xmax,F(13,:,I,J),'k','LineWidth',2);%ylim([0 0.05]);
        title(['S/F: nc=', num2str(nci(I)),' K2=', num2str(K2j(J))]);
        Count=Count+1;
    end
end

figure('Position',[2000 400 1000 1000]);
Xmax=0.1;X=0:Xmax/20:Xmax;Count=1;
for I=1:Imax%dc
    nc=nci(I);
    for J=1:Jmax%ec
        K2=K2j(J);Knc=K2^nc;
        SigmaDNt=(4*X).^nt;trans=dn* SigmaDNt./(Knt+SigmaDNt);
        DNc=X.^nc;cis=dc* DNc./(Knc+DNc);
        subplot(Imax,Jmax,Count);
        plot(X,cis,'r',X,trans,'k','LineWidth',2);xlim([0 Xmax]);ylim([0 0.2]);
        title(['nc=', num2str(nc),' K2=', num2str(K2)]);
        Count=Count+1;
    end
end