%2 dimension 5 variables Phase Diagram
%EGF (E), full length Notch (F), Notch signal (S), Delta (D), AS-C (A)
%Step function
%full length Notch degrades by cis-inhibition and trans-activation
dn=0.25;%trans activation
n=1;%full length Notch production
kf=0.1;%full length Notch degradation (spontaneous)
kt=0.1;%full length Notch degradation by trans-activation

kd=1;ae=1;ad=1;
ea=10;de=5;ke=1;ks=1;
dx = 2;dx2=dx * dx;dt=0.01;Tmax = 800;
Xmax=25;Xsteps=1:Xmax;
Emax=0.2; Nmax=0.1;Dmax=0.2;
colorlimit_u=[0 1]; colorlimit_e=[0 Emax]; colorlimit_n=[0 Nmax];colorlimit_d=[0 Dmax];

dci=[1 3 4 5];%cis inhibition
ecj=[0.05 0.08 0.10 0.12];%threshold
Imax=size(dci,2);Jmax=size(ecj,2);Count=1;
figure('Position',[2000 500 1000 1000]);colormap jet;

for I=1:Imax%dc
    dc=dci(I);
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
        ec=ecj(J);
        for T=1: Tmax-1
            Etemp=E(:,:,I,J);Ftemp=F(:,:,I,J);Stemp=S(:,:,I,J);Dtemp=D(:,:,I,J);Atemp=A(:,:,I,J);
            ESH= zeros(Xmax,Xmax);SigmaD= zeros(Xmax,Xmax);DeltaE= zeros(Xmax,Xmax);
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
            S(:,:,I,J)= Stemp + dt*( dn* SigmaD.*Ftemp  -ks*Stemp); 
            F(:,:,I,J)= Ftemp + dt*( n*(1-Ftemp) - dc*(Dtemp>=ec).*Ftemp-kf*Ftemp -kt*SigmaD.*Ftemp );
            ESH=E(:,:,I,J)-S(:,:,I,J);ESH=ESH.*(ESH>0);
            A(:,:,I,J)=Atemp +dt*ea*ESH.*(1-Atemp );
        end    
    end
    for J=1:Jmax
        subplot(Imax,Jmax,Count);plot(1:Xmax,S(13,:,I,J),1:Xmax,D(13,:,I,J));%ylim=[0 Nmax];
        title(['S: dc=', num2str(dci(I)),' ec=', num2str(ecj(J))]);
        Count=Count+1;
    end
end

figure('Position',[2000 400 1000 1000]);
Xmax=0.12;X=0:Xmax/20:Xmax;Count=1;
for I=1:Imax%dc
    dc=dci(I);
    for J=1:Jmax%ec
        ec=ecj(J);       
        subplot(Imax,Jmax,Count);
        plot(X,dc*(X>=ec),'r',X,4*dn*X,'k');xlim([0 Xmax]);ylim([0 max(dci)]);
        title(['S: dc=', num2str(dc),' ec=', num2str(ec)]);
        Count=Count+1;
    end
end