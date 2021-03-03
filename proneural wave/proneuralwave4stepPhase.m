%2 dimension 4 variables Phase Diagram (Supplementary Fig1)
%EGF (E), Notch signal (N), Delta (D), AS-C (A)
%Step function
dn=0.25;dc=0.1;ec=0.07;
kd=1;ae=1;ad=1;
ea=100;de=1;ke=1;kn=1;Sigma=0;

dx = 2;dx2=dx * dx;dt=0.01;Tmax = 1000;
Xmax=25;Xsteps=1:Xmax;
Emax=0.1; Nmax=0.05;Dmax=0.1;Kmax=0.1;
colorlimit_u=[0 1]; colorlimit_e=[0 Emax]; colorlimit_n=[-0.01 Nmax];colorlimit_d=[0 Dmax];colorlimit_k=[0 Kmax];

E= zeros(Xmax,Xmax,Tmax);
N= zeros(Xmax,Xmax,Tmax);
D= zeros(Xmax,Xmax,Tmax);
A= zeros(Xmax,Xmax,Tmax);
K= zeros(Xmax,Xmax,Tmax);
Nseries=zeros(1,Xmax,10,10);
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

figure('Position',[300 300 1500 1500]);colormap jet;
I=1;alphaM=5;betaM=5;
for alpha=1:alphaM%dn & dc
    dc=0.05+alpha/50;%dt=0-0.2
    %dc=0;
    for beta=1:betaM%et and ec
        ec=0.03+beta/100;%et=0.03-0.07
        for T=1: Tmax-1
            Afront=A(13,13,T);       
            if Afront>0.5
                subplot(alphaM,betaM,I);
                plot(Xsteps, A(13,:,T),'k',Xsteps, N(13,:,T)*30,'r',Xsteps,zeros(1,Xmax),'k--','LineWidth',2);
                %plot(Xsteps, A(13,:,T),'b',Xsteps, D(13,:,T)*10,'g',Xsteps, N(13,:,T)*30,'r',Xsteps, E(13,:,T)*10,'k',Xsteps,zeros(1,Xmax),'k--','LineWidth',2);
                set(gca,'YDir','normal');xlim([1 Xmax]);ylim([-0.1 1]);
                title(strcat('dc=',num2str(dc),' ec=',num2str(ec)));
                I=I+1;
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
                    N(Y,X,T+1)= N(Y,X,T) + dt*( dn*( D(Y,X+1,T)+D(Y,X-1,T)+D(Y+1,X,T)+D(Y-1,X,T)) -dc*(D(Y,X,T)>=ec) -kn*N(Y,X,T) );                     
                end
            end
            for X=2: Xmax-1        
                N(1,X,T+1)= N(1,X,T) + dt*( dn*( D(1,X+1,T)+D(1,X-1,T)+D(2,X,T)) -dc*(D(1,X,T)>=ec) -kn*N(1,X,T) );  
                N(Xmax,X,T+1)= N(Xmax,X,T) + dt*( dn*( D(Xmax,X+1,T)+D(Xmax,X-1,T)+D(Xmax-1,X,T)) -dc*(D(Xmax,X,T)>=ec) -kn*N(Xmax,X,T)  );         
            end
            for Y=2: Xmax-1        
                N(Y,1,T+1)= N(Y,1,T) + dt*( dn*( D(Y+1,1,T)+D(Y-1,1,T)+D(Y,2,T)) -dc*(D(Y,1,T)>=ec) -kn*N(Y,1,T));  
                N(Y,Xmax,T+1)= N(Y,Xmax,T) + dt*( dn*( D(Y+1,Xmax,T)+D(Y-1,Xmax,T)+D(Y,Xmax-1,T)) -dc*(D(Y,Xmax,T)>=ec) -kn*N(Y,Xmax,T) );         
            end       
            N(1,1,T+1)= N(1,1,T) + dt*( dn*( D(1,2,T)+D(1,1,T)) -dc*(D(1,1,T)>=ec) -kn*N(1,1,T));  
            N(Xmax,1,T+1)= N(Xmax,1,T) + dt*( dn*( D(Xmax,2,T)+D(Xmax-1,1,T)) -dc*(D(Xmax,1,T)>=ec) -kn*N(Xmax,1,T) );         
            N(1,Xmax,T+1)= N(1,Xmax,T) + dt*( dn*( D(2,Xmax,T)+D(1,Xmax-1,T)) -dc*(D(1,Xmax,T)>=ec) -kn*N(1,Xmax,T) );
            N(Xmax,Xmax,T+1)= N(Xmax,Xmax,T) + dt*( dn*( D(Xmax,Xmax-1,T)+D(Xmax-1,Xmax,T)) -dc*(D(Xmax,Xmax,T)>=ec) -kn*N(Xmax,Xmax,T) ); 

            Temp(:,:,T)=(E(:,:,T) > N(:,:,T)).*(E(:,:,T)-N(:,:,T)-Sigma);
            A(:,:,T+1)=A(:,:,T)+dt*ea.*Temp(:,:,T).*(1-A(:,:,T));
        end
    end
end

Xmax=0.1;X=0:Xmax/20:Xmax;I=1;
figure('Position',[300 300 1500 1500]);colormap jet;
for alpha=1:alphaM%dn & dc
    dc=0.05+alpha/50;
    for beta=1:betaM%et and ec
        ec=0.03+beta/100;%0.05-0.1
        Yn=dn* 4 * X;
        Yc=dc*(X>=ec);
        subplot(alphaM,betaM,I);plot(X, Yc,'r',X, Yn,'k','LineWidth',2);
        set(gca,'YDir','normal');xlim([0 Xmax]);ylim([0 0.2]);
        title(strcat('dc=',num2str(dc),' ec=',num2str(ec))); 
        I=I+1;
    end
end