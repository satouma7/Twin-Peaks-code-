%2 dimension 4 variables Phase Diagram (Supplementary Fig2)
%EGF (E), Notch signal (N), Delta (D), AS-C (A)
%Hill function
dn=0.25;dc=0.25;
nt=1;kt=1;ktnt=kt^nt;%nt: Hill's co-efficient@kt:activation co-efficient  for trans-activation
nc=5;kc=0.1;kcnc=kc^nc;%nc: Hill's co-efficient@kc:activation co-efficient for cis-inhibition

kd=1;ae=1;ad=1;
ea=100;de=1;ke=1;kn=1;Sigma=0;
dx = 2;dx2=dx * dx;dt=0.01;Tmax = 1000;
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

figure('Position',[300 300 1500 1500]);colormap jet;
I=1;alphaM=6;betaM=5;
for alpha=1:alphaM%dn & dc
    nc=1+alpha;
    for beta=1:betaM%et and ec
        kc=0.06+beta/100;%0.05-0.1
        kcnc=kc^nc;
        for T=1: Tmax-1
            Afront=A(13,13,T);       
            if Afront>0.5
                subplot(alphaM,betaM,I);
                plot(Xsteps, A(13,:,T),'k',Xsteps, N(13,:,T)*30,'r',Xsteps,zeros(1,Xmax),'k--','LineWidth',2);
                %plot(Xsteps, A(13,:,T),'b',Xsteps, D(13,:,T)*10,'g',Xsteps, N(13,:,T)*30,'r',Xsteps, E(13,:,T)*10,'k',Xsteps,zeros(1,Xmax),'k--','LineWidth',2);
                set(gca,'YDir','normal');xlim([1 Xmax]);ylim([-0.1 1]);
                title(strcat('nc=',num2str(nc),' kc=',num2str(kc)));
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
    end
end

Xmax=0.1;X=0:Xmax/20:Xmax;I=1;
figure('Position',[300 300 1500 1500]);colormap jet;
for alpha=1:alphaM%dn & dc
    nc=1+alpha;
    for beta=1:betaM%et and ec
        kc=0.06+beta/100;%0.05-0.1
        kcnc=kc^nc;SigmaDNt=(4*X).^nt;DNc=X.^nc;
        Yn=dn* SigmaDNt./(ktnt+SigmaDNt);
        Yc=dc*DNc./(kcnc+DNc);
        subplot(alphaM,betaM,I);plot(X, Yc,'r',X, Yn,'k','LineWidth',2);
        set(gca,'YDir','normal');xlim([0 Xmax]);ylim([0 0.2]);
        title(strcat('nc=',num2str(nc),' kc=',num2str(kc)));   
        I=I+1;
    end
end
