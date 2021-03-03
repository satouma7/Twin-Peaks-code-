%Quantification of Bsh/Run-overlap (Fig.7g-i)
clear;filename='vps3_';

%Loading Bsh_control (X1,Y1)
Bsh1=readmatrix(strcat(filename,'Bsh_ctrl.csv'));
X1=Bsh1(:,1);Y1=Bsh1(:,2);N1=size(X1,1);
Width=ceil(max(X1)/10)*10;
Height=ceil(max(Y1)/100)*100;
Y1=Height-Y1;
Xave=round(sum(X1)/N1);Yave=round(sum(Y1)/N1);Ymin=min(Y1);
R1=zeros(N1,1);Rave=zeros(Yave,1);Rstd=zeros(Yave,1);

%Loading Run_control (V1,W1))
Run1=readmatrix(strcat(filename,'Run_ctrl.csv'));
V1=Run1(:,1);W1=Run1(:,2);M1=size(V1,1);
W1=Height-W1;

%Specifying origin (X0, Y0)
X0=(sum(X1)+sum(V1))/(N1+M1);
for Y0=1:Yave
    for I=1:N1
        R1(I)=sqrt((X1(I)-X0)^2+(Y1(I)-Y0)^2);  
    end
    Rave(Y0)=mean(R1);
    Rstd(Y0)=std(R1);
end
[R0,Y0]=min(Rstd);

%Loading Bsh_mutant (X2,Y2) and Run_mutant (V2, W2)
%Radius R2 for Bsh_mutant points
Bsh2=readmatrix(strcat(filename,'Bsh_mut.csv'));
X2=Bsh2(:,1);Y2=Bsh2(:,2);N2=size(X2,1);
Y2=Height-Y2;
R2=zeros(N2,1);

%Radius S2 for Run_mutant points
Run2=readmatrix(strcat(filename,'Run_mut.csv'));
V2=Run2(:,1);W2=Run2(:,2);M2=size(V2,1);
W2=Height-W2;
S2=zeros(M2,1);

figure('Position',[1500 500 1000 1000]);
subplot(2,2,3);plot(X2,Y2,'r*',V2,W2,'b+',X0,Y0,'ko');xlim([0 Width]);ylim([0 Height]);
title(strcat(filename,' BshRun-mutant'));

%angles of all mutant cells Bsh(X2,Y2) and Run(V2,W2)
Theta0=zeros(N2+M2,1);
Theta0(1:N2)=atan((Y2-Y0)./(X2-X0));
Theta0(N2+1:N2+M2)=atan((W2-Y0)./(V2-X0));

%Specifying the range of analysis for the control area
ThetaD=max(Theta0)-min(Theta0);
Tmax=-pi/2+ThetaD;

%Specifying the control coordinates for Bsh(X3, Y3) and Run(V3, W3)
Theta1=atan((Y1-Y0)./(X1-X0));J=1;
for I=1:N1
    if Theta1(I)<Tmax
        X3(J,1)=X1(I);
        Y3(J,1)=Y1(I);
        J=J+1;
    end
end
Theta2=atan((W1-Y0)./(V1-X0));K=1;
for I=1:M1
    if Theta2(I)<Tmax
        V3(K,1)=V1(I);
        W3(K,1)=W1(I);
        K=K+1;
    end
end

%angles of control cells Bsh(X3,Y3) and Run(V3,W3)
N3=size(X3,1);M3=size(V3,1);
R3=zeros(N3,1);S3=zeros(M3,1);
for I=1:N3
    R3(I)=sqrt((X3(I)-X0)^2+(Y3(I)-Y0)^2);  
end
for I=1:M3
    S3(I)=sqrt((V3(I)-X0)^2+(W3(I)-Y0)^2);  
end

subplot(2,2,1);plot(X1,Y1,'y*',V1,W1,'c+',X0,Y0,'ko',X3,Y3,'r*',V3,W3,'b+');xlim([0 Width]);ylim([0 Height]);
title(strcat(filename,' BshRun-ctrl'));

%Drawing the Bsh and Run domains in mutant and control areas
mutBsh=convhull(X2,Y2);mutRun=convhull(V2,W2);
ctrlBsh=convhull(X3,Y3);ctrlRun=convhull(V3,W3);
ctrlBshPoly=polyshape(X3(ctrlBsh),Y3(ctrlBsh));
ctrlRunPoly=polyshape(V3(ctrlRun),W3(ctrlRun));
mutBshPoly=polyshape(X2(mutBsh),Y2(mutBsh));
mutRunPoly=polyshape(V2(mutRun),W2(mutRun));
ctrlInter=intersect(ctrlBshPoly,ctrlRunPoly);
mutInter=intersect(mutBshPoly,mutRunPoly);

subplot(2,2,2);plot(ctrlInter,'FaceColor','m');title(['intersection:',num2str(area(ctrlInter)),' Bsh:',num2str(N2),' Run:',num2str(M2)]);hold on;
plot(X3(ctrlBsh),Y3(ctrlBsh),'r',X3,Y3,'r*',V3(ctrlRun),W3(ctrlRun),'b',V3,W3,'b+');hold off;

subplot(2,2,4);plot(mutInter,'FaceColor','m');title(['intersection:',num2str(area(mutInter)),' Bsh:',num2str(N3),' Run:',num2str(M2)]);hold on;
plot(X2(mutBsh),Y2(mutBsh),'r',X2,Y2,'r*',V2(mutRun),W2(mutRun),'b',V2,W2,'b+');hold off;