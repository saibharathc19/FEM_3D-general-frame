function [GStiffnessBC,Kff,Kfc,Kcf,Kcc,X] =FrameStiffnessAssign(Kadd,DEG,NN,RestrainDoF,FreeDoF)


[t,d]=size(RestrainDoF);
[t,n]=size(FreeDoF);

GStiffnessBC=zeros(DEG*NN,DEG*NN);

Kcc=zeros(d,d);
Kfc=zeros(((DEG*NN)-d),d);
Kff=zeros(((DEG*NN)-d),((DEG*NN)-d));
Kcf=zeros(d,((DEG*NN)-d));

K=Kadd;
 
K(:,FreeDoF(1,:))=[];
K(FreeDoF(1,:),:)=[];
Kcc=K;
K=Kadd;
K(:,RestrainDoF(1,:))=[];
K(RestrainDoF(1,:),:)=[];
Kff=K;
K=Kadd;
K(:,RestrainDoF(1,:))=[];
K(FreeDoF(1,:),:)=[];
Kcf=K;
K=Kadd;
K(:,FreeDoF(1,:))=[];
K(RestrainDoF(1,:),:)=[];
Kfc=K;
Kff;
Kcc;
Kcf;
Kfc;
GStiffnessBC=[Kff Kfc;Kcf Kcc];
X=[Kff Kfc;Kcf Kcc];







