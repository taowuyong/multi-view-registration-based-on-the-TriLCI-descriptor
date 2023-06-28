function [ T,PC2t] = coarseregistration( PC1,PC2,pcloud1,pcloud2,pr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[n m]=size(PC2);
keypointcloud1 = pcdownsample(pcloud1,'gridAverage',7*pr);
keypoint1=keypointcloud1.Location;
keypointcloud2 = pcdownsample(pcloud2,'gridAverage',7*pr);
keypoint2=keypointcloud2.Location;
RR=15*pr;
[n1 m1]=size(keypoint1);
[idx,dist]=rangesearch(PC1,keypoint1,RR);
MV1=[];
MDV1=[];
for i=1:n1
    KNN=PC1(idx{i},:);
    d=dist{i};
    [V] = LRF_TriLCI(KNN,RR,d,keypoint1(i,:));
    MV1=[MV1;V];
    [DV] = Tri_LCI(KNN,keypoint1(i,:),V,RR);
    MDV1=[MDV1;DV];
end
[n2 m2]=size(keypoint2);
[idx,dist]=rangesearch(PC2,keypoint2,RR);
MV2=[];
MDV2=[];
for i=1:n2
    KNN=PC2(idx{i},:);
    d=dist{i};
    [V] = LRF_TriLCI(KNN,RR,d,keypoint2(i,:));
    MV2=[MV2;V];
    [DV] = Tri_LCI(KNN,keypoint2(i,:),V,RR);
    MDV2=[MDV2;DV];
end
[idxx distt]=knnsearch(MDV1,MDV2,'k',2);
Mmatch=[];
for i=1:n2
    if distt(i,1)/distt(i,2)<=0.9
        match=[i idxx(i,1)];
        Mmatch=[Mmatch;match];
    end
end
[n3 m3]=size(Mmatch);
%3D transformation estimation method in this paper
%%false correspondence removal via two GCs
h0=0;
for i=1:n3
    Seed=Mmatch(i,:);
    p1=keypoint2(Seed(1),:);
    p2=keypoint1(Seed(2),:);
    V1=MV2(3*Seed(1)-2:3*Seed(1),:);
    V2=MV1(3*Seed(2)-2:3*Seed(2),:);
    R12=V2*V1';
    inlier=[];
    for j=1:n3
        p3=keypoint2(Mmatch(j,1),:);
        p4=keypoint1(Mmatch(j,2),:);
        d1=abs(norm(p1-p3)-norm(p2-p4));
        V3=MV2(3*Mmatch(j,1)-2:3*Mmatch(j,1),:);
        V4=MV1(3*Mmatch(j,2)-2:3*Mmatch(j,2),:);
        R34=V4*V3';
        d2=abs(norm(p1*R12'-p2)-norm(p3*R34'-p4));
        if d1<3*pr && d2<50*pr        %参数
            inlier=[inlier;Mmatch(j,:)];
        end
    end
    [h l]=size(inlier);
    if h>h0
        Cinlier=inlier;
        h0=h;
    end
end
%%transformation estimation based on robust estimation
if h0>3;
    A=keypoint2(Cinlier(:,1),:);
    Y=keypoint1(Cinlier(:,2),:);
    W0=ones(1,h0);
    R0=zeros(3,3);
    t0=zeros(1,3);
    for i=1:10000
        uA=[W0*A(:,1)/sum(W0) W0*A(:,2)/sum(W0) W0*A(:,3)/sum(W0)];
        uY=[W0*Y(:,1)/sum(W0) W0*Y(:,2)/sum(W0) W0*Y(:,3)/sum(W0)];
        H=zeros(3);
        for j=1:h0
            H=H+(A(j,:)-uA)'*W0(j)*(Y(j,:)-uY);
        end
        [U S V]=svd(H);
        D=diag([1 1 det(U*V')]);
        R=V*D*U';
        t=uY-uA*R';
        Res=Y-A*R'-ones(h0,1)*t;
        clear dd
        for j=1:h0
            dd(j)=norm(Res(j,:));
        end
        MAD=1.483*median(abs(dd-median(dd)*ones(1,h0)));
        Lb=(dd-median(dd)*ones(1,h0))/MAD;
        for j=1:h0
            k0=1.5;               %参数
            k1=2;
            if abs(Lb(j))<=k0;
                W(j)=1;
            elseif abs(Lb(j))>k0 && abs(Lb(j))<k1;
                W(j)=(k0/abs(Lb(j)))*((k1-abs(Lb(j)))/(k1-k0))*((k1-abs(Lb(j)))/(k1-k0));
            elseif abs(Lb(j))>=k1;
                W(j)=0;
            end
        end
        if norm(R-R0)+norm(t-t0)<0.00001
            break
        end
        R0=R;
        t0=t;
        W0=W;
    end
else
    R=eye(3);
    t=zeros(1,3);
end
PC2t=PC2*R'+ones(n,1)*t;
T=[R t';zeros(1,3) 1];
end

