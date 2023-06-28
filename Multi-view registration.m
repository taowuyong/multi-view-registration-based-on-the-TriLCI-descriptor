%%comparison of two multi-view registration methods
%%shape growing-based method
list=dir(['E:\compile document\matlab\data\office22\','*.ply']);
kk=length(list);
% pr=0.49491739;
pr=0.006956074;
tic;
Mn=[];
for i=1:kk
    str= strcat ('E:\compile document\matlab\data\office22\', list(i).name);
    pcloud=pcread(str);
    PC=pcloud.Location;
    [n m]=size(PC);
    Mn=[Mn n];
end
[ma id]=max(Mn);
Rn=1:kk;
idid=find(Rn==id);
Rn(idid)=[];
str1= strcat ('E:\compile document\matlab\data\office22\', list(id).name);
pcloud1=pcread(str1);
PC1=pcloud1.Location;
eval(['PCt',num2str(id),'=PC1;']);
keypointcloud1 = pcdownsample(pcloud1,'gridAverage',7*pr);
keypoint1=keypointcloud1.Location;
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
for count1=1:100
    [n5 m5]=size(Rn);
    AdPC=[];
    Adcount=[];
    Adkeypoint=[];
    AdMV=[];
    AdMDV=[];
    for count=1:m5
        str2= strcat ('E:\compile document\matlab\data\office22\', list(Rn(count)).name);
        pcloud2=pcread(str2);
        PC2=pcloud2.Location;
%         plot3(PC1(:,1),PC1(:,2),PC1(:,3),'.b','MarkerSize',1);
%         hold on;
%         plot3(PC2(:,1),PC2(:,2),PC2(:,3),'.r','MarkerSize',1);
%         set(gca,'DataAspectRatio',[1 1 1]);
%         axis off
        [n m]=size(PC2);
        keypointcloud2 = pcdownsample(pcloud2,'gridAverage',7*pr);
        keypoint2=keypointcloud2.Location;
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
        if h0>3
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
                clear W
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
%         plot3(PC1(:,1),PC1(:,2),PC1(:,3),'.b','MarkerSize',1);
%         hold on;
%         plot3(PC2t(:,1),PC2t(:,2),PC2t(:,3),'.r','MarkerSize',1);
%         set(gca,'DataAspectRatio',[1 1 1]);
%         axis off
        [n4 m4]=size(PC1);
        [idx1 dist1]=knnsearch(PC1,PC2t,'k',1);
        iddist1=find(dist1<3*pr);
        overlap=length(iddist1)/min(n,n4);
        if overlap>0.5
            [Tg,PC2t,RMS,overlap1] = partialoverlapICP3(PC1,PC2t,pr);
            T=Tg*T;
            R=T(1:3,1:3);
            t=T(1:3,4)';
            eval(['PCt',num2str(Rn(count)),'=PC2t;']);
            AdPC=[AdPC;PC2t];
            Adcount=[Adcount count];
            keypoint2t=keypoint2*R'+ones(n2,1)*t;
            [idx2 dist2]=knnsearch(keypoint1,keypoint2t,'k',1);
            iddist2=find(dist2<7*pr);
            keypoint2t(iddist2,:)=[];
            Adkeypoint=[Adkeypoint;keypoint2t];
            MidMV2=[];
            for j=1:length(iddist2)
                idMV2=[3*iddist2(j)-2:3*iddist2(j)];
                MidMV2=[MidMV2 idMV2];
            end
            MV2(MidMV2,:)=[];
            [n6 m6]=size(MV2);
            MV2t=[];
            for j=1:n6/3
                MV2t=[MV2t;R*MV2(3*j-2:3*j,:)];
            end
            AdMV=[AdMV;MV2t];
            MDV2(iddist2,:)=[];
            AdMDV=[AdMDV;MDV2];
        end
    end
%     plot3(PC1(:,1),PC1(:,2),PC1(:,3),'.b','MarkerSize',1);
%     hold on;
%     plot3(AdPC(:,1),AdPC(:,2),AdPC(:,3),'.r','MarkerSize',1);
%     set(gca,'DataAspectRatio',[1 1 1]);
%     axis off
    Rn(Adcount)=[];
    PC1=[PC1;AdPC];
    keypoint1=[keypoint1;Adkeypoint];
    MV1=[MV1;AdMV];
    MDV1=[MDV1;AdMDV];
    if isempty(Rn) || isempty(Adcount)
        break;
    end
end
tt=toc;
for i=1:kk
    Idx=find(Rn==i);
    if ~isempty(Idx)
        continue;
    end
    eval(['PC','=PCt',num2str(i),';']);
    plot3(PC(:,1),PC(:,2),PC(:,3),'.','MarkerSize',1);
    set(gca,'DataAspectRatio',[1 1 1]);
    axis off
    hold on;
end









%%connected graph algorithm
list=dir(['E:\compile document\matlab\data\office22\','*.ply']);
kk=length(list);
% pr=0.49491739;
pr=0.006956074;
tic;
Mn=[];
for i=1:kk
    str= strcat ('E:\compile document\matlab\data\office22\', list(i).name);
    pcloud=pcread(str);
    PC=pcloud.Location;
    eval(['pcloud',num2str(i),'=pcloud;']);
    [n m]=size(PC);
    Mn=[Mn n];
end
[ma id]=max(Mn);
Rn=1:kk;
Rn(id)=[];
ReRn=id;
for count1=1:100
    eval(['pcloudd1','=pcloud',num2str(ReRn(count1)),';']);
    PCd1=pcloudd1.Location;
    [n2 m2]=size(Rn);
    ReRnC=[];
    Adcount=[];
    for count=1:m2
        eval(['pcloudd2','=pcloud',num2str(Rn(count)),';']);
        PCd2=pcloudd2.Location;
%         plot3(PCd1(:,1),PCd1(:,2),PCd1(:,3),'.b','MarkerSize',1);
%         hold on;
%         plot3(PCd2(:,1),PCd2(:,2),PCd2(:,3),'.r','MarkerSize',1);
%         set(gca,'DataAspectRatio',[1 1 1]);
%         axis off
        [ T,PCd2t] = coarseregistration( PCd1,PCd2,pcloudd1,pcloudd2,pr );
        [idx1 dist1]=knnsearch(PCd1,PCd2t,'k',1);
        iddist1=find(dist1<3*pr);
        [n3 m3]=size(PCd1);
        [n4 m4]=size(PCd2);
        overlap=length(iddist1)/min(n3,n4);
        if overlap>0.5
            [Tg,PCd2t,RMS,overlap1] = partialoverlapICP3(PCd1,PCd2t,pr);
            T=Tg*T;
            PCd2tCloud = pointCloud(PCd2t);
            eval(['pcloud',num2str(Rn(count)),'=PCd2tCloud;']);
            ReRnC=[ReRnC Rn(count)];
            Adcount=[Adcount count];
        end
    end
    Rn(Adcount)=[];
    ReRn=[ReRn ReRnC];
    [n5 m5]=size(ReRn);
    if count1==m5 || isempty(Rn)
        break;
    end
end
tt=toc;
for i=1:kk
    Idx=find(Rn==i);
    if ~isempty(Idx)
        continue;
    end
    eval(['pcloud','=pcloud',num2str(i),';']);
    PC=pcloud.Location;
    plot3(PC(:,1),PC(:,2),PC(:,3),'.','MarkerSize',1);
    set(gca,'DataAspectRatio',[1 1 1]);
    axis off
    hold on;
end
