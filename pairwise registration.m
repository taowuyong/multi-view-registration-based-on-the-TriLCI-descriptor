%%pairwise registration 
list=dir(['D:\compile document\matlab\data\Registration parasaurolophus_high\','*.ply']);
kk=length(list);
pr=0.49491739;
for i=1:kk-1
    str1= strcat ('D:\compile document\matlab\data\Registration parasaurolophus_high\', list(i).name);
    pcloud1=pcread(str1);
    PC1=pcloud1.Location;
    for j=i+1:kk
        str2= strcat ('D:\compile document\matlab\data\Registration parasaurolophus_high\', list(j).name);
        pcloud2=pcread(str2);
        PC2=pcloud2.Location;
%         plot3(PC1(:,1),PC1(:,2),PC1(:,3),'.b','MarkerSize',1);
%         hold on;
%         plot3(PC2(:,1),PC2(:,2),PC2(:,3),'.r','MarkerSize',1);
%         set(gca,'DataAspectRatio',[1 1 1]);
%         axis off
        [ T,PC2t] = coarseregistration( PC1,PC2,pcloud1,pcloud2,pr );
%         plot3(PC1(:,1),PC1(:,2),PC1(:,3),'.b','MarkerSize',1);
%         hold on;
%         plot3(PC2t(:,1),PC2t(:,2),PC2t(:,3),'.r','MarkerSize',1);
%         set(gca,'DataAspectRatio',[1 1 1]);
%         axis off
        name1=list(i).name;
        name2=list(j).name;
        strout= strcat ('D:\文本\多视配准\TOLDI registration result T-rex\',name2(1:end-4),'-',name1(1:end-4),'.txt');
        fid=fopen(strout,'wt');
        fprintf(fid,'%f %f %f %f\n',T');
        fclose(fid);
    end
end