% Preprocessing 
clear all
SubID= {'BXY_LR','SZH_LR','LH_LR','ZHB_LR','WWJ_LR','QQ_LR','QL_LR','XXD_LR','ZZK_LR','ZJN_LR','WJK_LR','WLW_LR','ZZY_LR','GF_LR','YZY_LR'};
Time = linspace(0.25,1.25,50);
SR = 50;  % Temporal sampling rate
TFSpec = zeros(25,50,length(SubID));
for Si = 1:length(SubID)
    Dat = [];
    File = dir(fullfile(pwd,SubID{Si},strcat(SubID{Si},'_*.mat')));
    for Ri = 1:length(File)
        load(fullfile(pwd,SubID{Si},File(Ri).name));
        Dat = [Dat;List];
    end
    %% Condition
    TB = unique(Dat(:,1));
    Acc = zeros(length(TB),2);
    Corr_RT = Acc;
    RT = Acc;
%     %% Filter extreme RT values
%     MeanRT = mean(Dat(:,5));
%     SD = std(Dat(:,5));
%     Dat(Dat(:,5)>MeanRT+3*SD,5)=0;
%     Dat(Dat(:,5)<MeanRT-3*SD,5)=0;
       
    for Ti = 1:length(TB)
        Tmp = Dat(Dat(:,1)==TB(Ti),:);     
        Acc(Ti,1) = 1-sum(Tmp(Tmp(:,7)==0,4))/sum(Tmp(:,7)==0);
        RT(Ti,1) = mean(Tmp(Tmp(:,7)==0&Tmp(:,5)~=0,5));
        Corr_RT(Ti,1) = mean(Tmp(Tmp(:,7)==0&Tmp(:,4)==0&Tmp(:,5)~=0,5));
        Acc(Ti,2) = 1-sum(Tmp(Tmp(:,7)==1,4))/sum(Tmp(:,7)==1);
        RT(Ti,2) = mean(Tmp(Tmp(:,7)==1&Tmp(:,5)~=0,5));
        Corr_RT(Ti,2) = mean(Tmp(Tmp(:,7)==1&Tmp(:,4)==0&Tmp(:,5)~=0,5));
    end
    save (strcat(SubID{Si},'.mat'),'Acc','RT','Corr_RT');
    
end