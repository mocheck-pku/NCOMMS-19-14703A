% Spectral analyses of the original time course
close all
clear all
SubID= {'BXY_LR','SZH_LR','LH_LR','ZHB_LR','WWJ_LR','QQ_LR','QL_LR','XXD_LR','ZZK_LR','ZJN_LR','WJK_LR','WLW_LR','ZZY_LR','GF_LR','YZY_LR'};
Time = linspace(0.25,1.25,50);
SR = 50;
NF = SR/2; % Nyquist frequency
FreqResol = 0.25; % Frequency domain resolution
Freqs = 0:FreqResol:NF;
NFreqs = length(Freqs); % Number of frequency points
GrpSpecs = zeros(length(SubID),length(Freqs),2);
GrpPhase = zeros(length(SubID),length(Freqs),2);
SigSource = 'Acc';
GSig = zeros(50,2);
List1 = [];
List2 = [];
for Si = 1:length(SubID)
    Sig = zeros(SR,2);
    load (fullfile(pwd,strcat(SubID{Si},'.mat')));
    eval(strcat('Sig = ',SigSource));
    GSig = GSig+Sig;
    List1 = [List1;Sig(:,1)'];
    List2 = [List2;Sig(:,2)'];
    [GrpSpecs(Si,:,1),GrpPhase(Si,:,1)] = FFT(Sig(:,1),0.25,50);
    [GrpSpecs(Si,:,2),GrpPhase(Si,:,2)] = FFT(Sig(:,2),0.25,50);    
end
figure (1)
hold on
if length(SubID)>=2
    subplot (2,1,1);
    imagesc(squeeze(GrpSpecs(:,:,1)));
    subplot (2,1,2);
    imagesc(squeeze(GrpSpecs(:,:,2)));
else    
    plot(Freqs,squeeze(GrpSpecs(1,:,1)),'b');
    plot(Freqs,squeeze(GrpSpecs(1,:,2)),'r');
end
% Ang = circ_dist(squeeze(GrpPhase(:,17,1)),squeeze(GrpPhase(:,17,2)));
for Si = 1:length(SubID)
    Ang(Si) = circ_dist(circ_mean(squeeze(GrpPhase(Si,17:18,1))),circ_mean(squeeze(GrpPhase(Si,17:18,2))));
end

figure(3)
circ_plot(Ang,'hist',[],12,false,true,'linewidth',2,'color','r');
save (strcat(SigSource,'_SpecAnalyses_Sub.mat'),'GrpSpecs');

figure(4)
hold on
xlim([min(Time),max(Time)]);
plot(Time,ft_preproc_polyremoval(GSig(:,1)'./length(SubID),1),'b');
plot(Time,ft_preproc_polyremoval(GSig(:,2)'./length(SubID),1),'r');
plot(Time,smooth(ft_preproc_polyremoval(GSig(:,1)'./length(SubID),1),3),'y');
plot(Time,smooth(ft_preproc_polyremoval(GSig(:,2)'./length(SubID),1),3),'g');

figure (5)
hold on
xlim([min(Time),max(Time)]);
plot(Time,GSig(:,1)'./length(SubID)-ft_preproc_polyremoval(GSig(:,1)'./length(SubID),1),'b');
plot(Time,GSig(:,2)'./length(SubID)-ft_preproc_polyremoval(GSig(:,2)'./length(SubID),1),'r');



MaxIt = 1000;
Distrib = zeros(MaxIt,2);
DisSpec1 = [];
DisSpec2 = [];

for It = 1:MaxIt
    for Si = 1:length(SubID)
        Sig = zeros(SR,2);
        load (fullfile(pwd,strcat(SubID{Si},'.mat')));
        eval(strcat('Sig = ',SigSource,';'));
        FC = fft(Sig(:,1));
        Mag = abs(FC);
        PS = angle(FC);
        PS = PS(randperm(length(PS)));
        IFC = Mag.*cos(PS)+1i*Mag.*sin(PS);
        TmpL = real(ifft(IFC));
        [FC,~,Freqs] = FFT(TmpL,0.25,50);
        TSSpec(1,Si,:) = FC;
        
        FC = fft(Sig(:,2));
        Mag = abs(FC);
        PS = angle(FC);
        PS = PS(randperm(length(PS)));
        IFC = Mag.*cos(PS)+1i*Mag.*sin(PS);
        TmpR = real(ifft(IFC));
        [FC,~,Freqs] = FFT(TmpR,0.25,50);
        TSSpec(2,Si,:) = FC;
    end
    X = mean(squeeze(TSSpec(1,:,:)),1);
    Y = mean(squeeze(TSSpec(2,:,:)),1);
    Distrib(It,1) = max(X);
    Distrib(It,2) = max(Y);
    DisSpec1 = [DisSpec1;X];
    DisSpec2 = [DisSpec2;Y];
end

 A= sort(Distrib(:,1),'descend');
 B= sort(Distrib(:,2),'descend');
 
 
 %%  Luo's correction
 FF1 = zeros(101,1);
 for Fi = 1:101
     T = sort(DisSpec1(:,Fi),'descend');
     FF1(Fi) = T(25);
 end
 FF2 = zeros(101,1);
 for Fi = 1:101
     T = sort(DisSpec2(:,Fi),'descend');
     FF2(Fi) = T(25);
 end
 
figure (2)
hold on
plot(Freqs,mean(squeeze(GrpSpecs(:,:,1)),1),'b');
plot(Freqs,mean(squeeze(GrpSpecs(:,:,2)),1),'r');
plot(Freqs,max(FF1)*ones(length(Freqs)),'b');
plot(Freqs,max(FF2)*ones(length(Freqs)),'r');
