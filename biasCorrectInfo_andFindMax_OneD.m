function biasCorrectInfo_andFindMax_OneD(fileDateName,dispName,parameters)
parameterz = parameters;
datename = fileDateName(1:21);
%determine the parent folder and directories
    mdir = mfilename('fullpath');
        [~,b] = regexp(mdir,'Tracking\w*/');
            if isempty(b)
                [~,b] = regexp(mdir,'Tracking\w*\');
            end
    parentdir = mdir(1:b); %specifies folder in which all analysis is being done
    exportdir = strcat(parentdir,'Export'); %specifies where data is exported
    cd(exportdir);
    loaddir = strcat(parentdir,'Export'); %specifies where data is exported

    [~,b] = regexp(mdir,'/');
            if isempty(b)
                [~,b] = regexp(mdir,'\');
            end

    mfiledir =mdir(1:b(end)); %specifies location of matlab function file    
mfiledirz = mfiledir;
%load the file
    cd(mfiledir)
    flist = dir(strcat(parameters.scalartype,'highINFORMATIONzSELIMKHANOV*',datename,'*',num2str(parameters.iterations),'*.mat'));
    fname = char(flist.name);
    disp(fname)
    load(fname);
    
    %load metadata associated with the experiment (requires manual input if there is ambiguity)
    FileName = fileDateName;
    [a,~] = regexp(FileName,'_tracking');
    datequery = strcat(FileName(1:a-1),'*metaData.mat');
    cd(loaddir)
    filelist = dir(datequery);
    if length({filelist.name}) ==1
        metaData = load(char(filelist.name));
    else
        filename = uigetfile();
        metaData = load(filename);
    end
%determine the timeVector from metaData [dim 1 is scene#, dim 2 is each time frame]
    timeVec = metaData.timeVec;
    tvec = timeVec(1,:) ;
    tvecz = round(tvec-tvec(stimulationFrame+1));
%     selectedFeatures = [stimulationFrame-2:stimulationFrame+30];
%     selectedFeatures = [stimulationFrame-1:stimulationFrame+11];
    tm = round(tvec(selectedFeaturesVec)-tvec(stimulationFrame+1));
    
    %plot for figure
%     infoInUnitsOfBitsMatrix(iter,inputcyclenum,samplingSizecycle,selectedFeatures)
    
    parameters.confidencepercent = parameterz.confidencepercent;
    sampleSizeVec = parameters.percentageofsamplesize;
    infoBiasCorr = zeros(size(infoInUnitsOfBitsMatrix,2),size(infoInUnitsOfBitsMatrix,4));
    infoBiasCorrSE = zeros(size(infoInUnitsOfBitsMatrix,2),size(infoInUnitsOfBitsMatrix,4));
    infoBiasCorrPOS = infoBiasCorrSE;
    infoBiasCorrNEG = infoBiasCorrSE;
%     for d4 = 1:size(infoInUnitsOfBitsMatrix,4) %cycle through all "responses"
    parfor d4 = 1:size(infoInUnitsOfBitsMatrix,4) %cycle through all "responses"
        infoForResponse = squeeze(infoInUnitsOfBitsMatrix(:,:,:,d4)); %dimensions are (iterations,inputcyclenum,samplingSizecycle)
        ibias = zeros(1,size(infoInUnitsOfBitsMatrix,2));
        ibiasSE = zeros(1,size(infoInUnitsOfBitsMatrix,2));
        ibiasPOS = ibiasSE;
        ibiasNEG = ibiasSE;
        for d2 = 1:size(infoInUnitsOfBitsMatrix,2) %cycle through all signal input distributions
            infoForResponseGivenInputDistribution = squeeze(infoForResponse(:,d2,:)); %dimesions are (iterations, samplingSizecycle)
            si1 = size(infoInUnitsOfBitsMatrix,1);
            yint = zeros(1,si1);
            for d1 = 1:size(infoInUnitsOfBitsMatrix,1)
                x = sampleSizeVec.^-1;
                y = infoForResponseGivenInputDistribution(d1,:);
                mdl = fitlm(x,y);
                estms = mdl.Coefficients.Estimate;
                yint(d1) = estms(1);
            end
            yintmean = nanmean(yint);
            yintSTD = nanstd(yint);
            yintnan = yint(~isnan(yint));
            yintPOS = prctile(yintnan,50+(parameters.confidencepercent./2));
            yintNEG = prctile(yintnan,50-(parameters.confidencepercent./2));
            
            ibias(d2) = yintmean;
            ibiasSE(d2) = yintSTD;
            ibiasPOS(d2) = yintPOS;
            ibiasNEG(d2) = yintNEG;
        end
        infoBiasCorr(:,d4) = ibias; %store yint
        infoBiasCorrSE(:,d4) = ibiasSE;
        infoBiasCorrPOS(:,d4) = ibiasPOS;
        infoBiasCorrNEG(:,d4) = ibiasNEG;
    end

    ebarstd = zeros(1,size(infoInUnitsOfBitsMatrix,4));
    ebarPOS = ebarstd;
    ebarNEG = ebarstd;
    imax = zeros(1,size(infoInUnitsOfBitsMatrix,4));
    for d4 = 1:size(infoInUnitsOfBitsMatrix,4)
        i3 = squeeze(infoBiasCorr(:,d4));
        i3SE = squeeze(infoBiasCorrSE(:,d4));
        i3POS = squeeze(infoBiasCorrPOS(:,d4));
        i3NEG = squeeze(infoBiasCorrNEG(:,d4));
        %determine input distribution at which info is max
        [~,idx] = max(i3,[],1);
        ebarstd(d4) = i3SE(idx);
        ebarPOS(d4) = i3POS(idx);
        ebarNEG(d4) = i3NEG(idx);
        imax(d4) = i3(idx);
    end
    numberOfCells = length([PRofScellarray{:,1}]);
    nstr = ['n = ' num2str(numberOfCells)];

    
    if strcmp(datename,'2014_09_30 plate exp1')
%         lc = [0 0.4980 0];
        lc = [0 0.4470 0.7410];
    elseif strcmp(datename,'2017_01_30 plate exp2')
        lc = [0.6353 0.0784 0.1843];
    elseif strcmp(datename,'2017_03_13 plate exp4')
%         lc = [0 0.4470 0.7410];
        lc = [0 0.4980 0];
    elseif strcmp(datename,'2017_03_13 plate exp3')
%         lc = [0 0.4470 0.7410];
        lc = [0.4 0.0784 0.25];
    elseif strcmp(datename,'2017_04_17 plate exp2')
        lc = [0.5 0.1784 0.15];
    elseif strcmp(datename,'2017_04_17 plate exp3')
        lc = [0.0 0.7784 0.15];
    elseif strcmp(datename,'2017_04_17 plate exp4')
        lc = [0.7 0.0784 0.75];
    end
    
    stophere=1;
sortAssem1 = cellfun(@(x) sort(x),assembledFeatures,'UniformOutput',0);
% timeStrArray = cellfun(@(x) [num2str(tvecz(x(1))) ' x ' num2str(tvecz(x(2)))],sortAssem1,'UniformOutput',0);
% sortAssem1StrArray = cellfun(@(x) [num2str(x(1)) num2str(x(2))],sortAssem1,'UniformOutput',0);
% [uniqueAssemStrArray,uidx] = unique(sortAssem1StrArray);
uniqueAssemStrArray = sortAssem1;

    cd(mfiledirz)
    save(strcat(parameterz.scalartype,'biasCorrMaxInfo',datename,'-',num2str(dimensions),'-',num2str(iterations),'.mat'),'-v7.3');

%     f=figure(44);
%         b = bar(1:length(uidx),imax(uidx));hold on
%         h=gca;
%         h.XTick = 1:length(uidx);
%         h.XTickLabel = timeStrArray(uidx);
%         h.XTickLabelRotation =45;
%         e = errorbar(1:length(uidx),imax(uidx),ebarstd(uidx));


%     s = subplot(2,1,1);
%     %level
%         lsm = length(selectedFeaturesVec);
%         tidx = 1:lsm;
%         idx = 1:lsm;
%         e = errorbar(tm(tidx),imax(idx),ebarstd(idx));hold on
%         e.DisplayName = ['level'];
%         e.LineWidth = 1.5;
%         e.Color = [0.6353 0.0784 0.1843];
%     %fold change
%         lsm = length(selectedFeaturesVec);
%         tidx = 1:lsm;
%         idx = [1:lsm]+lsm-1;
%         e = errorbar(tm(tidx),imax(idx),ebarstd(idx));hold on
%         e.DisplayName = ['fold-change'];
%         e.LineWidth = 1.5;
%         e.Color = [0 0.4471 0.7412];
%     %graph labels
%         ylim([0 1.4]);
%         s.Title.String = {'channel capacity of NG-Smad3 response'; [dispName ', ' nstr]};
%         s.XLabel.String = 'minutes after Tgfbeta addition';
%         s.YLabel.String = {'Channel Capacity, bits';'(with bias correction)';'+/- STDev'};
%         e.LineWidth = 1.5;
%         xlim([min(tm) max(tm)]);
%         ylim([0 1.4]);
%         l = legend('show');
%         f.Position = [1422 343 785 459];


    
%         f=gcf;
%         for h=f.Children'
%             if strcmp(h.Type,'legend')
%                 h.delete;
%             end
%         end
%         for h=f.Children'
%         l = legend(h,'show');
%         l.EdgeColor = 'w';
%         h.Units = 'pixels';
%         h.Position = [144.0000 64 528.0000 348.0250];
%         h.FontSize = 12;
%         h.FontName = 'helvetica';
%         h.Color = [0.97 0.97 0.97];
%         h.GridLineStyle = '--';
%         h.GridColor = 'k';
%         h.Box = 'off';
%         h.XGrid = 'on';
%         h.YGrid = 'on';
%         h.LineWidth = 1.5;
%         h.XColor = 'k';
%         h.YColor = 'k';
%         h.XLim = [-40 120];
%         h.Units = 'normalized';
%         end
%         f.Color = [1 1 1];
%         
%         %change directory to image export folder
%         cd(mfiledir)
%         savedirname = 'figureexport';
%         dirlist = dir(savedirname);
%         if isempty(dirlist)
%             mkdir(savedirname)
%         end
%         olddir = pwd;
%         cd(savedirname)
%         
%         filestr = fileDateName(1:21);
%         filename = [filestr ' biasCorrected-compare FC vs level' '.fig'];
%         saveas(f,filename);
%         filename = [filestr ' biasCorrected-compare FC vs level' '.png'];
%         saveas(f,filename);
%         
%         cd(olddir)
% close all


    
    
%     figure(33)
%     info = mean(squeeze(infoInUnitsOfBitsMatrix(:,1,:,1)),1);
%     plot((1./percentageofsamplesize),info);hold on
%     info = mean(squeeze(infoInUnitsOfBitsMatrix(:,1,:,2)),1);
%     plot((1./percentageofsamplesize),info);
% 
%     info = mean(squeeze(infoInUnitsOfBitsMatrix(:,1,:,2)),1);
%     infostd = nanstd(squeeze(infoInUnitsOfBitsMatrix(:,1,:,2)),1);
%     errorbar((1./percentageofsamplesize),info,infostd,'LineStyle','none','Color','k')
%     info = mean(squeeze(infoInUnitsOfBitsMatrix(:,1,:,1)),1);
%     infostd = nanstd(squeeze(infoInUnitsOfBitsMatrix(:,1,:,1)),1);
%     errorbar((1./percentageofsamplesize),info,infostd,'LineStyle','none','Color','k')
% 
%     ylim([0 2])
%     xlim([0 2])
%     xlabel('1/samplesize')
%     ylabel('estimated mutual information')
%     title('jacknife sampling without replacement to determine bias due to sample size')
end