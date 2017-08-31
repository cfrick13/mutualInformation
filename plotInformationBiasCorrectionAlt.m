function plotInformationBiasCorrectionAlt(fileDateName,dispName,parameters,numberOfPlots,plotNum)
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
addpath([parentdir 'Colormaps']);

%load the file
    cd(mfiledir)
    flist = dir(strcat('biasCorrMaxInfo*',datename,'*',num2str(parameters.iterations),'*.mat'));
    fname = char(flist.name);
    load(fname);
    
    
    
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
    
sortAssem1 = cellfun(@(x) sort(x),assembledFeatures,'UniformOutput',0);
timeStrArray = cellfun(@(x) [num2str(tvecz(x(1))) ' x ' num2str(tvecz(x(2)))],sortAssem1,'UniformOutput',0);


sortAssem1StrArray = cellfun(@(x) [num2str(x(1)) num2str(x(2))],sortAssem1,'UniformOutput',0);
[uniqueAssemStrArray,uidx] = unique(sortAssem1StrArray);

t1array = cellfun(@(x) x(1),sortAssem1,'UniformOutput',1);
t2array = cellfun(@(x) x(2),sortAssem1,'UniformOutput',1);
tvecz1array = tvecz(t1array(uidx));
tvecz2array = tvecz(t2array(uidx));
uTvals = unique(tvecz1array);

infoimage = nan(length(uTvals),length(uTvals));
infoimagestd = nan(length(uTvals),length(uTvals));
infoimageneg = infoimagestd;
infoimagepos = infoimagestd;
for i = 1:length(uidx)
    idx = uidx(i);
    infoval = imax(idx);
    ebstdval = ebarstd(idx);
    ebnegval = ebarneg(idx);
    ebposval = ebarpos(idx);
    t1val = tvecz1array(i);
    t2val = tvecz2array(i);
    a = find(uTvals == t1val);
    b = find(uTvals == t2val);
    infoimage(a,b) = infoval;
    infoimagestd(a,b) = ebstdval;
    infoimageneg(a,b) = ebnegval;
    infoimagepos(a,b) = ebposval;
%     infoimage(b,a) = infoval;
end

zeroIDX = find(uTvals==0);
ivec = [-1 1 3 7 9 11]+zeroIDX;
cutoff=4;
cutoffidx = ivec<(zeroIDX+cutoff);
cmap1 = colormap(winter(sum(cutoffidx)+1));
cmap2 = colormap(autumn(sum(~cutoffidx)*2));
f=  figure(202);
f.Color = 'w';
f.Position = [158 537 1402 441];
    s = subplot(1,numberOfPlots,plotNum);


    cyc1=0;
    cyc2=0;
    for i = ivec
        if i<zeroIDX+cutoff
            cyc1=cyc1+1;
            ecolor = cmap1(cyc1,:);
            mark = 's';
        else
            cyc2=cyc2+1;
            ecolor = cmap2(cyc2,:);
            mark = 'o';
        end
        infovec = infoimage(i,:);
        stdvec = infoimagestd(i,:);
        infoneg = infoimageneg(i,:);
        infopos = infoimagepos(i,:);
        xvec = uTvals;
        e = errorbar(xvec,infovec,infovec-infoneg,infopos-infovec);hold on
        e.Marker = mark;
        e.MarkerSize = 8;
        e.MarkerFaceColor = ecolor;
        e.DisplayName = num2str(uTvals(i));
        e.Color = ecolor./2;
        e.LineWidth =1;
    end
    s.YLim =[0 1.5];
    s.XLim =[-20 80];
    s.Color = [0.9 0.9 0.9];
    s.XGrid = 'on';
    s.YGrid = 'on';
    s.Title.String = {'2D infomation'; dispName};
    s.XLabel.String = 'timepoint 2, minutes';
    s.YLabel.String = 'channel capacity, bits';
    l = legend('show');
    l.Title.String = {'timepoint 1', 'minutes'};
    l.Color = 'w';
    l.Location = 'southeast';
    l.FontSize = 8;
    
    


xvec = 1:length(uidx);


    f=figure(44);
    subplot(2,numberOfPlots,plotNum);

%  
%     for tcombo = 1:length(tcombovec)
%     end
    tidx = tvecz1array==0;
        bplot = bar(xvec(tidx),imax(tidx));hold on
        h=gca;
        h.XTick = 1:length(uidx);
        h.XTickLabel = timeStrArray(uidx);
        h.XTickLabelRotation =45;
        eplot = errorbar(xvec(tidx),imax(tidx),ebarstd(tidx));
        eplot.LineStyle = 'none';

    s = subplot(2,numberOfPlots,plotNum+numberOfPlots);
        imagesc(infoimage);
        s.CLim = [0 1.5];
        colormap('viridis')
        s.YDir = 'normal';
        s.XTick = 1:length(uTvals);
        s.XTickLabel = uTvals;
        s.YTick = 1:length(uTvals);
        s.YTickLabel = uTvals;
        

stophere=1;
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