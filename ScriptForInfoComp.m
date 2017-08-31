function ScriptForInfoComp
% close all
parameters=struct();
% parameters.tpoints = [-2 -1 0 1 2 3 4 5 6 7 8 10 12 14 16 20 24] ;
parameters.tpoints = [-2 0 4 8 12 16 20 24 28 32 36 40 44 48 52] ;
parameters.iterations = 20;
parameters.percentageofsamplesize = 0.7:0.05:0.95;
parameters.number_of_PofS =50;
parameters.confidencepercent = 90;
colorlist = {[0 0.5 0], [0.4 0.4 0], [0.5 0 0],    [0.3 0 0.2],    [0 0.3 0.2]};
cycle=0;
for scalartype =    {'fc',      'diff' ,    'abs',          'int',          'fcint'}
% for scalartype =    {'int',          'fcint'}
cycle=cycle+1;
% for scalartype = {'int','fcint'} 
parameters.scalartype = char(scalartype);
parameters.color = colorlist{cycle};
if regexp(char(scalartype),'int')
    parameters.tpoints = [0 4 8 12 16 20 24 28 32 36 40 44 48 52] ;%for int and fcint you cannot do <0 time points because all values will be nan
end
    
%% quantify mutual information
% fileDateName = '2017_04_17 plate exp2_tracking_export.mat';
% MutualInformationSelimkhanovMethodIterations_newparfor_OneD(fileDateName,parameters)
% fileDateName = '2017_04_17 plate exp3_tracking_export.mat';
% MutualInformationSelimkhanovMethodIterations_newparfor_OneD(fileDateName,parameters)
% fileDateName = '2017_04_17 plate exp4_tracking_export.mat';
% MutualInformationSelimkhanovMethodIterations_newparfor_OneD(fileDateName,parameters)
% % % % % fileDateName = '2014_09_30 plate exp1_tracking_export.mat';
% % % % MutualInformationSelimkhanovMethodIterations_newparfor_OneD(fileDateName,parameters)
% % % % fileDateName = '2017_01_30 plate exp2_tracking_export.mat';
% % % % MutualInformationSelimkhanovMethodIterations_newparfor_OneD(fileDateName,parameters)
% % % % fileDateName = '2017_03_13 plate exp4_tracking_export.mat';
% % % % MutualInformationSelimkhanovMethodIterations_newparfor_OneD(fileDateName,parameters)
% % % % fileDateName = '2017_03_13 plate exp3_tracking_export.mat';
% % % % MutualInformationSelimkhanovMethodIterations_newparfor_OneD(fileDateName,parameters)
% % % % close all
% % % 
% % % % fileDateName = '2017_04_17 plate exp2_tracking_export.mat';
% % % % dispName = 'endogenous (day1)';
% % % % plotInformation(fileDateName,dispName,parameters)
% % % % fileDateName = '2017_04_17 plate exp3_tracking_export.mat';
% % % % dispName = 'endogenous (day2)';
% % % % plotInformation(fileDateName,dispName,parameters)
% % % % fileDateName = '2017_04_17 plate exp4_tracking_export.mat';
% % % % dispName = 'endogenous (day3)';
% % % plotInformation(fileDateName,dispName,parameters)
% % 

%% bias correct
% % determine the maximum information using biasCorrection algorithm
% fileDateName = '2017_04_17 plate exp2_tracking_export.mat';
% dispName = 'endogenous, day1';
% biasCorrectInfo_andFindMax_OneD(fileDateName,dispName,parameters)
% fileDateName = '2017_04_17 plate exp3_tracking_export.mat';
% dispName = 'endogenous, day2';
% biasCorrectInfo_andFindMax_OneD(fileDateName,dispName,parameters)
% fileDateName = '2017_04_17 plate exp4_tracking_export.mat';
% dispName = 'endogenous, day3';
% biasCorrectInfo_andFindMax_OneD(fileDateName,dispName,parameters)


%% make the plots
numberOfPlots = 3;
plotNum =1;
    fileDateName = '2017_04_17 plate exp2_tracking_export.mat';
    dispName = 'endogenous, day1';
    plotInformationBiasCorrectionAlt_OneD(fileDateName,dispName,parameters,numberOfPlots,plotNum)
plotNum=2;
    fileDateName = '2017_04_17 plate exp3_tracking_export.mat';
    dispName = 'endogenous, day2';
    plotInformationBiasCorrectionAlt_OneD(fileDateName,dispName,parameters,numberOfPlots,plotNum)
plotNum=3;
    fileDateName = '2017_04_17 plate exp4_tracking_export.mat';
    dispName = 'endogenous, day3';
    plotInformationBiasCorrectionAlt_OneD(fileDateName,dispName,parameters,numberOfPlots,plotNum)
% fileDateName = '2014_09_30 plate exp1_tracking_export.mat';
% dispName = 'exogenous, day2';
% plotInformationBiasCorrectionAlt(fileDateName,dispName,parameters)
% fileDateName = '2017_01_30 plate exp2_tracking_export.mat';
% dispName = 'rep1 endogenous, day2';
% plotInformationBiasCorrectionAlt(fileDateName,dispName,parameters)
% fileDateName = '2017_03_13 plate exp4_tracking_export.mat';
% dispName = 'endogenous-cmv, day1';
% plotInformationBiasCorrectionAlt(fileDateName,dispName,parameters)
% fileDateName = '2017_03_13 plate exp3_tracking_export.mat';
% dispName = 'rep1 endogenous, day1';
% plotInformationBiasCorrectionAlt(fileDateName,dispName,parameters)
end



end