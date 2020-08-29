%% Clean workspace and load data
close all
clear all
clc

load('FinsAblationEscape.mat')

%% Initiate matrix

nbStimulus = 10;
nbIdx =84;
% rows -> stimulus
% column -> fish
FishID= zeros(nbStimulus,nbIdx)*NaN;
FishTrial= zeros(nbStimulus,nbIdx)*NaN;
FishGeno= zeros(nbStimulus,nbIdx)*NaN;
Roll= zeros(nbStimulus,nbIdx)*NaN;
Cbend= zeros(nbStimulus,nbIdx)*NaN;
Latency= zeros(nbStimulus,nbIdx)*NaN;
BendTime= zeros(nbStimulus,nbIdx)*NaN;
BeatFreq= zeros(nbStimulus,nbIdx)*NaN;
BeatTime= zeros(nbStimulus,nbIdx)*NaN;
BoutDuration = zeros(nbStimulus,nbIdx)*NaN;
NumberOfOscillations= zeros(nbStimulus,nbIdx)*NaN;
Distance= zeros(nbStimulus,nbIdx)*NaN;
Speed= zeros(nbStimulus,nbIdx)*NaN;
mTBF= zeros(nbStimulus,nbIdx)*NaN;
TBF2 = zeros (nbStimulus, nbIdx)*NaN;

Counterbend= zeros(nbStimulus,nbIdx)*NaN;
RollIntegral  =  zeros(nbStimulus,nbIdx)*NaN;
RollSum=  zeros(nbStimulus,nbIdx)*NaN;
LongRoll = zeros(nbStimulus,nbIdx)*NaN;
RolloverOnset = zeros(nbStimulus,nbIdx)*NaN;
RollingOnset= zeros(nbStimulus,nbIdx)*NaN;
RolloverNb = zeros(nbStimulus,nbIdx)*NaN;

EscapeType = zeros(nbStimulus,nbIdx)*NaN;
EscapeRatio = zeros(nbStimulus,nbIdx);

%% Calculation for each trial
datasetPerBout3 = datasetPerBout([datasetPerBout(:).NumberOfOscillations]>1);

datasetPerBout1 = datasetPerBout3(([datasetPerBout3(:).BoutStart]>70) & ([datasetPerBout3(:).BoutStart]<100)); %(50ms)
%%
RolloverSeq = [1,1,1,1];%%10 ms
RollingSeq = (1);
for i = 1:length(datasetPerBout1)
            % Find New BoutStart
            TailAngle=57.2958*[datasetPerBout1(i).TailAngle_smoothed];
            
            % find the Start point by calculating diff(TailAngle)>5, find the transient start point
            diff_Angle=find(abs(diff(TailAngle))>3);   
            if isempty(diff_Angle)
                diff_Angle=1;
            end           
            % Set the threshold of movement at 5 degree for all movement
            New_MovePosition=find(abs(TailAngle)>5);           
            if isempty(New_MovePosition)
                New_MovePosition=1;
            end      
            % datasetPerBout New Parameter
            datasetPerBout1(i).Escape= 1;
            datasetPerBout1(i).New_Timing=[diff_Angle(1):1:New_MovePosition(end)];
            datasetPerBout1(i). New_BoutStart = diff_Angle(1)+ datasetPerBout1(i).BoutStart;
            datasetPerBout1(i). Latency = ((diff_Angle(1)+ datasetPerBout1(i).BoutStart)-78)/400*1000;
            datasetPerBout1(i).Cbend = abs(datasetPerBout1(i).AmplitudeOfAllBends(1));
            datasetPerBout1(i).BendTime = (datasetPerBout1(i).Bend_Timing(1)-diff_Angle(1)+1)*2.5;
            datasetPerBout1(i).BeatTime = (datasetPerBout1(i).Bend_Timing(2)-datasetPerBout1(i).Bend_Timing(1))*2.5;
            datasetPerBout1(i).BeatFreq = 1000/((datasetPerBout1(i).Bend_Timing(2)-datasetPerBout1(i).Bend_Timing(1))*2.5*2);%TBF1
            datasetPerBout1(i).Counterbend = abs(datasetPerBout1(i).AmplitudeOfAllBends(2));
            
            if datasetPerBout1(i).Cbend <60
                 datasetPerBout1(i).EscapeType  = 0; %S escape
             else
                 datasetPerBout1(i).EscapeType = 1; % C escape
             end
            datasetPerBout1(i).RollIntegral = mean (datasetPerBout1(i).rolloverProbability);  
           % datasetPerBout1(i).RollSum = sum (datasetPerBout1(i).rolloverProbability([datasetPerBout1(i).rolloverProbability]>0.8))/length ((datasetPerBout1(i).rolloverProbability([1:end])>0.8));  
            datasetPerBout1(i).RollSum =  2.5 * sum ([datasetPerBout1(i).rolloverProbability]>0.80);  
           
            datasetPerBout1(i).TBF2 = datasetPerBout1(i).InstantaneousTBF(2);  
            datasetPerBout1(i).mTBF = (datasetPerBout1(i).NumberOfOscillations)/ (datasetPerBout1(i).BoutDuration);  
            datasetPerBout1(i).rollover2 = ([datasetPerBout1(i).rolloverProbability]>0.8);
            %DL rolling start
            datasetPerBout1(i).RolloverOnset = NaN;            
            datasetPerBout1(i).RolloverStop = NaN;
             if max([datasetPerBout1(i).rollover2])<1
                datasetPerBout1(i).DLRoll = 0;
                datasetPerBout1(i).RollStart = NaN;              
                datasetPerBout1(i).LongRoll = 0; 
                datasetPerBout1(i).LongRollOnset = NaN;
                datasetPerBout1(i).DLRolling = NaN;  
                datasetPerBout1(i).RolloverTime = NaN;

             else
               RollingSequence = strfind([datasetPerBout1(i).rollover2],RollingSeq); 
               RolloverSequence = strfind([datasetPerBout1(i).rollover2],RolloverSeq);
                  if  isempty (RollingSequence)
               datasetPerBout1(i).DLRolling = NaN;
               datasetPerBout1(i).DLRoll = 0;
               datasetPerBout1(i).RollStart = NaN;  
               datasetPerBout1(i).LongRoll = 0; 
               datasetPerBout1(i).LongRollOnset = NaN; 
               datasetPerBout1(i).RolloverTime = NaN;
                  else
               datasetPerBout1(i).DLRolling = RollingSequence;    
               datasetPerBout1(i).RollStart = 2.5* RollingSequence(1);         
                      if isempty (RolloverSequence)
               datasetPerBout1(i).DLRoll = 1;                           
               datasetPerBout1(i).LongRollOnset = NaN;
               datasetPerBout1(i).LongRoll = 0;
               datasetPerBout1(i).RolloverTime = NaN;
    
                      else
               datasetPerBout1(i).DLRoll = 2;                           
               datasetPerBout1(i).LongRollOnset = RolloverSequence (1) *2.5;
               datasetPerBout1(i).LongRoll = 1;
               datasetPerBout1(i).RolloverTime = RolloverSequence;
                %%DL rollover onset and stop
               RollGap = find(diff(RolloverSequence)>1);
               NbRollover   = length(RollGap);
               datasetPerBout1(i).RolloverOnset = [datasetPerBout1(i).RolloverTime(1), datasetPerBout1(i).RolloverTime(RollGap+1)]; 
               datasetPerBout1(i).RolloverStop = [datasetPerBout1(i).RolloverTime(RollGap)+5,datasetPerBout1(i).RolloverTime(end)+5]; 
               datasetPerBout1(i).RolloverDuration = 2.5*( datasetPerBout1(i).RolloverStop - datasetPerBout1(i).RolloverOnset+1);
                      end
                  end    
     
             end
            if  isnan (datasetPerBout1(i). RolloverOnset)
                datasetPerBout1(i). RolloverNb = 0;
                datasetPerBout1(i). RolloverStart = NaN;
            else
           datasetPerBout1(i). RolloverNb = length(datasetPerBout1(i).RolloverOnset);
           datasetPerBout1(i). RolloverStart = 2.5*datasetPerBout1(i).RolloverOnset(1);
%           datasetPerBout1(i). RolloverStart = 1.58*datasetPerBout1(i).RolloverOnset(1);
          
            end
             
           
           if (datasetPerBout1(i).Condition) < 37
               datasetPerBout1(i).Trial = 1;
           end

            if ((datasetPerBout1(i).Condition) > 36  & (datasetPerBout1(i).Condition) <49)
                 datasetPerBout1(i).Trial = 2;
            end
            if ((datasetPerBout1(i).Condition) > 48  & (datasetPerBout1(i).Condition) <67)
                 datasetPerBout1(i).Trial = 3;
            end            
            if (datasetPerBout1(i).Condition) > 66  
                 datasetPerBout1(i).Trial = 4;
            end            
%             
end



datasetPerBout1= datasetPerBout1([datasetPerBout1(:).New_BoutStart]>78); %response after stimulus
datasetPerBout1= datasetPerBout1([datasetPerBout1(:).BendTime]>0);% response after stimulus
datasetEscape= datasetPerBout1([datasetPerBout1(:).Cbend]>0);%all escapes here
datasetEscape= datasetEscape([datasetEscape(:).DLRoll]<3);
datasetEscape = datasetEscape(([datasetEscape(:).TotalDistance]>0) & ([datasetEscape(:).TotalDistance]<25)); %(30ms, 19+127 = 146)

%%
%to remove fish that only C-start escaped once or twice to 10 stimulus
for i = 1:84  
    
    if length(datasetEscape([datasetEscape.Condition]==i) )<3
    datasetEscape([datasetEscape.Condition]==i)=[];
    end
end
%%
datasetOutput = struct([]);
for i = 1 :length(datasetEscape)
     datasetOutput(i).CLUTCH = datasetEscape(i).Trial;
    datasetOutput(i).FISH = datasetEscape(i).Condition;
    datasetOutput(i).CONDITION = datasetEscape(i).Genotype;
    datasetOutput(i).TRIAL = datasetEscape(i).Nstim;
    datasetOutput(i).Latency = datasetEscape(i).Latency;
    datasetOutput(i).BoutDuration = datasetEscape(i).BoutDuration;
    datasetOutput(i).BoutDistance = datasetEscape(i).TotalDistance;
    datasetOutput(i).Speed = datasetEscape(i).Speed;
    datasetOutput(i).NbOsc = datasetEscape(i).NumberOfOscillations;
    datasetOutput(i).mTBF = datasetEscape(i).mTBF;    
    datasetOutput(i).Cbend = datasetEscape(i).Cbend;
    datasetOutput(i).Counterbend = datasetEscape(i).Counterbend;
    datasetOutput(i).BendTime = datasetEscape(i).BendTime;
    datasetOutput(i).BeatFreq = datasetEscape(i).BeatFreq;
    datasetOutput(i).RollIntegral = datasetEscape(i).RollIntegral ;
    datasetOutput(i).RollISum = datasetEscape(i).RollSum ;
    datasetOutput(i).RolloverNb = datasetEscape(i).RolloverNb ;
    datasetOutput(i).LongRollOnset = datasetEscape(i).RolloverStart ;
    datasetOutput(i).RollingOnset = datasetEscape(i).RollStart ;
end   
%   
% ParameterTable = struct2table(datasetOutput);
% ToCoef1 = ParameterTable(:, (5:21));
% ToCoef = table2array (ToCoef1);
% [R,P,RL,RU]  = corrcoef(ToCoef);
% figure;
% subplot(1,2,1)
% h1 = heatmap(R) ; 
% 
% subplot(1,2,2)
% h2 = heatmap(P) ;

 %%   
    
    
Fish_ID= unique([datasetEscape.Condition]);

% Control
idx_Control = find([datasetEscape.Genotype] == 0);
fish_Control = unique([datasetEscape(idx_Control).Condition]);

% rostral KA ablation
idx_KA = find([datasetEscape.Genotype] == 1);
fish_KA =  unique([datasetEscape(idx_KA).Condition]);

% caudal KA ablation
idx_cKA = find([datasetEscape.Genotype] == 2);
fish_cKA =  unique([datasetEscape(idx_cKA).Condition]);

%% Calculation for each Stimulus

for n=1:nbStimulus
    idx=[];
    datasetStim = datasetEscape([datasetEscape(:).Nstim]== n);% datasetStim: extract dataset for each Stimulus
 
    Fish= unique([datasetStim.Condition]);
 
    for i = 1:length(Fish)
       idx{Fish(i)}=find([datasetStim.Condition]== Fish(i)); 
           
          for l= 1:length(idx{Fish(i)})  
           EscapeRatio(n, Fish(i)) = 1;
           EscapeType (n , Fish(i)) = datasetStim(idx{Fish(i)}).EscapeType;
           FishTrial (n , Fish(i)) = datasetStim(idx{Fish(i)}).Trial;
           FishID(n , Fish(i)) = datasetStim(idx{Fish(i)}).Condition;
           FishGeno(n , Fish(i)) = datasetStim(idx{Fish(i)}).Genotype;
           Latency(n , Fish(i))= datasetStim(idx{Fish(i)}).Latency;
           Distance(n , Fish(i)) = datasetStim(idx{Fish(i)}).TotalDistance;
           BoutDuration (n , Fish(i))= datasetStim(idx{Fish(i)}).BoutDuration;
           NumberOfOscillations(n , Fish(i)) = datasetStim(idx{Fish(i)}).NumberOfOscillations;
           
           Cbend_Amplitude = abs([datasetStim(idx{Fish(i)}).Bend_Amplitude(1)]);    
           Cbend_Pos= find(abs([datasetStim(idx{Fish(i)}).Bend_Amplitude(:)])== Cbend_Amplitude);   
           Counterbend(n , Fish(i)) = 57.2958*(abs([datasetStim(idx{Fish(i)}(l)).Bend_Amplitude(Cbend_Pos+1)]));
           Cbend(n , Fish(i)) = datasetStim(idx{Fish(i)}).Cbend;
%            Roll(n , Fish(i)) = datasetStim(idx{Fish(i)}).DLRoll;
           RollIntegral(n , Fish(i)) = datasetStim(idx{Fish(i)}).RollIntegral;
           RollSum(n , Fish(i)) = datasetStim(idx{Fish(i)}).RollSum;
%            LongRoll(n , Fish(i)) = datasetStim(idx{Fish(i)}).LongRoll;
           RolloverOnset(n , Fish(i)) = datasetStim(idx{Fish(i)}).RolloverStart;
           RollingOnset(n , Fish(i)) = datasetStim(idx{Fish(i)}).RollStart;
           RolloverNb (n , Fish(i)) = datasetStim(idx{Fish(i)}).RolloverNb ;
           BendTime(n , Fish(i)) = datasetStim(idx{Fish(i)}).BendTime;
           BeatTime(n , Fish(i)) = datasetStim(idx{Fish(i)}).BeatTime;
           BeatFreq(n , Fish(i)) = datasetStim(idx{Fish(i)}).BeatFreq;

            mTBF  (n , Fish(i))=  datasetStim(idx{Fish(i)}).mTBF; 
           Speed(n , Fish(i)) = datasetStim(idx{Fish(i)}(l)).Speed;
          
          end

            
        end
        
end

%%
output = struct( 'Fish_ID', [], 'FishGeno',[],'FishTrial',[], ...
       'meanLatency',[],'meanCbend',[],'meanBendTime',[],'meanDistance',[],'meanBoutDuration',[],'meanSpeed',[],'meanNumberOfOscillations',[],...
       'meanTBFm',[],'meanTBF1',[],'meanCounterBend', []);
         
output.Fish_ID = (nanmean(FishID))';
output.FishGeno = (nanmean(FishGeno))';
output.FishTrial = (nanmean(FishTrial))';
output.EscapeRatio = (nanmean(EscapeRatio))';

output.CstartRatio =[nanmean(EscapeType)]';

output.meanLatency =  (nanmean(Latency))';
output.meanBendTime = (nanmean(BendTime))';
output.meanBeatFreq = (nanmean(BeatFreq))';
output.meanCbend =  (nanmean (Cbend))';
output.meanBoutDuration = (nanmean (BoutDuration))';
output.meanNumberOfOscillations = (nanmean ( NumberOfOscillations))'; % 
output.meanDistance  =  (nanmean (Distance))';
output.meanSpeed   =    (nanmean (Speed))';
output.meanTBFm =   (nanmean(mTBF))';
output.meanTBF1 =   (nanmean (TBF2))';
output.meanCounterBend = (nanmean (Counterbend))';
output.meanRollIntegral = (nanmean (RollIntegral))';
output.meanRollDuration = (nanmean (RollSum))';
output.meanRolloverNb = (nanmean (RolloverNb))';
output.meanRolloverOnset = (nanmean (RolloverOnset))';
output.meanRollingOnset = (nanmean (RollingOnset))';

T = struct2table (output);


exlfilename = 'Escape_fins_ablation.xlsx';
writetable(T,exlfilename,'Sheet',1);


%%
%plot escape parameters mean per fish
f1 = figure(1);

title ('All trials mean '); hold on;
subplot(3,4,1)
scatter_meanPerFish_ControlvsAblationRvsAblationC(Latency,fish_Control,fish_KA, fish_cKA);
%ylim([0.1 0.4]);
ylabel('Latency (ms)','FontSize',10);
hold off;

subplot(3,4,2)
scatter_meanPerFish_ControlvsAblationRvsAblationC(Cbend,fish_Control,fish_KA, fish_cKA);
%ylim([0.1 0.4]);
ylabel('Cbend (degree)','FontSize',10);
hold off;

subplot(3,4,3)
scatter_meanPerFish_ControlvsAblationRvsAblationC(BendTime,fish_Control,fish_KA, fish_cKA);
%ylim([0.1 0.4]);
ylabel('BendTime (ms)','FontSize',10);
hold off;

subplot(3,4,4)
scatter_meanPerFish_ControlvsAblationRvsAblationC(BeatFreq,fish_Control,fish_KA, fish_cKA);
%ylim([0.1 0.4]);
ylabel('TBF1(Hz)','FontSize',10);
hold off;

subplot(3,4,5)
scatter_meanPerFish_ControlvsAblationRvsAblationC(mTBF,fish_Control,fish_KA, fish_cKA);
%ylim([0.1 0.4]);
ylabel('mTBF(Hz)','FontSize',10);
hold off;

subplot(3,4,6)
scatter_meanPerFish_ControlvsAblationRvsAblationC(BoutDuration,fish_Control,fish_KA, fish_cKA);
%ylim([0.1 0.4]);
ylabel('BoutDurations (sec)','FontSize',10);
hold off;

subplot(3,4,7)
scatter_meanPerFish_ControlvsAblationRvsAblationC(Distance,fish_Control,fish_KA, fish_cKA);
% ylim([5 20]);
ylabel('Distance (mm)','FontSize',10);
% ytickangle(90);
hold off;
subplot(3,4,8)
scatter_meanPerFish_ControlvsAblationRvsAblationC(Speed,fish_Control,fish_KA, fish_cKA);
% ylim([5 20]);
ylabel('Speed','FontSize',10);
% ytickangle(90);
hold off;

subplot(3,4,9)
scatter_meanPerFish_ControlvsAblationRvsAblationC(NumberOfOscillations,fish_Control,fish_KA, fish_cKA);
% ylim([5 20]);
ylabel('NumberOfOscillations','FontSize',10);
hold off;

subplot(3,4,10 )
scatter_meanPerFish_ControlvsAblationRvsAblationC(EscapeRatio,fish_Control,fish_KA, fish_cKA);hold on;

ylabel('EscapeRatio','FontSize',10);
hold off;

subplot(3,4,11 )
scatter_meanPerFish_ControlvsAblationRvsAblationC(EscapeType,fish_Control,fish_KA, fish_cKA);hold on;

ylabel('Escape cstart','FontSize',10);
hold off;
% title ('Escape');
saveas(f1,['mean per fish 1.fig'])

%%
%plot rolling mean per fish
f2 = figure(2);
subplot(2,4,1)
scatter_meanPerFish_ControlvsAblationRvsAblationC(RolloverOnset,fish_Control,fish_KA, fish_cKA);

ylabel('RolloverOnset (ms)','FontSize',10);

hold off;

subplot(2,4,2)
scatter_meanPerFish_ControlvsAblationRvsAblationC(RollSum,fish_Control,fish_KA, fish_cKA);
ylabel('Rolling Duration','FontSize',10);
hold off;

subplot(2,4,3)
scatter_meanPerFish_ControlvsAblationRvsAblationC(RolloverNb ,fish_Control,fish_KA, fish_cKA);
ylabel('LongRollNb (n)','FontSize',10);
hold off;
subplot(2,4,4)
scatter_meanPerFish_ControlvsAblationRvsAblationC(RollIntegral ,fish_Control,fish_KA, fish_cKA);
ylabel('RollIntegral','FontSize',10);
hold off;
saveas(f2,['rolling mean per fish 1.fig'])


%%
% per trial (habituation)
f3= figure(3);
subplot(3,4,1)
plot_meanPerStimulus(Latency,fish_Control,fish_KA, fish_cKA);
ylabel(" mean Latency ");hold off;
subplot(3,4,2)
plot_meanPerStimulus(Cbend,fish_Control,fish_KA, fish_cKA);hold on;
ylabel(" mean Cbend ");hold on;
subplot(3,4,3)
plot_meanPerStimulus(BendTime,fish_Control,fish_KA, fish_cKA);
ylabel(" BendTime ");hold off;
subplot(3,4,4)
plot_meanPerStimulus(BeatFreq,fish_Control,fish_KA, fish_cKA);
ylabel(" mean TBF1 ");hold off;

subplot(3,4,5)
plot_meanPerStimulus(BoutDuration,fish_Control,fish_KA, fish_cKA);hold on;
ylabel("mean Duration (s) ");hold on;

subplot(3,4,6)
plot_meanPerStimulus(NumberOfOscillations,fish_Control,fish_KA, fish_cKA);hold on;
ylabel("mean Oscillations ");hold on;


subplot(3,4,7)
plot_meanPerStimulus(Distance,fish_Control,fish_KA, fish_cKA);hold on;
ylabel(" mean Distance (mm) ");hold on;

subplot(3,4,8)
plot_meanPerStimulus(Speed,fish_Control,fish_KA, fish_cKA);hold on;
ylabel(" mean Speed (mm/s) ");hold on;
subplot(3,4,9)
plot_meanPerStimulus(Counterbend,fish_Control,fish_KA, fish_cKA);hold on;
ylabel(" mean Counterbend (deg) ");hold on;
subplot(3,4,10)
plot_meanPerStimulus(RollSum,fish_Control,fish_KA, fish_cKA);
ylabel(" mean RollSum ");hold off;

subplot(3,4,11)
plot_meanPerStimulus(RollIntegral,fish_Control,fish_KA, fish_cKA);
ylabel(" mean RollIntegral ");hold off;

subplot(3,4,12)
plot_meanPerStimulus(RolloverNb,fish_Control,fish_KA, fish_cKA);
ylabel(" mean RolloverNb ");hold off;

%%

datasetDLRoll = datasetEscape([datasetEscape.DLRoll] ==2);
f4= figure(4); 
RolloverStartCT = [datasetDLRoll([datasetDLRoll.Genotype]==0).RolloverStart];
RolloverStartR = [datasetDLRoll([datasetDLRoll.Genotype]==1).RolloverStart];
RolloverStartC = [datasetDLRoll([datasetDLRoll.Genotype]==2).RolloverStart];
xaxisRollCT = ones(1,length(RolloverStartCT));
xaxisRollR = ones(1,length(RolloverStartR))*2;
xaxisRollC = ones(1,length(RolloverStartC))*3; 
group = [xaxisRollCT'; xaxisRollR'; xaxisRollC'];

boxplot([RolloverStartCT';RolloverStartR'; RolloverStartC'],group );hold on;
plot (1, mean(RolloverStartCT), 'bp');hold on;
plot (2, mean(RolloverStartR), 'bp');hold on;
plot (3, mean(RolloverStartC), 'bp');hold on;
grid off;

xlim([0.5 3.5]);
xticks(1:3);    
xticklabels({'Control','AblationR', 'AblationC'});

