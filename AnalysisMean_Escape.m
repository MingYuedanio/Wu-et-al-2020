%% Clean workspace and load data
close all
clear all
clc

load('KA_PTU_ablation_5Clutches_AllData_20200301.mat')

%% Initiate matrix

nbStimulus = 10;
nbIdx =104;


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
NewDuration = zeros(nbStimulus,nbIdx)*NaN;
NumberOfOscillations= zeros(nbStimulus,nbIdx)*NaN;
Distance= zeros(nbStimulus,nbIdx)*NaN;
Speed= zeros(nbStimulus,nbIdx)*NaN;
mSpeed= zeros(nbStimulus,nbIdx)*NaN;
iSpeed= zeros(nbStimulus,nbIdx)*NaN;
mNewSpeed= zeros(nbStimulus,nbIdx)*NaN;
mTBF= zeros(nbStimulus,nbIdx)*NaN;
mNewTBF= zeros(nbStimulus,nbIdx)*NaN;
miTBF= zeros(nbStimulus,nbIdx)*NaN;
TBF2 = zeros (nbStimulus, nbIdx)*NaN;

Counterbend= zeros(nbStimulus,nbIdx)*NaN;
Thirdbend= zeros(nbStimulus,nbIdx)*NaN;
ThirdbendvsCbend= zeros(nbStimulus,nbIdx)*NaN;
RollIntegral  =  zeros(nbStimulus,nbIdx)*NaN;
RollSum=  zeros(nbStimulus,nbIdx)*NaN;
LongRoll = zeros(nbStimulus,nbIdx)*NaN;
RolloverOnset = zeros(nbStimulus,nbIdx)*NaN;
RolloverNb = zeros(nbStimulus,nbIdx)*NaN;

EscapeType = zeros(nbStimulus,nbIdx)*NaN;
EscapeRatio = zeros(nbStimulus,nbIdx);

%% Calculation for each Stimulus
datasetPerBout2 = datasetPerBout([datasetPerBout(:).ManualRoll]<2);
datasetPerBout3 = datasetPerBout2([datasetPerBout2(:).NumberOfOscillations]>1.5);
 datasetPerBout1 = datasetPerBout3([datasetPerBout3(:).NStim]<11);
datasetPerBout1 = datasetPerBout1(([datasetPerBout1(:).BoutStart]>120) & ([datasetPerBout1(:).BoutStart]<158)); %(50ms)
datasetPerBout1= datasetPerBout1([datasetPerBout1(:).Condition]~=75);
RolloverSeq = [1,1,1,1,1,1];%%10 ms
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
            datasetPerBout1(i). New_BoutStart = diff_Angle(1)+ datasetPerBout1(i).BoutStart+1;
            datasetPerBout1(i). Latency = ((diff_Angle(1)+ datasetPerBout1(i).BoutStart)-126)/634*1000;
            datasetPerBout1(i).NewDuration=length([diff_Angle(1):1:New_MovePosition(end)])/634;
            datasetPerBout1(i).mSpeed=datasetPerBout1(i).TotalDistance/datasetPerBout1(i).BoutDuration;
            datasetPerBout1(i).iSpeed=nanmedian (datasetPerBout1(i).InstantaneousSpeed);
            datasetPerBout1(i).mNewSpeed=datasetPerBout1(i).TotalDistance/datasetPerBout1(i).NewDuration;
            datasetPerBout1(i).Cbend = abs(datasetPerBout1(i).AmplitudeOfAllBends(1));
            datasetPerBout1(i).BendTime = (datasetPerBout1(i).Bend_Timing(1)-diff_Angle(1))*1.58;
            datasetPerBout1(i).BeatTime = (datasetPerBout1(i).Bend_Timing(2)-datasetPerBout1(i).Bend_Timing(1))*1.58;
            datasetPerBout1(i).BeatFreq = 1000/((datasetPerBout1(i).Bend_Timing(2)-datasetPerBout1(i).Bend_Timing(1))*1.58*2);%TBF1
            datasetPerBout1(i).Counterbend = abs(datasetPerBout1(i).AmplitudeOfAllBends(2));
            
            if datasetPerBout1(i).Cbend <60
                 datasetPerBout1(i).EscapeType  = 0; %S escape
             else
                 datasetPerBout1(i).EscapeType = 1; % C escape
             end
            datasetPerBout1(i).Thirdbend = abs(datasetPerBout1(i).AmplitudeOfAllBends(3));
            datasetPerBout1(i).bendRatio = (datasetPerBout1(i).Counterbend)/(datasetPerBout1(i).Cbend) ;
            datasetPerBout1(i).CCbendRatio = (datasetPerBout1(i).Cbend)/ (datasetPerBout1(i).Counterbend);
            datasetPerBout1(i).RollIntegral = mean (datasetPerBout1(i).rolloverProbability);  
           % datasetPerBout1(i).RollSum = sum (datasetPerBout1(i).rolloverProbability([datasetPerBout1(i).rolloverProbability]>0.8))/length ((datasetPerBout1(i).rolloverProbability([1:end])>0.8));  
            datasetPerBout1(i).RollSum =  1.58 * sum ([datasetPerBout1(i).rolloverProbability]>0.80);  
           
            datasetPerBout1(i).TBF2 = datasetPerBout1(i).InstantaneousTBF(2);  
            datasetPerBout1(i).mTBF = (datasetPerBout1(i).NumberOfOscillations)/ (datasetPerBout1(i).BoutDuration);  
            datasetPerBout1(i).mNewTBF = 3000/ (4*(datasetPerBout1(i).BendTime + datasetPerBout1(i).BeatTime));
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
               datasetPerBout1(i).RollStart = 1.58* RollingSequence(1);         
                      if isempty (RolloverSequence)
               datasetPerBout1(i).DLRoll = 1;                           
               datasetPerBout1(i).LongRollOnset = NaN;
               datasetPerBout1(i).LongRoll = 0;
               datasetPerBout1(i).RolloverTime = NaN;
    
                      else
               datasetPerBout1(i).DLRoll = 2;                           
               datasetPerBout1(i).LongRollOnset = RolloverSequence (1) *1.58;
               datasetPerBout1(i).LongRoll = 1;
               datasetPerBout1(i).RolloverTime = RolloverSequence;
                %%DL rollover onset and stop
               RollGap = find(diff(RolloverSequence)>1);
               NbRollover   = length(RollGap);
               datasetPerBout1(i).RolloverOnset = [datasetPerBout1(i).RolloverTime(1), datasetPerBout1(i).RolloverTime(RollGap+1)]; 
               datasetPerBout1(i).RolloverStop = [datasetPerBout1(i).RolloverTime(RollGap)+5,datasetPerBout1(i).RolloverTime(end)+5]; 
               datasetPerBout1(i).RolloverDuration = 1.58*( datasetPerBout1(i).RolloverStop - datasetPerBout1(i).RolloverOnset+1);
                      end
                  end    
     
             end
            if  isnan (datasetPerBout1(i). RolloverOnset)
                datasetPerBout1(i). RolloverNb = 0;
                datasetPerBout1(i). RolloverStart = NaN;
            else
           datasetPerBout1(i). RolloverNb = length(datasetPerBout1(i).RolloverOnset);
           datasetPerBout1(i). RolloverStart = 1.58*datasetPerBout1(i).RolloverOnset(1);
            end
             
           
           if (datasetPerBout1(i).Condition) < 17
               datasetPerBout1(i).Trial = 1;
           end

            if ((datasetPerBout1(i).Condition) > 16  & (datasetPerBout1(i).Condition) <33)
                 datasetPerBout1(i).Trial = 2;
            end
            if ((datasetPerBout1(i).Condition) > 32  & (datasetPerBout1(i).Condition) <53)
                 datasetPerBout1(i).Trial = 3;
            end            
            if ((datasetPerBout1(i).Condition) > 52  & (datasetPerBout1(i).Condition) <75)
                 datasetPerBout1(i).Trial = 4;
            end            
%             
            if (datasetPerBout1(i).Condition) > 74
                   datasetPerBout1(i).Trial = 5;
            end
end



datasetPerBout1= datasetPerBout1([datasetPerBout1(:).New_BoutStart]>125);
datasetPerBout1= datasetPerBout1([datasetPerBout1(:).BendTime]>0);
datasetEscape= datasetPerBout1([datasetPerBout1(:).Cbend]>0);
datasetEscape= datasetEscape([datasetEscape(:).DLRoll]<3);

%%
%to remove fish that only C-start escaped once or twice to 10 stimulus
for i = 1:104   
    
    if length(datasetEscape([datasetEscape.Condition]==i) )<1
    datasetEscape([datasetEscape.Condition]==i)=[];
    end
end
%%
datasetOutput = struct([]);
for i = 1 :length(datasetEscape)
    datasetOutput(i).CLUTCH = datasetEscape(i).Trial;
    datasetOutput(i).FISH = datasetEscape(i).Condition;
    datasetOutput(i).CONDITION = datasetEscape(i).Genotype;
    datasetOutput(i).TRIAL = datasetEscape(i).NStim;
    datasetOutput(i).Latency = datasetEscape(i).Latency;
    datasetOutput(i).BoutDuration = datasetEscape(i).BoutDuration;
    datasetOutput(i).BoutDistance = datasetEscape(i).TotalDistance;
    datasetOutput(i).NbOsc = datasetEscape(i).NumberOfOscillations;
    datasetOutput(i).mTBF = datasetEscape(i).mTBF;    
    datasetOutput(i).mNewTBF = datasetEscape(i).mNewTBF; 
    datasetOutput(i).Speed = datasetEscape(i).mSpeed;
    datasetOutput(i).Cbend = datasetEscape(i).Cbend;
    datasetOutput(i).BendTime = datasetEscape(i).BendTime;
    datasetOutput(i).BeatFreq = datasetEscape(i).BeatFreq;
    datasetOutput(i).Counterbend = datasetEscape(i).Counterbend;
    datasetOutput(i).Thirdbend = datasetEscape(i).Thirdbend;

    datasetOutput(i).bendRatio = datasetEscape(i).bendRatio;
    datasetOutput(i).RollIntegral = datasetEscape(i).RollIntegral ;
    datasetOutput(i).RollISum = datasetEscape(i).RollSum ;
    datasetOutput(i).RolloverNb = datasetEscape(i).RolloverNb ;
    datasetOutput(i).RolloverOnset = datasetEscape(i).RolloverStart ;
end   
%   
ParameterTable = struct2table(datasetOutput);
ToCoef1 = ParameterTable(:, (5:21));
ToCoef = table2array (ToCoef1);
[R,P,RL,RU]  = corrcoef(ToCoef);
figure;
subplot(1,2,1)
h1 = heatmap(R) ; 

subplot(1,2,2)
h2 = heatmap(P) ;

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
    datasetStim = datasetEscape([datasetEscape(:).NStim]== n);% datasetStim: extract dataset for each Stimulus
 
    Fish= unique([datasetStim.Condition]);
 
    for i = 1:length(Fish)
       idx{Fish(i)}=find([datasetStim.Condition]== Fish(i)); 
           
          for l= 1:length(idx{Fish(i)})  
           EscapeRatio(n, Fish(i)) = 1;
           EscapeType (n , Fish(i)) = datasetStim(idx{Fish(i)}).EscapeType;
           FishTrial (n , Fish(i)) = datasetStim(idx{Fish(i)}).Trial;
           FishID(n , Fish(i)) = datasetStim(idx{Fish(i)}).Condition;
           FishGeno(n , Fish(i)) = datasetStim(idx{Fish(i)}).Genotype;
           Latency(n , Fish(i))= ((datasetStim(idx{Fish(i)}).New_BoutStart- 126)/635)* 1000;
           Distance(n , Fish(i)) = datasetStim(idx{Fish(i)}).TotalDistance;
           BoutDuration (n , Fish(i))= datasetStim(idx{Fish(i)}).BoutDuration;
           NewDuration (n , Fish(i))= datasetStim(idx{Fish(i)}).NewDuration;
           NumberOfOscillations(n , Fish(i)) = datasetStim(idx{Fish(i)}).NumberOfOscillations;
           
           Cbend_Amplitude = abs([datasetStim(idx{Fish(i)}).Bend_Amplitude(1)]);    
           Cbend_Pos= find(abs([datasetStim(idx{Fish(i)}).Bend_Amplitude(:)])== Cbend_Amplitude);   
           Counterbend(n , Fish(i)) = 57.2958*(abs([datasetStim(idx{Fish(i)}(l)).Bend_Amplitude(Cbend_Pos+1)]));
           Cbend(n , Fish(i)) = datasetStim(idx{Fish(i)}).Cbend;
           Roll(n , Fish(i)) = datasetStim(idx{Fish(i)}).DLRoll;
           RollIntegral(n , Fish(i)) = datasetStim(idx{Fish(i)}).RollIntegral;
           RollSum(n , Fish(i)) = datasetStim(idx{Fish(i)}).RollSum;
           LongRoll(n , Fish(i)) = datasetStim(idx{Fish(i)}).LongRoll;
           RolloverOnset(n , Fish(i)) = datasetStim(idx{Fish(i)}).RolloverStart;
           RolloverNb (n , Fish(i)) = datasetStim(idx{Fish(i)}).RolloverNb ;
           BendTime(n , Fish(i)) = datasetStim(idx{Fish(i)}).BendTime;
           BeatTime(n , Fish(i)) = datasetStim(idx{Fish(i)}).BeatTime;
           BeatFreq(n , Fish(i)) = datasetStim(idx{Fish(i)}).BeatFreq;

           Thirdbend(n , Fish(i)) = 57.2958*(abs([datasetStim(idx{Fish(i)}(l)).Bend_Amplitude(Cbend_Pos+2)]));
           ThirdbendvsCbend(n , Fish(i)) = abs(datasetStim(idx{Fish(i)}(l)).Bend_Amplitude(Cbend_Pos+2)/datasetStim(idx{Fish(i)}(l)).Bend_Amplitude(Cbend_Pos));
          
           miTBF(n , Fish(i))= median([datasetStim(idx{Fish(i)}(l)).InstantaneousTBF]);
           mTBF  (n , Fish(i))=  datasetStim(idx{Fish(i)}).mTBF;
           mNewTBF  (n , Fish(i))=  datasetStim(idx{Fish(i)}).mNewTBF;
           TBF2  (n , Fish(i))=  datasetStim(idx{Fish(i)}).TBF2; 
           Speed(n , Fish(i)) = datasetStim(idx{Fish(i)}(l)).Speed;
           iSpeed(n , Fish(i)) = datasetStim(idx{Fish(i)}(l)).iSpeed;
           mSpeed(n , Fish(i)) = (datasetStim(idx{Fish(i)}(l)).TotalDistance)/(datasetStim(idx{Fish(i)}(l)).BoutDuration);
           mNewSpeed(n , Fish(i)) = (datasetStim(idx{Fish(i)}(l)).TotalDistance)/(datasetStim(idx{Fish(i)}(l)).NewDuration);
           
          end

            
        end
        
end

%%
output = struct( 'Fish_ID', [], 'FishGeno',[],'FishTrial',[], ...
       'meanLatency',[],'meanCbend',[],'meanBendTime',[],'meanDistance',[],'meanBoutDuration',[],'meanSpeed',[],'meanNumberOfOscillations',[],...
       'meanTBFm',[],'meanTBF1',[],'meanCounterBend', [],'meanThirdBend', [])
         
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
output.meanNewTBF =   (nanmedian(mNewTBF))'; % 
output.meanTBF1 =   (nanmean (TBF2))';
output.meanCounterBend = (nanmean (Counterbend))';
output.meanThirdBend = (nanmean (Thirdbend))'; 
output.meanThirdbendvsCbend = (nanmean (ThirdbendvsCbend))';

output.meanRollIntegral = (nanmean (RollIntegral))';
output.meanRollDuration = (nanmean (RollSum))';
output.meanRolloverNb = (nanmean (RolloverNb))';
output.meanRolloverOnset = (nanmean (RolloverOnset))';



%%
T = struct2table (output);
outputtable = rmmissing(T);

% exlfilename = 'AllTrialsOutput20200312_escape.xlsx';
% writetable(outputtable,exlfilename,'Sheet',1);


%%
%plot normal escape mean per fish
f4 = figure(4)

title ('All trials mean 1'); hold on;
subplot(2,4,1)
scatter_meanPerFish_ControlvsAblationRvsAblationC(BoutDuration,fish_Control,fish_KA, fish_cKA);
%ylim([0.1 0.4]);
ylabel('BoutDurations (sec)','FontSize',10);
hold off;

subplot(2,4,2)
scatter_meanPerFish_ControlvsAblationRvsAblationC(Distance,fish_Control,fish_KA, fish_cKA);
% ylim([5 20]);
ylabel('Distance (mm)','FontSize',10);
% ytickangle(90);
hold off;
subplot(2,4,3)
scatter_meanPerFish_ControlvsAblationRvsAblationC(Speed,fish_Control,fish_KA, fish_cKA);
% ylim([5 20]);
ylabel('Speed','FontSize',10);
% ytickangle(90);
hold off;
subplot(2,4,4)
scatter_meanPerFish_ControlvsAblationRvsAblationC(iSpeed,fish_Control,fish_KA, fish_cKA);
% ylim([5 20]);
ylabel('iSpeed','FontSize',10);
hold off;
subplot(2,4,5)
scatter_meanPerFish_ControlvsAblationRvsAblationC(iSpeed,fish_Control,fish_KA, fish_cKA);

ylabel('iSpeed','FontSize',10);

hold off;

subplot(2,4,6)
scatter_meanPerFish_ControlvsAblationRvsAblationC(Cbend,fish_Control,fish_KA, fish_cKA);

ylabel('Cbend','FontSize',10);

hold off;

subplot(2,4, 7)
scatter_meanPerFish_ControlvsAblationRvsAblationC(BendTime,fish_Control,fish_KA, fish_cKA);

ylabel('TimeToPeak (msec)','FontSize',10);

hold off;

subplot(2,4,8)
scatter_meanPerFish_ControlvsAblationRvsAblationC(BeatFreq,fish_Control,fish_KA, fish_cKA);

ylabel('BeatFreq (Hz)','FontSize',10);

hold off;


saveas(f4,['mean per fish 1.fig'])
%%
%plot rolling mean per fish
f5 = figure(5);
subplot(2,4,1)
scatter_meanPerFish_ControlvsAblationRvsAblationC(RolloverNb ,fish_Control,fish_KA, fish_cKA);
ylabel('RolloverNb (n)','FontSize',10);
hold off;

subplot(2,4,2)
scatter_meanPerFish_ControlvsAblationRvsAblationC(RolloverOnset,fish_Control,fish_KA, fish_cKA);
ylabel('Rollover Onset ','FontSize',10);
hold off;

subplot(2,4, 3)
scatter_meanPerFish_ControlvsAblationRvsAblationC(Roll,fish_Control,fish_KA, fish_cKA);
ylabel('DLRoll','FontSize',10);
hold off;


subplot(2,4,4)
scatter_meanPerFish_ControlvsAblationRvsAblationC(RollSum,fish_Control,fish_KA, fish_cKA);
ylabel('RollSum','FontSize',10);
hold off;

subplot(2,4, 5)
scatter_meanPerFish_ControlvsAblationRvsAblationC(RollIntegral,fish_Control,fish_KA, fish_cKA);
ylabel('RollIntegral','FontSize',10);
hold off;
saveas(f5,['rollover mean per fish 1.fig'])
%%
f6 = figure(6)

subplot(2,4, 1)
scatter_meanPerFish_ControlvsAblationRvsAblationC(mNewTBF,fish_Control,fish_KA, fish_cKA);hold on;
% ylim([30 70]);
ylabel('mean NewTBF(Hz)','FontSize',10);
hold off;


subplot(2,4, 2)
scatter_meanPerFish_ControlvsAblationRvsAblationC(mTBF,fish_Control,fish_KA, fish_cKA);hold on;
% ylim([30 70]);
ylabel('mean mTBF(Hz)','FontSize',10);
hold off;

subplot(2,4, 3)
scatter_meanPerFish_ControlvsAblationRvsAblationC(TBF2,fish_Control,fish_KA, fish_cKA);hold on;
% ylim([30 70]);
ylabel('mean TBF2(Hz)','FontSize',10);
hold off;
subplot(2,4,4 )

scatter_meanPerFish_ControlvsAblationRvsAblationC(Counterbend,fish_Control,fish_KA, fish_cKA);hold on;
% ylim([30 120]);

ylabel('Counter bend (degree)','FontSize',10);
hold off;

% 
% 
subplot(2,4,5 )

scatter_meanPerFish_ControlvsAblationRvsAblationC(Thirdbend,fish_Control,fish_KA, fish_cKA);hold on;
% ylim([30 120]);

ylabel('Thirdbend (degree)','FontSize',10);
hold off;

subplot(2,4,6 )
scatter_meanPerFish_ControlvsAblationRvsAblationC(ThirdbendvsCbend,fish_Control,fish_KA, fish_cKA);hold on;

ylabel('3rd bend vs C-bend','FontSize',10);
hold off;

subplot(2,4,7 )
scatter_meanPerFish_ControlvsAblationRvsAblationC(EscapeRatio,fish_Control,fish_KA, fish_cKA);hold on;

ylabel('EscapeRatio','FontSize',10);
hold off;
subplot(2,4,8 )
scatter_meanPerFish_ControlvsAblationRvsAblationC(EscapeType,fish_Control,fish_KA, fish_cKA);hold on;

ylabel('Escape cstart','FontSize',10);
hold off;
% title ('Escape');
saveas(f5,['mean per fish 2.fig'])

%%
% per trial (habituation)
f7 = figure(7);
subplot(3,4,1)
plot_meanPerStimulus(BoutDuration,fish_Control,fish_KA, fish_cKA);hold on;
ylabel("mean Duration (s) ");hold on;

subplot(3,4,2)
plot_meanPerStimulus(NumberOfOscillations,fish_Control,fish_KA, fish_cKA);hold on;
ylabel("mean Oscillations ");hold on;


subplot(3,4,3)
plot_meanPerStimulus(Distance,fish_Control,fish_KA, fish_cKA);hold on;
ylabel(" mean Distance (mm) ");hold on;

subplot(3,4,4)
plot_meanPerStimulus(Speed,fish_Control,fish_KA, fish_cKA);hold on;
ylabel(" mean Speed (mm/s) ");hold on;
 
 
subplot(3,4,5)
plot_meanPerStimulus(Latency,fish_Control,fish_KA, fish_cKA);
ylabel(" mean Latency ");hold off;

subplot(3,4,6)
plot_meanPerStimulus(Cbend,fish_Control,fish_KA, fish_cKA);hold on;
ylabel(" mean Cbend ");hold on;

subplot(3,4,7)
plot_meanPerStimulus(TBF2,fish_Control,fish_KA, fish_cKA);hold on;

ylabel(" mean TBF2 "); hold on;

subplot(3,4,8)
plot_meanPerStimulus(BeatFreq,fish_Control,fish_KA, fish_cKA);
ylabel(" mean beatfreq ");hold off;


subplot(3,4,9)
plot_meanPerStimulus(BendTime,fish_Control,fish_KA, fish_cKA);
ylabel(" BendTime ");hold off;


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
%trial

% f7= figure(7)
% subplot(4,4,1)
% plot_meanPerStimulusAll(Latency);
% ylabel(" mean Latency ");hold off;
% 
% subplot(4,4,2)
% plot_meanPerStimulusAll(Cbend);hold on;
% ylabel(" mean Cbend ");hold on;
% 
% subplot(4,4,3)
% plot_meanPerStimulusAll(BendTime);
% ylabel(" mean BendTime ");hold off;
% 
% subplot(4,4,4)
% plot_meanPerStimulusAll(Speed);hold on;
% ylabel(" mean Speed (mm/s) ");hold on;
% subplot(4,4,5)
% plot_meanPerStimulus(Latency,fish_Control,fish_KA, fish_cKA);
% ylabel(" mean Latency ");hold off;
% subplot(4,4,6)
% plot_meanPerStimulus(Cbend,fish_Control,fish_KA, fish_cKA);hold on;
% ylabel(" mean Cbend ");hold on;
% subplot(4,4,7)
% plot_meanPerStimulus(BendTime,fish_Control,fish_KA, fish_cKA);
% ylabel(" mean BendTime ");hold off;
% % subplot(4,4,7)
% % plot_meanPerStimulus(RollIntegral,fish_Control,fish_KA, fish_cKA);
% % ylabel(" mean RollDeepLearn ");hold off;
% % subplot(2,7,3)
% % plot_meanPerStimulus(Counterbend,fish_Control,fish_KA, fish_cKA);hold on;
% % ylabel(" mean Counter bend "); hold on;
% subplot(4,4,8)
% plot_meanPerStimulus(Speed,fish_Control,fish_KA, fish_cKA);hold on;
% ylabel(" mean Speed (mm/s) ");hold on;
% % subplot(4,4,8)
% % plot_meanPerStimulus(Roll,fish_Control,fish_KA, fish_cKA);hold on;
% % ylabel(" mean Roll (mm/s) ");hold on;
% subplot(4,4,9)
% plot_meanPerStimulusAll(NumberOfOscillations);hold on;
% ylabel("mean Oscillations ");hold on;
% 
% subplot(4,4,10)
% plot_meanPerStimulusAll(Distance);hold on;
% ylabel(" mean Distance (mm) ");hold on;
% 
% subplot(4,4,11)
% plot_meanPerStimulusAll(BoutDuration);hold on;
% ylabel("mean Duration (s) ");hold on;
% 
% subplot(4,4,12)
% plot_meanPerStimulusAll(BeatFreq);hold on;
% ylabel("mean  TBF1 (Hz) ");hold on;
% 
% 
% subplot(4,4,13)
% plot_meanPerStimulus(NumberOfOscillations,fish_Control,fish_KA, fish_cKA);hold on;
% ylabel("mean Oscillations ");hold on;
% 
% subplot(4,4,14)
% plot_meanPerStimulus(Distance,fish_Control,fish_KA, fish_cKA);hold on;
% ylabel(" mean Distance (mm) ");hold on;
% 
% subplot(4,4,15)
% plot_meanPerStimulus(BoutDuration,fish_Control,fish_KA, fish_cKA);hold on;
% ylabel("mean Duration (s) ");hold on;
% 
% subplot(4,4,16)
% plot_meanPerStimulus(BeatFreq,fish_Control,fish_KA, fish_cKA);hold on;
% ylabel("mean  TBF1 (Hz) ");hold on;
% 


f8 = figure(8)

subplot(2,4,1)
scatter2_meanPerFish_ControlvsAblationRvsAblationC(Cbend,fish_Control,fish_KA, fish_cKA);hold on;
% ylim([30 70]);
ylabel('S Start ratio','FontSize',10);
hold off;
subplot(2,4,2)
scatter2_meanPerFish_ControlvsAblationRvsAblationC(Latency,fish_Control,fish_KA, fish_cKA);hold on;
% ylim([30 70]);
ylabel('C Start ratio','FontSize',10);
hold off;

subplot(2,4,3)
scatter2_meanPerFish_ControlvsAblationRvsAblationC(BendTime,fish_Control,fish_KA, fish_cKA);hold on;
% ylim([30 70]);
ylabel('EscapeRatio','FontSize',10);
hold off;

subplot(2,4,4)
scatter2_meanPerFish_ControlvsAblationRvsAblationC(EscapeType,fish_Control,fish_KA, fish_cKA);hold on;
ylim([0 15]);
ylabel('Latency (msec)','FontSize',10);
hold off;

subplot(2,4,5)
scatter2_meanPerFish_ControlvsAblationRvsAblationC(Speed,fish_Control,fish_KA, fish_cKA);hold on;
% ylim([30 70]);
ylabel('mean bend rise','FontSize',10);
hold off;
% subplot(2,4,6)
% scatter2_meanPerFish_ControlvsAblationRvsAblationC(BendTime,fish_Control,fish_KA, fish_cKA);hold on;
% % ylim([30 70]);
% ylabel('median bend rise','FontSize',10);
% hold off;

% 
% %%
% 
% f10 = figure(10)
% subplot(3,4,1)
% plot_Distribution_all(BoutDuration,fish_Control,fish_KA, fish_cKA);hold on;
% xlabel("Duration ");hold on;
% % ax1 = gca;
% % ax1.XTick = [...];
% %     
% subplot(3,4,2)
% plot_Distribution_all(NumberOfOscillations,fish_Control,fish_KA, fish_cKA);hold on;
% xlabel("Oscillations ");hold on;
% 
% 
% subplot(3,4,3)
% plot_Distribution_all(Distance,fish_Control,fish_KA, fish_cKA);hold on;
% xlabel(" Distance (mm) ");hold on;
% 
% subplot(3,4,4)
% plot_Distribution_all(Speed,fish_Control,fish_KA, fish_cKA);hold on;
% xlabel(" Speed  ");hold on;
%  
%  
% subplot(3,4,5)
% plot_Distribution_all(mTBF,fish_Control,fish_KA, fish_cKA);hold on;
% 
% xlabel("mTBF  ");hold on;
% % ax5 = gca;
% % ax5.XTick = [...];
% 
% subplot(3,4,6)
% plot_Distribution_all(Cbend,fish_Control,fish_KA, fish_cKA);hold on;
% xlabel(" Cbend ");hold on;
% 
% subplot(3,4,7)
% plot_Distribution_all(Counterbend,fish_Control,fish_KA, fish_cKA);hold on;
% xlabel(" Counter bend "); hold on;
% 
% subplot(3,4,8)
% plot_Distribution_all(Thirdbend,fish_Control,fish_KA, fish_cKA);
% xlabel(" Thirdbend ");hold off;
% 
% subplot(3,4,9)
% plot_Distribution_all(Latency,fish_Control,fish_KA, fish_cKA);
% xlabel(" Latency ");hold off;
% 
% subplot(3,4,10)
% plot_Distribution_all(BeatTime,fish_Control,fish_KA, fish_cKA);
% xlabel(" BeatTime ");hold off;
% 
% subplot(3,4,11)
% plot_Distribution_all(BeatFreq,fish_Control,fish_KA, fish_cKA);
% xlabel(" BeatFreq ");hold off;
% 
% subplot(3,4,12)
% plot_Distribution_all(TBF2,fish_Control,fish_KA, fish_cKA);
% xlabel(" TBF2 ");hold off;


%%
datasetDLRoll  = datasetPerBout1([datasetPerBout1.DLRoll] >1);
RolloverDurationCT = [];
RolloverDurationR = [];
RolloverDurationC = [];
for i = 1:length(datasetDLRoll)

    RolloverDurationCurrent = datasetDLRoll(i).RolloverDuration(1:end);
     if    datasetDLRoll(i).Genotype ==0   
    RolloverDurationCT = [RolloverDurationCT RolloverDurationCurrent];
     elseif datasetDLRoll(i).Genotype == 1
             
    RolloverDurationR = [RolloverDurationR RolloverDurationCurrent];
         else
    RolloverDurationC = [RolloverDurationC RolloverDurationCurrent];
         end
         
   
end

figure;
Dist_G0 = histc(RolloverDurationCT, (0:3.16:100));
Density_G0= Dist_G0/length(RolloverDurationCT<100);
Dist_G1 = histc(RolloverDurationR, (0:3.16:100));
Density_G1= Dist_G1/length(RolloverDurationR<100);
Dist_G2 = histc(RolloverDurationC, (0:3.16:100));
Density_G2= Dist_G2/length(RolloverDurationC<100);

plot(Density_G0', 'k-o','MarkerSize', 2);hold on;
plot(Density_G1', 'r-o','MarkerSize', 2);hold on;
plot(Density_G2', 'g-o','MarkerSize', 2);
box off;













