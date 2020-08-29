
%% Clean workspace and load data
close all
clear 
clc
load ('CSFcN_ablation_Exploration.mat');

%% Initiate variables 

%  add trial info;    
for i = 1:length(datasetPerBout);
if (datasetPerBout(i).Condition) < 21
   datasetPerBout(i).NTrial = 5;
end

if ((datasetPerBout(i).Condition) > 20  & (datasetPerBout(i).Condition) <36)
     datasetPerBout(i).NTrial = 1;
end
if ((datasetPerBout(i).Condition) > 35  & (datasetPerBout(i).Condition) <53)
     datasetPerBout(i).NTrial = 2;
end
if ((datasetPerBout(i).Condition) >52  & (datasetPerBout(i).Condition) <75)
     datasetPerBout(i).NTrial = 4;
end
if (datasetPerBout(i).Condition) > 74
       datasetPerBout(i).NTrial = 3;
end

datasetPerBout(i).FirstBend  = max(abs(datasetPerBout(i).AmplitudeOfAllBends));
datasetPerBout(i).MedianBend  = median(abs(datasetPerBout(i).AmplitudeOfAllBends));
datasetPerBout(i).MeanBend  = mean(abs(datasetPerBout(i).AmplitudeOfAllBends));
datasetPerBout(i).MeanTBF  = (datasetPerBout(i).NumberOfOscillations)/(datasetPerBout(i).BoutDuration);
datasetPerBout(i).MaxTBF  = max(datasetPerBout(i).InstantaneousTBF);%[datasetBigBout(idx_BigBout{Fish_temp(i)}(h)).InstantaneousTBF]
datasetPerBout(i).HeadingDirection  = 57.2958*(abs(datasetPerBout(i).Heading(end)-datasetPerBout(i).Heading(1)));
% datasetPerBout(i).Deltaheading = abs(datasetPerBout(i).Polynomial_Deltaheading);

end

%%
%change the code here to select bout type to analyze
datasetPerBout1 = datasetPerBout;% all bouts
%  datasetPerBout1 = datasetPerBout([datasetPerBout(:).FirstBend]>30);%turn
%datasetPerBout1 = datasetPerBout([datasetPerBout(:).FirstBend]<30); %forward
%  datasetPerBout1 = datasetPerBout([datasetPerBout(:).HeadingDirection]>30);
%%
  
output = struct( 'NTrial', [], 'Fish_ID', [] )

window = [0:9600:48000];
Fish = unique([datasetPerBout1.Condition]);
for i = 1:length(Fish);
    i
datasetPerFish1(i) =  datasetPerFish([datasetPerFish(:).Condition]==Fish(i));%

end

for i=1:length(Fish);
    i
    if datasetPerFish1(i).TotalDistanceInActiveMotion <100   %%select swimmer >100mm
        datasetPerFish1(i).Swimmer = 0;
    else
        datasetPerFish1(i).Swimmer = 1;
    end
end
datasetPerFish2 =  datasetPerFish1([datasetPerFish1(:).Swimmer]==1);
Fish_temp=Fish ([datasetPerFish1.Swimmer] == 1); % Good Swimmer >100mm
%%


for i=1:length(Fish_temp);
    i
     % index for IBI
    index{Fish_temp(i)}= find(~([datasetPerBout1(:).Condition]-Fish_temp(i)));
    index{Fish_temp(i)}(1)=[];   
   
    Bouts= index{Fish_temp(i)};
    BoutsIBI_window1= index{Fish_temp(i)}(find([datasetPerBout1(index{Fish_temp(i)}).BoutStart]< window(2)));
    BoutsIBI_window2= index{Fish_temp(i)}(find([datasetPerBout1(index{Fish_temp(i)}).BoutStart]> window(2)&[datasetPerBout1(index{Fish_temp(i)}).BoutStart]< window(3)));
    BoutsIBI_window3= index{Fish_temp(i)}(find([datasetPerBout1(index{Fish_temp(i)}).BoutStart]> window(3)&[datasetPerBout1(index{Fish_temp(i)}).BoutStart]< window(4)));
    BoutsIBI_window4= index{Fish_temp(i)}(find([datasetPerBout1(index{Fish_temp(i)}).BoutStart]> window(4)&[datasetPerBout1(index{Fish_temp(i)}).BoutStart]< window(5)));
    BoutsIBI_window5= index{Fish_temp(i)}(find([datasetPerBout1(index{Fish_temp(i)}).BoutStart]> window(5)));
    
    
    
    % all index
    allindex{Fish_temp(i)} = find(~([datasetPerBout1(:).Condition]-Fish_temp(i)));


    Bouts_BigOscllations= allindex{Fish_temp(i)}(find([datasetPerBout1(allindex{Fish_temp(i)}).NumberOfOscillations]> 0));
    datasetBigBout =datasetPerBout1(Bouts_BigOscllations);    
    idx_BigBout{Fish_temp(i)}= find(~([datasetBigBout(:).Condition]-Fish_temp(i)));

        %calculate 
    IBI {Fish_temp(i)}  =     [datasetPerBout1(Bouts).InstantaneousIBI];
    IBI_window1 {Fish_temp(i)} = [datasetPerBout1(BoutsIBI_window1).InstantaneousIBI];
    IBI_window2{Fish_temp(i)}  = [datasetPerBout1(BoutsIBI_window2).InstantaneousIBI];
    IBI_window3{Fish_temp(i)}  = [datasetPerBout1(BoutsIBI_window3).InstantaneousIBI];
    IBI_window4 {Fish_temp(i)} = [datasetPerBout1(BoutsIBI_window4).InstantaneousIBI];
    IBI_window5 {Fish_temp(i)} = [datasetPerBout1(BoutsIBI_window5).InstantaneousIBI];
    
    NTrial{Fish_temp(i)} = [datasetPerBout1(Bouts_BigOscllations).NTrial];
    FishGenotype{Fish_temp(i)} = [datasetPerBout1(Bouts_BigOscllations).Genotype];
    TimeBout_pre{Fish_temp(i)} = [datasetPerBout1(Bouts_BigOscllations).BoutStart];
    
    BoutDuration{Fish_temp(i)}=[datasetPerBout1(Bouts_BigOscllations).BoutDuration];
    medianDuration{Fish_temp(i)}=nanmedian(BoutDuration{Fish_temp(i)});
    
    NumberOfOscillations{Fish_temp(i)}=[datasetPerBout1(Bouts_BigOscllations).NumberOfOscillations];
    medianNbOscillations {Fish_temp(i)}= nanmedian(NumberOfOscillations{Fish_temp(i)}); 
    meanNbOscillations {Fish_temp(i)}= nanmean(NumberOfOscillations{Fish_temp(i)}); 
        
    BoutDistance{Fish_temp(i)}=[datasetPerBout1(Bouts_BigOscllations).TotalDistance];
    medianBoutDistance{Fish_temp(i)}=nanmedian(BoutDistance{Fish_temp(i)});
    meanBoutDistance{Fish_temp(i)}=nanmean(BoutDistance{Fish_temp(i)});
    
    Speed{Fish_temp(i)}=[datasetPerBout1(Bouts_BigOscllations).Speed];
    medianSpeed{Fish_temp(i)}= nanmedian(Speed{Fish_temp(i)});
    meanSpeed{Fish_temp(i)}= nanmean(Speed{Fish_temp(i)});
    
    nBout{Fish_temp(i)}= length( Bouts);
    nBout_window1{Fish_temp(i)}= length( BoutsIBI_window1);
    nBout_window2{Fish_temp(i)}= length( BoutsIBI_window2);
    nBout_window3{Fish_temp(i)}= length( BoutsIBI_window3);
    nBout_window4{Fish_temp(i)}= length( BoutsIBI_window4);
    nBout_window5{Fish_temp(i)}= length( BoutsIBI_window5);
      
    
    if length(idx_BigBout{Fish_temp(i)})==0;
        TailAnglePerBout{Fish_temp(i)}=0;
        Bend_AmplitudePerBout{Fish_temp(i)}= 0;
        Pre_AmpMedianPerBout{Fish_temp(i)}=0;
        Pre_iAmpPerBout{Fish_temp(i)}=0;
        TBF{Fish_temp(i)}= 0;
        TBFmedianPerBout{Fish_temp(i)}= 0;
        
    else
        
       %calculate arraies of one bout.    
        for h=1:length(idx_BigBout{Fish_temp(i)});
            
            display(['currently processing fish ' num2str(i)])
            display(['currently processing bout number ' num2str(h)])
            
            % need to convert in degrees... otherwise : TailAngle{Fish(i)}{h}=[datasetPerBout(allindex{Fish(i)}(h)).TailAngle_smoothed]';
          
            TailAnglePerBout{Fish_temp(i)}{h}=57.2958*[datasetBigBout(idx_BigBout{Fish_temp(i)}(h)).TailAngle_smoothed]';
            Bend_TimingPerBout{Fish_temp(i)}{h}=[datasetBigBout(idx_BigBout{Fish_temp(i)}(h)).Bend_Timing]; 
            Bend_AmplitudePerBout{Fish_temp(i)}{h} = [datasetBigBout(idx_BigBout{Fish_temp(i)}(h)).AmplitudeOfAllBends];
            FirstAmpPerBout{Fish_temp(i)}(h)=abs(Bend_AmplitudePerBout{Fish_temp(i)}{h}(1));
            MaxAmpPerBout{Fish_temp(i)}(h)=max(abs(Bend_AmplitudePerBout{Fish_temp(i)}{h}));
            MeanAmpPerBout{Fish_temp(i)}(h)=mean(abs(Bend_AmplitudePerBout{Fish_temp(i)}{h}));
            MedianAmpPerBout{Fish_temp(i)}(h)=median(abs(Bend_AmplitudePerBout{Fish_temp(i)}{h}));
            RoutineTurnBigBout{Fish_temp(i)}(h) = (FirstAmpPerBout{Fish_temp(i)}(h)>30);
            
            
           
            TBF{Fish_temp(i)}{h} = [datasetBigBout(idx_BigBout{Fish_temp(i)}(h)).InstantaneousTBF];
            maxTBFPerBout{Fish_temp(i)}(h)   = nanmax( TBF{Fish_temp(i)}{h});
            meanTBFPerBout{Fish_temp(i)}(h)  = (datasetBigBout(idx_BigBout{Fish_temp(i)}(h)).NumberOfOscillations)/(datasetBigBout(idx_BigBout{Fish_temp(i)}(h)).BoutDuration);
            medianTBFPerBout{Fish_temp(i)}(h)= nanmedian(TBF{Fish_temp(i)}{h});
        end      
    end
    
    for m = 1:length(allindex{Fish_temp(i)});
            display(['currently processing fish ' num2str(i)])
            display(['currently processing bout number ' num2str(m)])
         HeadingPerBout{Fish_temp(i)}{m} = [datasetPerBout1(allindex{Fish_temp(i)}(m)).Heading];
         deltaHeadingPerBout{Fish_temp(i)}(m)= 57.2958* abs(HeadingPerBout{Fish_temp(i)}{m}(end)-HeadingPerBout{Fish_temp(i)}{m}(1));
         Abend{Fish_temp(i)}(m) = abs([datasetPerBout1(allindex{Fish_temp(i)}(m)).AmplitudeOfAllBends(1)]);
         RoutineTurnHeading {Fish_temp(i)}(m) = (deltaHeadingPerBout{Fish_temp(i)}(m) >70);
         RoutineTurnTail{Fish_temp(i)}(m) = (Abend{Fish_temp(i)}(m)>30);
    end   
        

    output(i).Fish_ID=Fish_temp(i);
    output(i).NTrial=unique(NTrial{Fish_temp(i)});
    output(i).FishGenotype=unique(FishGenotype{Fish_temp(i)}); 
    output(i).TotalDistance = datasetPerFish2(i).TotalDistanceActiveMotionAndGliding;
    output(i).mIBI=median(IBI{Fish_temp(i)},'omitnan');
    output(i).mIBI_window1=median(IBI_window1{Fish_temp(i)},'omitnan');
    output(i).mIBI_window2=median(IBI_window2{Fish_temp(i)},'omitnan');
    output(i).mIBI_window3=median(IBI_window3{Fish_temp(i)},'omitnan');
    output(i).mIBI_window4=median(IBI_window4{Fish_temp(i)},'omitnan');
    output(i).mIBI_window5=median(IBI_window5{Fish_temp(i)},'omitnan');
    
    output(i).meanBoutDuration_pre=mean(BoutDuration{Fish_temp(i)},'omitnan');
    output(i).meanNumberOfOscillations=mean(NumberOfOscillations{Fish_temp(i)},'omitnan');
    output(i).meanSpeed=mean(Speed{Fish_temp(i)},'omitnan');
    output(i).meanBoutDistance_pre=mean(BoutDistance{Fish_temp(i)},'omitnan');
 
    output(i).medianBoutDuration_pre=median(BoutDuration{Fish_temp(i)},'omitnan');
    output(i).medianNumberOfOscillations=median(NumberOfOscillations{Fish_temp(i)},'omitnan');
    output(i).medianSpeed=median(Speed{Fish_temp(i)},'omitnan');
    output(i).medianBoutDistance_pre=median(BoutDistance{Fish_temp(i)},'omitnan');
    
    output(i).BoutRate=nBout{Fish_temp(i)}/300;
    output(i).nBout=nBout{Fish_temp(i)};
    output(i).BoutRate_window1=nBout_window1{Fish_temp(i)}/60;
    output(i).BoutRate_window2=nBout_window2{Fish_temp(i)}/60;
    output(i).BoutRate_window3=nBout_window3{Fish_temp(i)}/60;
    output(i).BoutRate_window4=nBout_window4{Fish_temp(i)}/60;
    output(i).BoutRate_window5=nBout_window5{Fish_temp(i)}/60;     

    output(i).medianTBF=nanmedian(medianTBFPerBout{Fish_temp(i)});
    output(i).meanTBF=nanmedian (meanTBFPerBout{Fish_temp(i)});
    output(i).meanmaxTBF=nanmean (maxTBFPerBout{Fish_temp(i)});
    output(i).medianmaxTBF=nanmedian (maxTBFPerBout{Fish_temp(i)});    
    output(i).meanmedianTBF=mean(medianTBFPerBout{Fish_temp(i)},'omitnan');
    output(i).meanmeanTBF=mean (meanTBFPerBout{Fish_temp(i)});   
    output(i).deltaHeadingPerBout = nanmean (deltaHeadingPerBout{Fish_temp(i)});
    output(i).routineTurnRatioHead = sum (RoutineTurnHeading {Fish_temp(i)})/length( RoutineTurnHeading {Fish_temp(i)});
    output(i).forwardSwimRatioHead = 1- sum (RoutineTurnHeading {Fish_temp(i)})/length( RoutineTurnHeading {Fish_temp(i)});
    
    output(i).routineTurnRatioTail = sum (RoutineTurnTail {Fish_temp(i)})/length( RoutineTurnTail {Fish_temp(i)});
    output(i).forwardSwimRatioTail = 1- sum (RoutineTurnTail {Fish_temp(i)})/length( RoutineTurnTail {Fish_temp(i)});    
    
    output(i).medianmaxAmplitude_pre=median(MaxAmpPerBout{Fish_temp(i)},'omitnan');
    output(i).medianmeanAmplitude_pre=median(MeanAmpPerBout{Fish_temp(i)},'omitnan');
    output(i).medianmedianAmplitude_pre=median(MedianAmpPerBout{Fish_temp(i)},'omitnan');
    
    output(i).meanmaxAmplitude_pre=mean(MaxAmpPerBout{Fish_temp(i)},'omitnan');
    output(i).meanmeanAmplitude_pre=mean(MeanAmpPerBout{Fish_temp(i)},'omitnan');
    output(i).meanmedianAmplitude_pre=mean(MedianAmpPerBout{Fish_temp(i)},'omitnan');     

    i= i+1;

end;
%%
% 
% T = struct2table (output);
% exlfilename = 'AllTrialsSlowSwimOutput20190909.xlsx';
% writetable(T,exlfilename,'Sheet',2);
%% median Calculation
    
    G2medianBoutDuration_pre = [output([output.FishGenotype]==2). medianBoutDuration_pre];
    G2meanBoutDuration_pre   = [output([output.FishGenotype]==2). meanBoutDuration_pre];
 
    G2medianNumOfOsc=[output([output.FishGenotype]==2). medianNumberOfOscillations];
    G2meanNumOfOsc=[output([output.FishGenotype]==2). meanNumberOfOscillations];
 
    
    G2medianBoutDistance_pre=[output([output.FishGenotype]==2). medianBoutDistance_pre];
    G2meanBoutDistance_pre=[output([output.FishGenotype]==2). meanBoutDistance_pre];
    G2TotalDistance = [output([output.FishGenotype]==2).TotalDistance];
    
    G2medianSpeed_pre=[output([output.FishGenotype]==2). medianSpeed];
    G2meanSpeed_pre=[output([output.FishGenotype]==2). meanSpeed];

    % Calculate  TBF 
    G2medianTBF_pre =       [output([output.FishGenotype]==2). medianTBF];
    G2meanTBF_pre =         [output([output.FishGenotype]==2). meanTBF];
    G2meanmaxTBF =          [output([output.FishGenotype]==2). meanmaxTBF];
    G2medianmaxTBF =          [output([output.FishGenotype]==2). medianmaxTBF];
    G2meanmedianTBF_pre =   [output([output.FishGenotype]==2). meanmedianTBF];
    G2meanmeanTBF_pre =     [output([output.FishGenotype]==2). meanmeanTBF];
    G2medianIBIm_pre =      [output([output.FishGenotype]==2). mIBI];
    G2medianIBIm_window1 =  [output([output.FishGenotype]==2). mIBI_window1];
    G2medianIBIm_window2 =  [output([output.FishGenotype]==2). mIBI_window2];
    G2medianIBIm_window3 =  [output([output.FishGenotype]==2). mIBI_window3];
    G2medianIBIm_window4 =  [output([output.FishGenotype]==2). mIBI_window4];
    G2medianIBIm_window5 =  [output([output.FishGenotype]==2). mIBI_window5];
 
    G2nBoutRate_pre=[output([output.FishGenotype]==2). BoutRate];
    G2nBout=[output([output.FishGenotype]==2). nBout];
    G2nBoutRatewindow1= [output([output.FishGenotype]==2). BoutRate_window1];
    G2nBoutRatewindow2= [output([output.FishGenotype]==2). BoutRate_window2];
    G2nBoutRatewindow3= [output([output.FishGenotype]==2). BoutRate_window3];
    G2nBoutRatewindow4= [output([output.FishGenotype]==2). BoutRate_window4];
    G2nBoutRatewindow5= [output([output.FishGenotype]==2). BoutRate_window5];
    
     % Calculate median Amplitude 
    G2medianmedianAmpmitude_pre=[output([output.FishGenotype]==2). medianmedianAmplitude_pre];
    G2medianmaxAmpmitude_pre=[output([output.FishGenotype]==2). medianmaxAmplitude_pre];
    G2medianmeanAmpmitude_pre=[output([output.FishGenotype]==2). medianmeanAmplitude_pre];
    G2meanmaxAmpmitude_pre=[output([output.FishGenotype]==2). meanmaxAmplitude_pre]; 
    G2meanmeanAmpmitude_pre=[output([output.FishGenotype]==2).meanmeanAmplitude_pre];
    G2meanmedianAmpmitude_pre=[output([output.FishGenotype]==2). meanmedianAmplitude_pre];    
    
    G2routineTurnRatioHead =[output([output.FishGenotype]==2). routineTurnRatioHead];%heading
    G2routineTurnRatioTail = [output([output.FishGenotype]==2). routineTurnRatioTail];%tail
    
    
    
   % Rostral
    G1medianBoutDuration_pre = [output([output.FishGenotype]==1). medianBoutDuration_pre];
    G1meanBoutDuration_pre   = [output([output.FishGenotype]==1). meanBoutDuration_pre];
 
    G1medianNumOfOsc=[output([output.FishGenotype]==1). medianNumberOfOscillations];
    G1meanNumOfOsc=[output([output.FishGenotype]==1). meanNumberOfOscillations];
 
    
    G1medianBoutDistance_pre=[output([output.FishGenotype]==1). medianBoutDistance_pre];
    G1meanBoutDistance_pre=[output([output.FishGenotype]==1). meanBoutDistance_pre];
    G1TotalDistance = [output([output.FishGenotype]==1).TotalDistance];
    % Calculate Speed
    G1medianSpeed_pre=[output([output.FishGenotype]==1). medianSpeed];
    G1meanSpeed_pre=[output([output.FishGenotype]==1). meanSpeed];
  
    
    % Calculate  TBF 
    G1medianTBF_pre =       [output([output.FishGenotype]==1). medianTBF];
    G1meanTBF_pre =         [output([output.FishGenotype]==1). meanTBF];
    G1meanmaxTBF =          [output([output.FishGenotype]==1). meanmaxTBF];
    G1medianmaxTBF =          [output([output.FishGenotype]==1). medianmaxTBF];
    G1meanmedianTBF_pre =   [output([output.FishGenotype]==1). meanmedianTBF];
    G1meanmeanTBF_pre =     [output([output.FishGenotype]==1). meanmeanTBF];
    G1medianIBIm_pre =      [output([output.FishGenotype]==1). mIBI];
    
    G1medianIBIm_window1 =  [output([output.FishGenotype]==1). mIBI_window1];
    G1medianIBIm_window2 =  [output([output.FishGenotype]==1). mIBI_window2];
    G1medianIBIm_window3 =  [output([output.FishGenotype]==1). mIBI_window3];
    G1medianIBIm_window4 =  [output([output.FishGenotype]==1). mIBI_window4];
    G1medianIBIm_window5 =  [output([output.FishGenotype]==1). mIBI_window5];
 
    G1nBoutRate_pre=[output([output.FishGenotype]==1). BoutRate];
    G1nBout=[output([output.FishGenotype]==1). nBout];
    G1nBoutRatewindow1= [output([output.FishGenotype]==1). BoutRate_window1];
    G1nBoutRatewindow2= [output([output.FishGenotype]==1). BoutRate_window2];
    G1nBoutRatewindow3= [output([output.FishGenotype]==1). BoutRate_window3];
    G1nBoutRatewindow4= [output([output.FishGenotype]==1). BoutRate_window4];
    G1nBoutRatewindow5= [output([output.FishGenotype]==1). BoutRate_window5];
    
     % Calculate median Amplitude 
    G1medianmedianAmpmitude_pre=[output([output.FishGenotype]==1). medianmedianAmplitude_pre];
    G1medianmaxAmpmitude_pre=[output([output.FishGenotype]==1). medianmaxAmplitude_pre];
    G1medianmeanAmpmitude_pre=[output([output.FishGenotype]==1). medianmeanAmplitude_pre];
    G1meanmaxAmpmitude_pre=[output([output.FishGenotype]==1). meanmaxAmplitude_pre]; 
    G1meanmeanAmpmitude_pre=[output([output.FishGenotype]==1).meanmeanAmplitude_pre];
    G1meanmedianAmpmitude_pre=[output([output.FishGenotype]==1). meanmedianAmplitude_pre];    
    
    G1routineTurnRatioHead =[output([output.FishGenotype]==1). routineTurnRatioHead];%heading
    G1routineTurnRatioTail = [output([output.FishGenotype]==1). routineTurnRatioTail];%tail
    
   % Control
    G0medianBoutDuration_pre = [output([output.FishGenotype]==0). medianBoutDuration_pre];
    G0meanBoutDuration_pre   = [output([output.FishGenotype]==0). meanBoutDuration_pre];
 
    G0medianNumOfOsc=[output([output.FishGenotype]==0). medianNumberOfOscillations];
    G0meanNumOfOsc=[output([output.FishGenotype]==0). meanNumberOfOscillations];
 
    
    G0medianBoutDistance_pre=[output([output.FishGenotype]==0). medianBoutDistance_pre];
    G0meanBoutDistance_pre=[output([output.FishGenotype]==0). meanBoutDistance_pre];
    G0TotalDistance = [output([output.FishGenotype]==0).TotalDistance];
    % Calculate Speed
    G0medianSpeed_pre=[output([output.FishGenotype]==0). medianSpeed];
    G0meanSpeed_pre=[output([output.FishGenotype]==0). meanSpeed];

    % Calculate  TBF 
    G0medianTBF_pre =       [output([output.FishGenotype]==0). medianTBF];
    G0meanTBF_pre =         [output([output.FishGenotype]==0). meanTBF];
    G0meanmaxTBF =          [output([output.FishGenotype]==0). meanmaxTBF];
    G0medianmaxTBF =          [output([output.FishGenotype]==0). medianmaxTBF];
    G0meanmedianTBF_pre =   [output([output.FishGenotype]==0). meanmedianTBF];
    G0meanmeanTBF_pre =     [output([output.FishGenotype]==0). meanmeanTBF];
    G0medianIBIm_pre =      [output([output.FishGenotype]==0). mIBI];
    G0medianIBIm_window1 =  [output([output.FishGenotype]==0). mIBI_window1];
    G0medianIBIm_window2 =  [output([output.FishGenotype]==0). mIBI_window2];
    G0medianIBIm_window3 =  [output([output.FishGenotype]==0). mIBI_window3];
    G0medianIBIm_window4 =  [output([output.FishGenotype]==0). mIBI_window4];
    G0medianIBIm_window5 =  [output([output.FishGenotype]==0). mIBI_window5];
 
    G0nBoutRate_pre=[output([output.FishGenotype]==0). BoutRate];
    G0nBout=[output([output.FishGenotype]==0). nBout];
    G0nBoutRatewindow1= [output([output.FishGenotype]==0). BoutRate_window1];
    G0nBoutRatewindow2= [output([output.FishGenotype]==0). BoutRate_window2];
    G0nBoutRatewindow3= [output([output.FishGenotype]==0). BoutRate_window3];
    G0nBoutRatewindow4= [output([output.FishGenotype]==0). BoutRate_window4];
    G0nBoutRatewindow5= [output([output.FishGenotype]==0). BoutRate_window5];
 %%   
    G0Bout = sum(G0nBout);
    G1Bout = sum(G1nBout);
    G2Bout = sum(G2nBout);
    %%
     % Calculate median Amplitude 
    G0medianmedianAmpmitude_pre=[output([output.FishGenotype]==0). medianmedianAmplitude_pre];
    G0medianmaxAmpmitude_pre=[output([output.FishGenotype]==0). medianmaxAmplitude_pre];
    G0medianmeanAmpmitude_pre=[output([output.FishGenotype]==0). medianmeanAmplitude_pre];
    G0meanmaxAmpmitude_pre=[output([output.FishGenotype]==0). meanmaxAmplitude_pre]; 
    G0meanmeanAmpmitude_pre=[output([output.FishGenotype]==0).meanmeanAmplitude_pre];
    G0meanmedianAmpmitude_pre=[output([output.FishGenotype]==0). meanmedianAmplitude_pre];    
    
    G0routineTurnRatioHead =[output([output.FishGenotype]==0). routineTurnRatioHead];%heading
    G0routineTurnRatioTail = [output([output.FishGenotype]==0). routineTurnRatioTail];%tail
%

%% subplot all parameters

G0nBoutRate = [mean(G0nBoutRatewindow1) mean(G0nBoutRatewindow2) mean(G0nBoutRatewindow3) mean(G0nBoutRatewindow4) mean(G0nBoutRatewindow5)] ;
G1nBoutRate = [mean(G1nBoutRatewindow1) mean(G1nBoutRatewindow2) mean(G1nBoutRatewindow3) mean(G1nBoutRatewindow4) mean(G1nBoutRatewindow5)];
G2nBoutRate = [mean(G2nBoutRatewindow1) mean(G2nBoutRatewindow2) mean(G2nBoutRatewindow3) mean(G2nBoutRatewindow4) mean(G2nBoutRatewindow5)];

G0error = [std(G0nBoutRatewindow1) std(G0nBoutRatewindow2) std(G0nBoutRatewindow3) std(G0nBoutRatewindow4) std(G0nBoutRatewindow5)]/ sqrt(38);
G1error = [std(G1nBoutRatewindow1) std(G1nBoutRatewindow2) std(G1nBoutRatewindow3) std(G1nBoutRatewindow4) std(G1nBoutRatewindow5)]/sqrt(40);
G2error = [std(G2nBoutRatewindow1) std(G2nBoutRatewindow2) std(G2nBoutRatewindow3) std(G2nBoutRatewindow4) std(G2nBoutRatewindow5)]/sqrt(15);



h1 = figure(1);
title ('bout rate through 5 min'); hold on;
errorbar(window(2:6),G0nBoutRate,G0error , 'k' );
hold on;
errorbar((window(2:6)-600),G1nBoutRate,G1error,'r' );
hold on;
errorbar((window(2:6)+600),G2nBoutRate, G2error,'g' );

hold off;
%%
h2=figure(2); 

boxplot_3groupes(G0routineTurnRatioTail,G1routineTurnRatioTail,G2routineTurnRatioTail); hold on;
ylim([0 1]);
ylabel('Tail Turn ratio','FontSize',10);
hold off;

saveas(h2,['Parameters overview.fig'])

%%


h3=figure(3); 

subplot(1,6,1)
boxplot_3groupes(G0meanBoutDistance_pre,G1meanBoutDistance_pre,G2meanBoutDistance_pre);hold on;
ylim([0 4]);
ylabel('mean Distance (mm)','FontSize',10);
hold off;

subplot(1,6,2)
boxplot_3groupes(G0meanNumOfOsc,G1meanNumOfOsc,G2meanNumOfOsc);hold on;
 ylim([1 6]);
ylabel('mean Oscillations','FontSize',10);
hold off;

subplot(1,6,3)
boxplot_3groupes(G0meanBoutDuration_pre,G1meanBoutDuration_pre,G2meanBoutDuration_pre);hold on;
 ylim([0 0.3]);
ylabel('mean Duration (sec)','FontSize',10);
hold off;


subplot(1,6,4)
boxplot_3groupes(G0meanSpeed_pre,G1meanSpeed_pre,G2meanSpeed_pre);hold on;
ylim([0 15]);
ylabel('mean Speed (mm/sec)','FontSize',10);
hold off;

subplot(1,6,5)
boxplot_3groupes(G0meanTBF_pre,G1meanTBF_pre,G2meanTBF_pre); hold on;
 ylim([15 25]);
ylabel('mean TBF (Hz)','FontSize',10);
hold off;
subplot(1,6,6)
boxplot_3groupes(G0meanmedianAmpmitude_pre,G1meanmedianAmpmitude_pre,G2meanmedianAmpmitude_pre); hold on;
 ylim([0 15]);
ylabel('mean of median Amplitude(degree)','FontSize',10);
hold off;

saveas(h3,['Parameters to plot.fig'])