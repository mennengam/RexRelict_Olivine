% RexRelict_Olivine.m
%
% A script adapted to calculate recrystallized grain size in a deformed
% polycrystalline olivine aggregate containing both recrystallized and relict
% grains. Recrystallized grains (small, low strain) are separated from
% relict (large, high strain) grains using the grain orientation spread
% (GOS). GOS is the average 'mis2mean' value for each grain, where mis2mean
% is the misorientation between each pixel in a grain, and the mean
% orientation of that grain. (Cross et al., 2017)
%
% To use with different data (i.e. different phases) use the MTEX import
% wizard and replace lines 52-61.
%

clc
clear all
close all force
f = 1;  % f is a counter to keep figure numbers in order

%% Summer 2020 samples

 Name = 'Webster Addie';
 sampleName = 'Peterson_Olivine_3 NC16-19 Area 4 Montaged Data 4 Montaged Map Data.cpr';
% Name = 'Deposit no. 9';
% sampleName = 'Peterson_ultramafic_summer2020 NC17-27A Area 1 Montaged Data 2 Montaged Map Data.cpr'; % Compatible with both .ctf and .cpr/crc formats
% Name = 'Mincey Mine';
% sampleName = 'Peterson_ultramafic_summer2020 NC17-28A Area 2 Montaged Data 1 Montaged Map Data.cpr'; % Compatible with both .ctf and .cpr/crc formats
% Name = 'Moores Knob';
% sampleName = 'Peterson_ultramafic_summer2020 NC17-30 Area 3 Montaged Data 3 Montaged Map Data.cpr'; % Compatible with both .ctf and .cpr/crc formats
% Name = 'Dark Ridge';
% sampleName = 'Peterson_ultramafic_summer2020 NC17-31B Area 4 Montaged Data 4 Montaged Map Data.cpr'; % Compatible with both .ctf and .cpr/crc formats


%% Buck Creek Samples
% NC13-9B
% sampleName = 'Peterson - olivine NC13-9B Area 5 Montaged Data 3 Montaged Map Data.cpr';
% NC13-10
% sampleName = 'Peterson - olivine NC13-10 Area 2 Montaged Data 1 Montaged Map Data.cpr'
% NC15-5
%  sampleName = 'Peterson - olivine NC15-5 Area 6 Montaged Data 7 Montaged Map Data.cpr';
% NC16-3
% sampleName = 'Peterson_Olivine_2 NC16-3 Area 1 Montaged Data 1 Montaged Map Data.cpr';
% NC16-6A
% sampleName = 'Peterson_Olivine_3 NC16-6a Area 3 Montaged Data 3 Montaged Map Data.cpr'
% NC16-12
% sampleName = 'Peterson_Olivine_3 NC16-12 Area 2 Montaged Data 2 Montaged Map Data.cpr'
% NC16-16
% sampleName = 'Peterson - olivine NC16-16 Area 7 Montaged Data 8 Montaged Map Data.cpr' 
% NC16-17
% sampleName = 'Peterson_Olivine_3 NC16-17 Area 1 Montaged Data 1 Montaged Map Data.cpr'

%% Specify Crystal and Specimen Symmetries. To change phase use import_wizard

% crystal symmetry
CS = {... 
  crystalSymmetry('mmm', [4.8 10 6], 'mineral', 'Forsterite', 'color', [0.53 0.81 0.98])};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');


%% Read Files

% path to files
pname = 'MATLAB';

% which files to be imported
fname = [pname '\' sampleName];


%% Import the Data

% create an EBSD variable containing the data
if sampleName(end-2:end) == 'ctf'
    ebsd = EBSD.load(fname,CS,'interface','ctf',...
      'convertEuler2SpatialReferenceFrame');
elseif sampleName(end-2:end) == 'cpr'
    ebsd = EBSD.load(fname,CS,'interface','crc',...
      'convertEuler2SpatialReferenceFrame');
end


%% Reduce the map area (useful for debugging - reduces computation time)
% Use this to test code

% ebsd = ebsd(inpolygon(ebsd,[0 0 0.25*max(ebsd.x) 0.25*max(ebsd.y)]));


%% Calculate grains

rawebsd = ebsd;
ebsd = ebsd('indexed');

% Filter data to remove low-quality data
ebsd = ebsd(ebsd.mad <1.2);

% Construct grains using a critical misorientation of 15 degrees for olivine

[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd,'angle',15*degree);

% Remove wild spikes (1 pixel grains)
% - this step drastically reduces computation time, mostly w.r.t. twin merging
 ebsd(grains(grains.grainSize == 1)).phase = 0;
 ebsd = ebsd('indexed');

% Reconstruct grains without wild-spikes

[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd,'angle',15*degree);

%% Plot EBSD pixel data
% 
 figure(f); f=f+1;
 plot(ebsd,ebsd.orientations)
     hold on
 plot(grains.boundary)

%% GOS thresholding

% Use a trade-off curve to find the cutoff between low and high GOS values
GOS_knee = tradeOff(grains.GOS/degree);
        xlabel('Number of grains (cumulative)')
        ylabel('Grain Orientation Spread ( ^\circ )')


%% Separate out relict grains 

% Find IDs of merged grains with high GOS values (> cutoff)
ids = unique(grains(grains.GOS/degree > GOS_knee).id);

relictGrains = grains(ismember(grains.id,ids));

% All other grains are those with low GOS values (i.e. recrystallized grains)
rexGrains = grains(~ismember(grains.id,ids));


%% Plot mis2mean (hot colours represent high internal misorientation)

 figure
 plot(ebsd,ebsd.mis2mean.angle/degree)
     hold on
 plot(relictGrains.boundary) % Plot relict grain boundaries
     colorbar
    
        
%% Plot relict and recrystallized grains

figure(f); f=f+1;
plot(relictGrains,'facecolor',[1 0.5 0.5]) % Relict grains = red
    hold on
plot(rexGrains,'facecolor',[0.5 0.5 1]) % Recrystallized grains = blue
    hold on
plot(grains.boundary)% Plot grain boundaries
      
        
%% Find merged grains that are bisected by the map border
% - we want to remove these for the grain size analysis

face_id = grains.boundary.hasPhaseId(0);
bordergrain_id = grains.boundary(face_id).grainId;
bordergrain_id(bordergrain_id==0) = []; % Remove zeros

bordergrains = grains(ismember(grains.id,bordergrain_id));
nonbordergrains = grains(~ismember(grains.id,bordergrain_id));


%% Get area-equivalent circle diameters for all grains

d = 2*equivalentRadius(nonbordergrains);


%% Get area-equivalent circle diameters for relict and recrystallized grains

% Find relict and recrystallized grains that aren't bisected by the map border
relictNonBorder = relictGrains(ismember(relictGrains.id,nonbordergrains.id));
rexNonBorder = rexGrains(ismember(rexGrains.id,nonbordergrains.id));

% Area equivalent circle diameters for relict and recrystallized grains
relictD = 2*equivalentRadius(relictNonBorder);
rexD = 2*equivalentRadius(rexNonBorder);


%% Get grain size statistics for the recrystallized grains
% Contains code from Signal Processing Toolbox

amean_low = mean(rexD); % Arithmetic mean
gmean_low = 10^(mean(log10(rexD))); % Geometric mean
rmsmean_low = rms(rexD); % Root mean square (RMS)
median_low = median(rexD); % Median
mode_low = mode(rexD); % Mode

a1std_low = std(rexD); % 1 standard deviation


%% Plot grain size histograms

edges = (0:0.075:2.5); % Set the histogram bin widths
loglim = [0 2.5 0 0.25]; % Set the histogram axis limits
    
figure(f); f=f+1;
 set(gcf,'units','normalized','position',[0.15 0.15 0.7 0.5])


% Plot grain size distribution for all grains
% subplot(1,2,1),
% histogram(log10(d),edges,'Normalization','probability',...
%     'facecolor',[0.5 0.5 0.5]);
%     axis(loglim)
%      xlabel('Grain size (\mum)')
%     ylabel('Relative frequency (%)')
%     axis(loglim)
%     set(gca,'xtick',(0:2))
%     set(gca,'xticklabel',10.^get(gca,'xtick'),'yticklabel',100.*get(gca,'ytick'))
    
    
% Plot separate normalized grain size distributions for relict and recrystallized grains
% subplot(1,2,2),
% histogram(log10(relictD),edges,'Normalization','probability',...
%     'facecolor',[1 0.2 0.2]); % Relict grain size histogram (red)
%     hold on
% histogram(log10(rexD),edges,'Normalization','probability',...
%     'facecolor',[0.2 0.2 1]); % Recrystallized grain size histogram (blue)
%     xlabel('Grain size (\mum)')
%     ylabel('Relative frequency (%)')
%     axis(loglim)
%     set(gca,'xtick',(0:2))
%     set(gca,'xticklabel',10.^get(gca,'xtick'),'yticklabel',100.*get(gca,'ytick'))
    
%% Plot combined histograms of relict and rex grains

histogram(relictD,'BinWidth',100,'FaceColor','r','EdgeAlpha',0.75)
    hold on
histogram(rexD,'BinWidth',100,'FaceColor','b','EdgeAlpha',0.75)
xlabel('Grain size (\mum)')
ylabel('Grains')
[t,s] = title(Name,'Grain Size Distribution');
s.FontSize = 12;
legend({'Relict','Recrystallized'},'Location','east')

%% Grain Size Histograms

figure(f); f=f+1;
histogram(rexD)
xline(median_low,'-r',{'Median'})
xlabel('Grain size (\mum)')
ylabel('Grains')
[t,s] = title(Name,'Recrystallized Grain Size Distribution');
s.FontSize = 12;

figure(f); f=f+1;
histogram(relictD)
xline(median_low,'-r',{'Median'})
xlabel('Grain size (\mum)')
ylabel('Grains')
[t,s] = title(Name,'Relict Grain Size Distribution');
s.FontSize = 12;