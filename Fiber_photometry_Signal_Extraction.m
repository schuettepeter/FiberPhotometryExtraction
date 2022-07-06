%% Run this script to extract fiber photometry signal from recording.
% (navigate to the 'sample recording' folder and run script.)

%% 
 dual_rec = 0;  

 minutes=[nan nan]; % enter minutes to be analyzed. minutes(1)==nan represents the first sample, and minutes(2)==nan is the last sample.
 %if minutes(n) isnt nan, then the minutes specified will be analyzed.
 %Example: minutes=[2 4], then minutes 2 to 4 will be analyzed
 %Example: minutes=[nan 5], then  minutes 0 to 5 will be analyzed
 %Example: minutes=[3 nan], then minutes 3 to end will be analyzed

 
 %% 
 
if dual_rec == 0
load('signal.mat');
end    

if dual_rec == 1
    if sigToUse == 1
        load('signal1.mat');
    elseif sigToUse == 2
        load('signal2.mat');
    end
end
 
 peak_rates=[];
 peak_vals=[];
    
 %Parameters for peak extraction
 num_std=3;
 min_peak_height=num_std;  %peaks larger than 3*std will be found
 Prominence=0.08;
 min_peak_distance=40;
 min_width=10;
 numfibers=1;

%x column 1: time
%x column 2: biobserve time pulse (every 5 seconds)
%x columns 3 to 7: safest to most dangerous zone in the environment.
%(safest is 3, and highre numbers are progressively more dangerous.)
%in the epm, 3,4 are closed arms, 5 is center, 6,7 are open arms
%in openfield, 3,4 are periphery and center, respectively
%in rat exposure, 3 is the farthest from the shock place, and 5 is the
%nearest. 6 is the shock location.

sig=sig';
ref=ref';

sig(1)=sig(2); sig(end)=sig(end-1);
ref(1)=ref(2);ref(end)=ref(end-1);

fs=1/10;
time=[0:fs:(length(sig)-1)*fs];

sig_norm=zeros(size(sig,1),size(sig,2));

a=1:length(time); %uses only the beginning of the recording for normalizing

f=find(sig(a)<mean(sig(a))+3*std(sig(a))); %does normalization ignoring points with very large peaks, which throw off the normalization
f=[a, f];

    p=polyfit(ref(1,f),sig(1,f),1);
    ref_scaled=ref(1,:)*p(1)+p(2);
    s = sig - ref_scaled;
    sig_norm=sig-ref_scaled+mean(sig);
    sig_norm=100*(sig_norm-mean(sig))/mean(sig);
    
    dfof=sig_norm; %now dfof is df/F in %
    

ref_norm=ref_scaled;
ref_norm=100*(ref_norm-mean(ref_norm))/mean(ref_norm);

min_peak_height=num_std*std(ref_norm(round(end/2):end));

% 
 f=find(dfof<-2);
 if length(f)>0.10*length(dfof)
     [counts,bins]=hist(dfof(f),50);
     baseline_correction=find(counts==max(counts));
     baseline_correction=bins(baseline_correction(1));
     dfof=dfof+abs(baseline_correction);
%      disp(['baseline correction=' num2str(baseline_correction)]);
 end
 
a=minutes(1);
b=minutes(2);

if isnan(minutes(1))
    a=1;
else
    a=minutes(1)*10*60;
end

if isnan(minutes(2))
b=length(time);
else
b=minutes(2)*10*60;
end

time2=time;
time=time(a:b);
dfof=dfof(a:b);
ref_norm=ref_norm(a:b);
 
 
[peak_val,peak_ind, peak_width, ~]=findpeaks(dfof,'minpeakheight',min_peak_height,'minpeakdistance',min_peak_distance,'minpeakwidth',min_width,'minpeakprominence',Prominence);%'threshold',0.001);

if dual_rec == 0
save('signal_extracted', '-regexp', '^(?!(folders|jj|nn|ans|a|b|i|burrow_norat_xp|burrow_xp|EPM|miniscope|mouse_toyrat|other_vid|shockgrid_xp|simpleRat_xp|fear_cond_chamber|setByParentScript)$).')
end

if dual_rec == 1
    if sigToUse == 1
        save('signal_extracted1', '-regexp', '^(?!(folders|jj|nn|ans|a|b|i|burrow_norat_xp|burrow_xp|EPM|miniscope|mouse_toyrat|other_vid|shockgrid_xp|simpleRat_xp|fear_cond_chamber|setByParentScript)$).')
    elseif sigToUse == 2
        save('signal_extracted2', '-regexp', '^(?!(folders|jj|nn|ans|a|b|i|burrow_norat_xp|burrow_xp|EPM|miniscope|mouse_toyrat|other_vid|shockgrid_xp|simpleRat_xp|fear_cond_chamber|setByParentScript)$).')
    end
end