function [OUT_p,REGS] = preproc_physio(physio_fn,save_dir,TR,nframes)
% Initial processing of physio data. Align with fMRI scan, extract
% "resp" (input to retro/rvhr), and initial heart beat detection
% (may need correcting later). Also write matlab-readable raw data
% ("TABLE"). All variables in OUT_p are aligned to fMRI.
% * TR: in msec
    
%
% Note: does not run correction on the fMRI data.
% Note: params are for biopac setup @ NIAAA 3T scanner @ NIH
% Note: must first convert *.acq files to *.svd with jacco's
% "biopac2svd.sh", see below
%
%
% version: 8/23/2017
% author:  catie
%

% First, must convert *.acq files to *.svd with jacco's code. Run
% this on the /DATA/<scandir>/phys/ folders, like:
%
%  ->  /raid/common/sleep/scripts/misc/biopac2svd.sh <path_to_dir>
%
% For some odd reason, it doesn't run for me, but did run for dante
% & jacco. Jacco figured out I need to run the following first:
% (see email from catie@tango, 8/21)
%    setenv IDL_PATH "+/misc/prog/rsi/idl:+/misc/imeel/dezwart/idl"
%    unsetenv IDL_STARTUP
  
% in case we need to run this manually in IDL for individual files:
% 
%fname="/raid/common/20170818_1/phys/catie2017-08-18T16_53_26.acq"
%d = read_biopac(fname,header=h)
%svd_file="test.svd"
%save_data,svd_file,d,header=h
    
% addpath('/misc/imeel/dezwart/matlab/');
% rmpath(genpath('~/matlab/chronux'));
% rmpath(genpath('~/matlab/plot_hht'));
% ----------------------------
% for testing
% ----------------------------
%physio_fn = 'sub_0006-mr_0015-eot_echo1.svd';

% ------------------------------------------ %
% parameters 
% ------------------------------------------ %
fs_phys = 2000; % Hz physio sampling (NIAAA)

% ------------------------------------------ %
% read files from IDL (svd) format,
% converted from biopac via Jacco's IDL script
% ------------------------------------------ %
if ~exist(physio_fn,'file')
    error('file not found');
end
[in_path,in_fname,ext] = fileparts(physio_fn);
data = load(physio_fn);
figure, plot(data, 'b');

% read physio file (NIAAA)
% ------------------------------------------ %
dt_phys = 1/fs_phys; % 2000 Hz 
    
% extract cardiac data & run peak detection
% after some initial filtering
% ------------------------------------------ %
card_dat = data(:,1);  % aligned with fMRI
% band-pass filter to help peak detection
fcut_BPF = [0.5,2];
Fn = fs_phys/2;
Wn = fcut_BPF/Fn; Nb = 2;
[B, A] = butter(Nb,Wn);
card_bpf =  filtfilt(B,A,double(card_dat));
% peak detection
card_rng = iqr(card_bpf); 
%[maxtab, mintab] = peakdet(card_bpf, 0.05*card_rng);
% card_trig_samples = maxtab(:,1);
minHeight = 0.05*card_rng;
minDist = 100+(fs_phys/2); % 04.01.17
[pks,locs] = findpeaks(card_bpf,'minpeakheight',minHeight,'minpeakdistance',minDist);
maxtab_c(:,1) = locs; maxtab_c(:,2) = pks;
% extract cardiac trigger times and approximate heart rate
card_trig_samples = locs;
card_trig_times = card_trig_samples*dt_phys;
% ibi & "instantaneous" hr
IBI = (diff(card_trig_times));
HR = (1./diff(card_trig_times))*60; %bpm

% check peak detection 
% ------------------------------------------ %
figure(3); clf; set(gcf,'color','w'); 
g1 = subplot(3,1,1); hold on;
plot(card_bpf); 
plot(maxtab_c(:,1),maxtab_c(:,2),'r.','markersize',20); % does well!
legend('cardiac data (band-pass filt)','detected beats');
xlabel('physio sample #');
% check heart rate as a function of time
g2 = subplot(3,1,2); hold on;
plot(card_trig_samples(1:end-1),HR,'b'); ylabel('beats per min'); title('cardiac rate');
hold on;
plot(card_trig_samples(1:end-1),HR,'c.'); % +dots
xlabel('physio sample #');
linkaxes([g1,g2],'x');
% IBI time series (by sample number - flag outliers)
g3 = subplot(3,1,3); hold on;
plot(IBI); 
xlabel('index'); ylabel('IBI');
drawnow;

% make output struct 
% ------------------------------------------ %
% note, all are aligned with fMRI
OUT_p.IBI_raw = IBI;
OUT_p.HR_raw = HR;
OUT_p.card_trig_times_s = card_trig_times; %sec
OUT_p.card_trig_samples = card_trig_samples;
OUT_p.card_dat = card_dat;
OUT_p.card_bpf = card_bpf;
OUT_p.card_pks = maxtab_c;
OUT_p.dt_phys = dt_phys;

% remove outlier beats?
% ------------------------------------------ %
%[OUT_p.IBI_clean,OUT_p.HR_clean,outliers] = despike_hr(OUT_p);
%OUT_p.outlier_beats = outliers;

% generate fMRI regressors
% ------------------------------------------ %
%REGS = build_fmri_regs(OUT_p,TR,nframes);

% save
% ------------------------------------------ %
if ~isempty(save_dir)
    if ~exist(save_dir,'dir')
        mkdir(save_dir);
    end
    new_name = [save_dir, '/', in_fname, '_physOUT.mat'];
    save(new_name,'OUT_p');
end





function [IBI_clean,HR_clean,outliers] = despike_hr(IN_p);
% for manually despiking 
    outliers = input('input vector of outliers: ');
    
    fig;
    subplot(211);
    IBI_clean = interp_ts(IN_p.IBI_raw, outliers, 1);
    subplot(212);
    HR_clean = interp_ts(IN_p.HR_raw, outliers, 1);

    
    
function REGS = build_fmri_regs(IN_p,TR,nframes)
% make some regressors sampled in TR bins
% hr: need smaller bins to see hrv?? check this.
    
    dt_phys = IN_p.dt_phys;
    resp = IN_p.resp;
    card_dat = IN_p.card_dat;
    card_bpf = IN_p.card_bpf;
    IBI_clean = IN_p.IBI_clean;
    
    % lookup table of IBI (denoised) v. cardiac trig
    % time (assigned to halfway between the respective beats)
    t_ibi = 0.5*(IN_p.card_trig_times_s(2:end) + IN_p.card_trig_times_s(1:end-1));
    assert(length(t_ibi)==length(IBI_clean))
    
    % sampling to match fMRI tr (center of each tr)
    Twin = 6; % sec windows (3s on either side)
    TR_s = TR/1000;
    t_fmri = (TR_s/2)+[0:TR_s:TR_s*nframes-TR_s];
    
    % make RV & HR regressors, as well as pulseox amplitude (stdev)
    rv = [];
    pa = [];
    hr = [];
    pa_bpf = [];
    for kk=1:nframes
        t = t_fmri(kk);
        
        % heart rate
        % ---------------------- %
        % get time bin centered at this TR
        t1 = max(0,t-Twin*0.5);
        t2 = min(TR_s*nframes,t+Twin*0.5);
        % find IBI's falling within this interval
        inds = intersect(find(t_ibi<=t2),find(t_ibi>=t1));
        hr(kk) = (60./median(IBI_clean(inds)));
        hrv(kk) = sqrt(mean(diff(IBI_clean(inds)).^2)); %rmssd
        
        % pulse amplitude
        % ---------------------- %
        if length(resp.wave)~=length(card_dat)
            error('resp & card sampled at different rates');
        else
            np = length(resp.wave);
        end
        % window (in samples)
        i1 = max(1,floor((t - Twin*0.5)/dt_phys)); 
        i2 = min(np, floor((t + Twin*0.5)/dt_phys));
        pa(kk) = std(card_dat(i1:i2));
        pa_bpf(kk) = std(card_bpf(i1:i2));
        % respiration variation
        % ---------------------- %
        rv(kk) = std(resp.wave(i1:i2));
    end
    
    % also convolve rv with RRF
    rv = rv-mean(rv);
    t = [0:TR:40-TR]; % 40-sec impulse response
    R = 0.6*(t.^2.1).*exp(-t/1.6) - 0.0023*(t.^3.54).*exp(-t/4.25); 
    R = R/max(R);
    rv_rrf = conv(rv,R);
    rv_rrf = rv_rrf(1:length(rv));
    % time derivative
    rv_rrf_d = diff(rv_rrf);
    rv_rrf_d = [rv_rrf_d(1), rv_rrf_d];
 
    % regressors for fmri
    % ---------------------- %
    REGS.rv = rv(:);
    REGS.rv_rrf = rv_rrf(:);
    REGS.rv_rrf_d = rv_rrf_d(:);
    REGS.pa = pa(:);
    REGS.pa_bpf = pa_bpf(:);
    REGS.hr = hr(:);
    REGS.hrv = hrv(:);
    
    fig;
    subplot(411); 
    plot(hr); title('heart rate');
    subplot(412);
    plot(hrv); title('heart rate variability (rmssd)');
    subplot(413);
    plot(pa); title('pulseox amplitude');
    hold on;
    plot(pa_bpf,'r'); title('pulseox amplitude - cardbpf');
    subplot(414);
    plot(rv); title('respiratory variation');
    xlabel('TR');