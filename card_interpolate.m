function card_interpolate(phys_name)
% Load physio function into workspace
phys_fn_0 = load(phys_name);
phys_fn = phys_fn_0.OUT_p;
phys_fix = phys_fn;

% Detect min peaks - they aren't in the OUT_p yet
card_rng = iqr(phys_fn.card_bpf); 
minHeight = 0.05*card_rng;
minDist = 100+((1/phys_fn.dt_phys)/2);
[phys_min_raw,loc] = findpeaks(-phys_fn.card_bpf,'minpeakheight',minHeight,...
    'minpeakdistance',minDist);

% Make sure the min peaks are synched with the max peaks
phys_min = zeros(length(phys_min_raw),1);
phys_min(1:length(phys_min)) = nan;
new_loc = zeros(length(phys_min_raw),1);
new_loc(1:length(new_loc)) = nan;
skip = 0;
phys_min(1) = phys_min_raw(1);
trav_cond = min(length(phys_fn.card_trig_samples), length(loc));
for i = 2:trav_cond
    while i+skip <= length(phys_min_raw) && loc(i+skip) < phys_fn.card_trig_samples(i-1)
        skip = skip+1;
    end
    if i+skip > length(phys_min_raw)
        break;
    end
    phys_min(i) = phys_min_raw(i+skip);
    new_loc(i) = loc(i+skip);
end

mark = 0;

% Run over the signal to check validity
for i = 1:length(phys_fn.card_trig_samples) - 2
    % Partition the data, then find outliers with IQR in each segment
    if mod(i, 200) == 1 && (length(phys_fn.card_trig_samples) - i) >= 200
        low = i;
        high = i + 199;
        if high > length(phys_fn.card_trig_samples)
            high = length(phys_fn.card_trig_samples);
        end
        
        % IBI_seg = phys_fn.IBI_clean(low:(high - 1));
        IBI_seg = phys_fn.IBI_raw(low:(high - 1));
        int_mu = mean(IBI_seg);
        int_quantiles = quantile(IBI_seg, 3);
        int_IQR = int_quantiles(3) - int_quantiles(1);
        int_up_bound = int_quantiles(3) + 1.5 * int_IQR;
        int_low_bound = int_quantiles(1) - 1.5 * int_IQR;
        
        peak_seg = phys_fn.card_dat(phys_fn.card_trig_samples(low:high));
        peak_quantiles = quantile(peak_seg, 3);
        peak_IQR = peak_quantiles(3) - peak_quantiles(1);
        peak_up_bound = peak_quantiles(3) + 2 * peak_IQR;
        peak_low_bound = peak_quantiles(1) - 2 * peak_IQR;
        
        min_seg = phys_min(low:high);
        min_quantiles = quantile(min_seg, 3);
        min_IQR = min_quantiles(3) - min_quantiles(1);
        min_up_bound = min_quantiles(3) + 2 * min_IQR;
        min_low_bound = min_quantiles(1) - 2 * min_IQR;
    end
    
    % Check whether a peak is an outlier, by its height and its interval
    if check_validity(phys_fn, phys_min, int_up_bound, int_low_bound, peak_up_bound,...
            peak_low_bound, min_up_bound, min_low_bound, i)
        if mark == 0 && i > 10
            mark = i;
        end     
    else
        % The end of a bad region, fill up the signal with the info
        % provided by peaks around it
        if mark ~= 0 && ~check_validity(phys_fn, phys_min, int_up_bound, ...
                int_low_bound, peak_up_bound, peak_low_bound, min_up_bound, ...
                min_low_bound, i+1) && ~check_validity(phys_fn, phys_min, int_up_bound, ...
                int_low_bound, peak_up_bound, peak_low_bound, min_up_bound, ...
                min_low_bound, i+2) && ~check_validity(phys_fn, phys_min, int_up_bound, ...
                int_low_bound, peak_up_bound, peak_low_bound, min_up_bound, ...
                min_low_bound, i+3)
            rep = round((phys_fn.card_trig_samples(i+1)-phys_fn.card_trig_samples(mark-3))/(int_mu/phys_fn.dt_phys));
            if rep == 0
                rep = 1;
            end
            std_int = ceil((phys_fn.card_trig_samples(i+1)-phys_fn.card_trig_samples(mark-3))/rep);
            overlap_seg = compute_seg(-phys_min, phys_fn, mark-3, i+1, rep, std_int);
            
            % Plot the interpolated signal for user to check
            figure;
            if i+10 < length(phys_fn.card_trig_samples)
                plot_length = i+10;
            else
                plot_length = length(phys_fn.card_trig_samples);
            end
            plot(phys_fn.card_bpf(phys_fn.card_trig_samples(mark-10):phys_fn.card_trig_samples(plot_length)), 'b');
            hold on;
            plot(overlap_seg(phys_fn.card_trig_samples(mark-10):phys_fn.card_trig_samples(plot_length)), 'r');
            x = input('Keep this change? (y/n)', 's');
            
            % If the user wants to keep the auto fix, replace the original
            % signal
            if (strcmp(x, 'y') || strcmp(x, 'Y'))
                phys_fix.card_bpf(phys_fn.card_trig_samples(mark-2):phys_fn.card_trig_samples(i)-rep) =...
                    overlap_seg(phys_fn.card_trig_samples(mark-2):phys_fn.card_trig_samples(i)-rep);
            end
            close all;
            mark = 0;
        end
    end
end

% Save the signal onto a new MAT file
card_rng = iqr(phys_fix.card_bpf); 
minHeight = 0.05*card_rng;
minDist = 100+((1/phys_fix.dt_phys)/2);
[~,final_max] = findpeaks(phys_fix.card_bpf,'minpeakheight',minHeight,...
    'minpeakdistance',minDist);
OUT_p = {};
OUT_p.dt_phys = phys_fix.dt_phys;
OUT_p.card_dat = phys_fix.card_dat;
OUT_p.card_bpf = phys_fix.card_bpf;
OUT_p.card_trig_samples = final_max;
OUT_p.card_trig_times_s = final_max*phys_fix.dt_phys;
OUT_p.IBI_clean = diff(OUT_p.card_trig_times_s);
[fp,fn,ext] = fileparts(phys_name);
fp = extractBefore(fp, '/preprocessed');

save([fp,'/fixed/',fn,'_fixed',ext],'OUT_p');
end

% The function used to check whether a peak is an outlier
% INVALID = TRUE
function is_invalid = check_validity(phys_fn, phys_min, int_up_bound, int_low_bound, peak_up_bound,...
            peak_low_bound, min_up_bound, min_low_bound, index)
is_invalid = (index <= length(phys_fn.IBI_raw) && ...
            (phys_fn.IBI_raw(index) > int_up_bound || ...
            phys_fn.IBI_raw(index) < int_low_bound)) || ...
            (index <= length(phys_min) && ~isnan(phys_min(index))&& ...
            (phys_min(index) > min_up_bound || phys_min(index) < min_low_bound)) || ...
            (phys_fn.card_dat(phys_fn.card_trig_samples(index)) > peak_up_bound...
            || phys_fn.card_dat(phys_fn.card_trig_samples(index)) < peak_low_bound);

end

% Interpolate the segment that contains bad data
function overlap_seg = compute_seg(phys_min, phys_fn, head, tail, rep, std_int)
tmp_seg = zeros(std_int, 1);
weight_arr = [0.046, 0.272, 0.682];
overlap_seg = zeros(length(phys_fn.card_bpf), 1);

% Extract signals on the left (Proven to be valid)
for i = 1:3
    component = phys_fn.card_bpf(phys_fn.card_trig_samples(head-4+i):...
        phys_fn.card_trig_samples(head-3+i));
    component = resize(component, std_int);
    tmp_seg = tmp_seg + weight_arr(i)*component;
end

% Equalize the end
tmp_seg = end_equalize_complete(tmp_seg(1), tmp_seg(1), tmp_seg);

% Overlay the waves as much as needed
for i = 1:rep
    overlap_seg((phys_fn.card_trig_samples(head)+(i-1)*std_int):...
        (phys_fn.card_trig_samples(head)+i*std_int-1)) = tmp_seg;
end

% Detect the minimums for convenience
minDist = 100+((1/phys_fn.dt_phys)/2);

[~,max_loc] = findpeaks(overlap_seg,'minpeakdistance',minDist);
[~,min_loc] = findpeaks(-overlap_seg,'minpeakdistance',minDist);

if (abs(length(max_loc)-length(min_loc)) > 2)
    disp('Error with previous data, skipping location');
    overlap_seg(:) = nan;
    overlap_seg(max_loc(1):(max_loc(length(max_loc))-1)) = phys_fn.card_bpf(max_loc(1):(max_loc(length(max_loc))-1));
    return
end

% Interpolate the max and the min
max_value_arr = lin_interpolate(phys_fn.card_bpf(phys_fn.card_trig_samples(head)), ...
    phys_fn.card_bpf(phys_fn.card_trig_samples(tail)), rep);
min_value_arr = lin_interpolate(phys_min(head), phys_min(tail), rep);

% Connect the ends seamlessly
for i = 1:length(min_loc)
    overlap_seg(max_loc(i):min_loc(i)) = end_equalize_half(max_value_arr(i), ...
    min_value_arr(i+1), overlap_seg(max_loc(i):min_loc(i)));
    
    overlap_seg(min_loc(i)+1:max_loc(i+1)-1) = end_equalize_half(min_value_arr(i+1), ...
    max_value_arr(i+1), overlap_seg(min_loc(i)+1:max_loc(i+1)-1));
end

overlap_seg(1:max_loc(1)-1) = nan;
overlap_seg(max_loc(length(max_loc)):length(overlap_seg)) = nan;

end

% Resize the signal wave by linearly interpolating, similar to image resize
function new_seg = resize(segment, newSize)
seg_size = size(segment);
seg_size = seg_size(1);
new_seg = interp1(1:seg_size,segment(:),linspace(1,seg_size,newSize));
new_seg = reshape(new_seg, [newSize, 1]);
end

function new_seg = end_equalize_complete(fore, post, segment)
head_diff = fore - segment(1);
new_seg = segment + head_diff;
overall_diff = post - new_seg(length(new_seg));
for i = 1:length(new_seg)
    new_seg(i) = new_seg(i) + (overall_diff*(double(i)/length(new_seg)));
end
end

% Adjust the function smoothly so that its ends match the correct value
function new_seg = end_equalize_half(fore, post, segment)
head_diff = fore - segment(1);
new_seg = segment + head_diff;
init_val = new_seg(1);
scale_val = (post - init_val)/(new_seg(length(new_seg)) - init_val);
new_seg = (new_seg - init_val) * scale_val + init_val;
end

% Can use another version - exp_interpolate
function value_arr = lin_interpolate(head, tail, rep)
value_arr = zeros(rep+1, 1);
overall_diff = tail - head;
for i = 1:length(value_arr)
    value_arr(i) = head + double(overall_diff)*(i-1)/rep;
end
end