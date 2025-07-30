% wmMGS1.m
%
% 528 0.75 s TRs (396 s)
% 
% behavioral stimulus presentation script for single-item 'mapping' MGS
% experiment (fMRI data from this used to estimate WM encoding model and/or
% Bayesian generative model)
%
% participant sees a single item, chosen from ~16 or 8 positions, makes
% memory-guided eye movement after 12 s delay (10 or 15 TRs). Slight
% modifications of ssPri's mapping task (that was 40 unique positions).
% Shooting for ~ 6 min runs if possible. 
%
% offset should be one of: 5.625/16.875 <for staggered runs/32 pos> OR
% 11.25 if just doing 16 positions offset from meridians
% CMRR 7T: offset for 8 locations is 22.5 with jitter -5.625 +5.625
%
% updated TCS 9/13/2017 - reacts to variable TR (for mb experiments, p.TR
% should usually be some # of TRs - here, 4x0.75=3; could also use shorter,
% but not all that necessary pending how task timing is set up)
% try:
% subj = 'test';ses=1;run = 5;offset = 22.5;wmMGS1_7T(subj,ses,run,offset)
%
% check TR trigger alignment with important expt event:
% rem(p.delay_start - p.expt_start, 0.75) the first delay locked to TR
% < 0.05 should be fine
% "trial should have started" marker: p.behind_by
% expt time recorded: p.end_expt-p.expt_start
%
% TCS 6/9/2017
%%%%%%%%%%%% CMRR 7T MGS %%%%%%%%%%%
% max eccentrity 9.04 x 20.34 deg.
% make a ellipse aperture
% TODO: stop recording w/ ESC
% TODO: fixation as dot within circle

function wmMGS1_7T(subj,ses,run,offset)

%debug
% subj = 'test';
% run = 0;
% offset = 0;
debugg = 1;
p.do_et = 0;
p.scanner = 1;
% screenRes = [1920 1080]; % in the behav room [1920 1080] NYU scanner room [1280 1024]
%%
p.expt_name = 'wmMGS1';
p.TR = 3; % 4x multiband, so measured TR is 0.75, but "TR" for stim is 3

p.subj = subj;
p.run = run;
p.ses = ses;
p.offset = offset;

if ~exist('./data','dir')
    mkdir('./data');
end

p.filename = sprintf('./data/%s_ses-%02.f_run-%02.f_task-%s.mat',p.subj,p.ses,p.run,p.expt_name);
if p.do_et == 1
    p.eyedatafile = sprintf('%s_M%02.f',p.subj(1:min(length(p.subj),3)),p.run);
end

p.rng_seed = cputime*1000;
rng(p.rng_seed);

%% ------ size of relevant stim features, etc ------ %
% the viewable area is about 9.04 x 20.34 deg, ecc should be in radius
p.wm_ecc_long = 9;     % deg, radius
p.wm_ecc_short = 3.5;
p.cue_size = 0.55; % deg
p.wm_size = 0.65/2;  % deg, size of WM dots half of the original for CMRR7T
p.aperture_size_hl = 10; % [or max ecc?] radius, half of 20.34
p.aperture_size_vl = 4.5; % radius, half of 9.04

%% EXP CONDITION
p.n_pos = 8; % and we'll use 2 offset angles
%%
p.fix_size_in  = 0.075; % radius, deg
p.fix_size_out = 0.30; % radius, deg
p.fix_pen = 1.5;

p.fix_size_mult = 1.1;%1.25; % for pre/post experiment, fix is bigger

% now we draw out (frameoval), cue (if end of trial), in (filloval) for fix

% ------- keyboard stuff --------------------------- %
if ismac == 1
    p.esc_key = KbName('escape'); % press this key to abort
else
    p.esc_key = KbName('esc');
end
p.start_key = [KbName('5%') KbName('5')];  % should be lower-case t at prisma? (or %5, which is top-row 5, or 5, which is numpad 5)
p.space = KbName('space');

% ------ color of relevant stim features ----------- %
darkgrey = 100*[1 1 1]; midgrey = 128*[1 1 1]; black = [1 1 1]*49; 
white = 255*[1 1 1];
p.bg_color  = darkgrey;%20*[1 1 1]
p.fix_color = white;%midgrey; blackgrey75*[1 1 1];% whitegrey[150 150 150];        % during trial/delay/etc

p.wm_color = white;%p.fix_color;

p.go_color = 0.8*midgrey;%p.fix_color;% 130*[1 1 1];%[255 255 255]; % when subj should choose, color of fix

p.dim_amt = 0.58;%0.85;%0.8; % multiply fix_color by this during ITI, start, end periods

% ------ conditions ------ %
% p.r_cond = [1]; % 1: R1, 2: R2_cued, 3: R2_choose
p.repetitions = 2; % 2 * 8 = 16 trials 
p.ntrials = p.repetitions * p.n_pos;

% this is the only "condition" for this experiment - so only saving this,
% no "conditions" field

p.wm_ang = linspace(0,360-360/p.n_pos,p.n_pos).' + p.offset; % CARTESIAN!!
p.wm_ang = repmat(p.wm_ang, p.repetitions, 1);
% (TODO: if we change # of positions and increase # of repetitions, which
% is unlikely, would need to slightly rework this)
% random within each repetition cycle
rnd_idx_1st = randperm(p.n_pos);
rnd_idx_2nd = randperm(p.n_pos)+8;
p.rnd_idx = [rnd_idx_1st rnd_idx_2nd];
p.wm_ang = p.wm_ang(p.rnd_idx, :);
p.jitter = (p.offset/4)*(-1 + 2.*rand(p.ntrials,1));
p.wm_ang = p.wm_ang + p.jitter;
 
%% ------ timing of trial events --------- %
% ITIs should be int_num * TR + 0.9 + 0.9 so that beginning of delay is
% locked to TR
p.warn_dur = 1.0;  % brightening of fix to 'warn' subject target is about to appear [due to long ITI]
p.targ_dur = 0.5;
p.delay_dur = 12;
p.cue_dur = 0.7; % as in previous reports [NOTE shorter than the wmChoose task]
p.feedback_dur = 0.8; %1.3; %0.8; 
 
% we know that n_trials should be divisible by 4 (quadrants), so we'll do
% 25% short, 50% medium, 75% long [2, 3, 4 TRs at 3 s TR 
% ORIGINAL
%p.itis = 0+p.TR*[2*ones(p.ntrials/4,1); 3*ones(p.ntrials/2,1); 4*ones(p.ntrials/4,1) ];
 
% DEPENDENT UPON TR: (fills from end of previous trial to next TR -
% p.TR-(p.feedback_dur+p.cue_dur), then adds n_TRs, then fills to warning
% cue (p.TR-(p.warn_dur+p.targ_dur)
p.itis = p.TR-(p.feedback_dur+p.cue_dur) + ...
         p.TR*[1*ones(p.ntrials/4,1); 2*ones(p.ntrials/2,1); 3*ones(p.ntrials/4,1) ] + ...
         p.TR-(p.warn_dur+p.targ_dur);
% % Add 2s more to the original ITI, 4.5, 7.5, 9.5
% p.itis = p.TR-(p.feedback_dur+p.cue_dur) + ...
%          p.TR+p.TR*[1*ones(p.ntrials/4,1); 2*ones(p.ntrials/2,1); 3*ones(p.ntrials/4,1) ] + ...
%          p.TR-(p.warn_dur+p.targ_dur);
p.itis = p.itis(randperm(length(p.itis)));     

% -------- timing of experiment events ------- %
p.start_wait = 1 * p.TR + p.TR-(p.warn_dur+p.targ_dur);%2.4 + 0.9; % after first trigger [after dummys]
p.end_wait = (p.TR-(p.feedback_dur+p.cue_dur))+2*p.TR; %0.9 + 4.8 + 0.6; % means minimum of 14.8 s after last trial event pre-ITI [adjusting to be integer # TRs!!!!]

p.trial_dur = p.warn_dur + p.targ_dur + p.delay_dur + p.cue_dur + p.feedback_dur + p.itis;
p.exp_dur = p.start_wait + p.end_wait + p.ntrials * (p.warn_dur + p.targ_dur + p.delay_dur + p.cue_dur + p.feedback_dur) + sum(p.itis);

% ------- things to save -------- %
%% wm target coordinates
p.wm_coords = [p.wm_ecc_long p.wm_ecc_short] .* [cosd(p.wm_ang) sind(p.wm_ang)];

% being a bit redundant here...
p.trial_start = nan(p.ntrials,1);
p.targ_start  = nan(p.ntrials,1);
p.delay_start = nan(p.ntrials,1);
p.resp_start  = nan(p.ntrials,1);
p.feedback_start = nan(p.ntrials,1);
p.iti_start   = nan(p.ntrials,1);
p.trial_end   = nan(p.ntrials,1);

p.behind_by = nan(p.ntrials,1); % keep track of timing errors, they should be basically tiny

% ------- Screen setup, optics --------- %

if p.scanner == 0
    p.resolution = [1920 1080];%[1280 1024][1280 800]; % desired resolution, compare to the 'actual' resolution below [this is for mbp laptop]
    p.screen_height = 37.5; %39.29; %37.5; % cm
    p.viewing_distance = 57; %184; %57; % cm

else
    p.resolution = [1920 1080]; % scanner
%     p.resolution = screenRes; % scanner at CMRR
    p.screen_height = 39.29;%39.29; %35; % cm
    p.viewing_distance = 184; %184;%63; % cm
end
p.refresh_rate = 120;
if p.scanner == 1
    %p.screen_height = p.screen_height * 1024/1080;
end

p.screen_width = p.screen_height * p.resolution(1)/p.resolution(2); %
p.screen_height_deg = 2*atan2d(p.screen_height/2,p.viewing_distance);
p.screen_width_deg  = 2*atan2d(p.screen_width/2, p.viewing_distance);
p.ppd = p.resolution(2)/p.screen_height_deg;  % used to convert rects, positions later on
p.center = p.resolution/2;  % could do offset centers, etc?

if debugg
    [w, p.scr_rect] = Screen('OpenWindow',max(Screen('Screens')),[0 0 0], ...
        [0 0 p.resolution]);
else
    [w, p.scr_rect] = Screen('OpenWindow',max(Screen('Screens')),[0 0 0]); HideCursor;
end

Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference','TextRenderer',1);

% check resolution
if p.scr_rect(3) ~= p.resolution(1) || p.scr_rect(4) ~= p.resolution(2)
    disp('Check resolution!!!');
    Screen('CloseAll');
    ShowCursor;
    return;
end

p.ifi = Screen('GetFlipInterval',w);
%(p.aperture_rect(3)-p.aperture_rect(1))/p.ppd is 2*radius
%p.aperture_rect = CenterRectOnPoint([0 0 2 2]*p.ppd*p.aperture_size,p.center(1),p.center(2));
aperutre_rect = [0 0 2 2].*[0 0 p.aperture_size_hl p.aperture_size_vl]*p.ppd;
p.aperture_rect = CenterRectOnPoint(aperutre_rect,p.center(1),p.center(2));
p.fix_rect_out  = CenterRectOnPoint([0 0 2 2] * p.ppd  * p.fix_size_out,p.center(1),p.center(2));
p.fix_rect_in   = CenterRectOnPoint([0 0 2 2] * p.ppd  * p.fix_size_in, p.center(1),p.center(2));

% --------- eyetracking ----------- %
if p.do_et == 1
    if p.scanner == 1
        Eyelink('SetAddress','192.168.1.5')
    end
    
    el=EyelinkInitDefaults(w);
    
    el.backgroundcolour=p.bg_color(1);  % TODO: fix this?
    el.calibrationtargetcolour=p.fix_color(1); % white
    el.calibrationtargetsize=2*2*p.wm_size; %%%%%% 1.3% percentage of the screen
    el.calibrationtargetwidth=0.5; %%% inner area 1% as percentage of screen
    el.msgfontcolour=p.fix_color(1);
    p.foregroundcolour=p.fix_color(1);
    % 192.168.1.5
    % sca
    EyelinkUpdateDefaults(el);

    Eyelink('Initialize','PsychEyelinkDispatchCallback') % initialises the eyetracker
   % SCANNER: right eye!!!!!!
    Eyelink('command','calibration_type=HV9'); % updating number of callibration dots
    
    %Customize calibration point
    width = 1920; height = 1080; p2p_dist=350; %%%% CHANGE for scanner????
    Eyelink('command', 'generate_default_targets = NO');
    Eyelink('command','calibration_samples = 9');
    Eyelink('command','calibration_sequence = 0,1,2,3,4,5,6,7,8');
    Eyelink('command','calibration_targets = %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d',...
        width/2,height/2,  width/2+p2p_dist,height/2,  width/2-p2p_dist,height/2,  width/2,height/2+p2p_dist,  width/2,height/2-p2p_dist,...
        width/2+p2p_dist,height/2+p2p_dist,  width/2+p2p_dist,height/2-p2p_dist,  width/2-p2p_dist,height/2-p2p_dist,  width/2-p2p_dist,height/2+p2p_dist);
    Eyelink('command','validation_samples = 9');
    Eyelink('command','validation_sequence = 0,1,2,3,4,5,6,7,8');
    Eyelink('command','validation_targets = %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d',...
        width/2,height/2,  width/2+p2p_dist,height/2,  width/2-p2p_dist,height/2,  width/2,height/2+p2p_dist,  width/2,height/2-p2p_dist,...
        width/2+p2p_dist,height/2+p2p_dist,  width/2+p2p_dist,height/2-p2p_dist,  width/2-p2p_dist,height/2-p2p_dist,  width/2-p2p_dist,height/2+p2p_dist);
    
    s=Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,AREA');% (,GAZERES,HREF,PUPIL,STATUS,INPUT');
    s=Eyelink('command', 'sample_rate = 1000');
    s=Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
    s=Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
    
    % make sure that we get gaze data from the Eyelink
    
    %------ calibrate the eye tracker --------
    EyelinkDoTrackerSetup(el);
    if s~=0
        error('link_sample_data error, status: ',s)
    end
    Eyelink('openfile',p.eyedatafile);
%     if p.scanner == 1
%         Eyelink('SetAddress','192.168.1.5')
%     end
%     
%     el=EyelinkInitDefaults(w);
%     
%     el.backgroundcolour=p.bg_color(1);  % TODO: fix this?
%     el.calibrationtargetcolour=p.fix_color(1);
%     el.calibrationtargetsize=2*p.wm_size;
%     el.calibrationtargetwidth=1;
%     el.msgfontcolour=p.fix_color(1);
%     p.foregroundcolour=p.fix_color(1);
%      % 192.168.1.5
%     % sca
%     EyelinkUpdateDefaults(el);
% 
%     
%     Eyelink('Initialize','PsychEyelinkDispatchCallback') % initialises the eyetracker
%    % SCANNER: right eye!!!!!!
%     Eyelink('command','calibration_type=HV13'); % updating number of callibration dots
%     s=Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,AREA');% (,GAZERES,HREF,PUPIL,STATUS,INPUT');
%     s=Eyelink('command', 'sample_rate = 1000');
%     s=Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
%     s=Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
%     
%     
%     
%     % make sure that we get gaze data from the Eyelink
%     
%     
%     
%     %------ calibrate the eye tracker --------
%     EyelinkDoTrackerSetup(el);
%     if s~=0
%         error('link_sample_data error, status: ',s)
%     end
%     Eyelink('openfile',p.eyedatafile);

end

% save key for XDAT labels
p.XDAT_labels = {'warn','targets','delay','response','feedback','ITI'};

draw_aperture = @() Screen('FillOval',w,p.bg_color,p.aperture_rect);

Screen('FillRect',w,[0 0 0]);
draw_aperture();

% NO instructions during scanning (big visual transient when they
% disappear)
%txt = 'Remember dot position precisely';
%Screen('TextSize', w, 30);
%DrawFormattedText(w,txt,'center',p.center(2)-4*p.ppd,p.fix_color);

%Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
%Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

Screen('DrawDots',w,[0;0],p.fix_size_mult*p.fix_size_out*p.ppd*2+p.fix_pen,p.dim_amt*p.fix_color,p.center,2);
Screen('DrawDots',w,[0;0],p.fix_size_mult*p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2);

Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.dim_amt*p.fix_color,p.center,2);
%Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
Screen('Flip',w);

% check for esc, space.... 

resp = 0;
while resp == 0
    [resp, ~] = checkForResp(p.start_key, p.esc_key);
    if resp == -1
        Screen('CloseAll'); ShowCursor;
        Eyelink('ShutDown'); 
        return;
    end
end
clear resp;
p.expt_start = GetSecs;

% blank screen - fix back to normal size
Screen('FillRect',w,[0 0 0]);
draw_aperture();
Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.dim_amt*p.fix_color,p.center,2); 
Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.dim_amt*p.fix_color,p.center,2); 
%Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);

Screen('Flip',w);

if p.do_et == 1
    Eyelink('Message','xDAT %i', 10);
    Eyelink('StartRecording'); % make 1 big edf file (save time)
end

% ------ initial wait time (check for esc) ---------

resp = 0;
while (GetSecs-p.expt_start) < p.start_wait
    [resp, ~] = checkForResp([], p.esc_key);
    if resp == -1
        sca;
        fprintf('ESC pressesd during post-trigger wait time\n');
        ShowCursor;
        return;
    end
end
clear resp;

for tt = 1:p.ntrials
    fprintf('Trial %i: start\n',...
                tt);
    % warning cue (XDAT 1) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % NOTE: I'm making this a "trial should have started" marker, to keep
    % timing from accumulating errors & to keep below code more readable
    trial_start = p.expt_start + p.start_wait + sum(p.trial_dur(1:(tt-1)));
    p.trial_start(tt) = GetSecs;
    p.behind_by(tt) = p.trial_start(tt)-trial_start;
    
    if p.do_et == 1

        Eyelink('Message','TarX %s', num2str(0));
        Eyelink('Message','TarY %s', num2str(0));    
        Eyelink('Message','xDAT %i',1); 
        Eyelink('command', 'record_status_message "TRIAL %d of %d"', tt, p.ntrials);
        
    end

    while GetSecs < trial_start + p.warn_dur
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        % fixation
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
                
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
        
    end
        
    % targets (XDAT 2) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('Trial %i: wm target: %f\n',...
                tt,p.wm_ang(tt,1));
    p.targ_start(tt) = GetSecs;
        
    if p.do_et == 1
        
        Eyelink('Message','xDAT %i',2);
        
    end
       
    while GetSecs < trial_start + p.warn_dur + p.targ_dur
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        % target 1
        Screen('DrawDots',w,p.ppd*[1;-1].*p.wm_coords(tt,:).', p.wm_size*p.ppd, p.wm_color, p.center, 2);
        
        % fixation
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        if isnan(p.targ_start(tt))
            p.targ_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
        
    end
    
    % delay (XDAT 3) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('Trial %i: delay 12 s\n',...
                tt);
    if p.do_et == 1
        Eyelink('Message','xDAT %i',3);    
    end
        
    while GetSecs < trial_start + p.warn_dur + p.targ_dur + p.delay_dur
        
        % draw aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
%         if p.conditions(tt,1) ~= 3
%             Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.choose_color,p.center,2);
%         else
%             Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.targ_colors{1}(tt,:),p.center,2);
%         end
        
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        if isnan(p.delay_start(tt)) % only do this once
            p.delay_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
        
    end
    
    % go cue (XDAT 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('Trial %i: respond\n',...
                tt);
    if p.do_et == 1
        Eyelink('Message','TarX %s', num2str(p.wm_coords(tt,1)));
        Eyelink('Message','TarY %s', num2str(p.wm_coords(tt,2)));
        Eyelink('Message','xDAT %i',4);
    end
    
    while GetSecs < trial_start + p.warn_dur + p.targ_dur + p.delay_dur + p.cue_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.go_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

        Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.go_color,p.center,2);
        %Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.bg_color,p.center,2); 

        Screen('Flip',w);
        
        if isnan(p.resp_start(tt))
            p.resp_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
    end
    
    % feedback (XDAT 5, tarx, tary) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('Trial %i: feedback\n',...
                tt);
    if p.do_et == 1
        Eyelink('Message','TarX %s', num2str(p.wm_coords(tt,1)));
        Eyelink('Message','TarY %s', num2str(p.wm_coords(tt,2)));
        % NOTE: incorrect on 1/3 trials
        Eyelink('Message','xDAT %i',5);
        
    end
    
    while GetSecs < trial_start + p.warn_dur + p.targ_dur + p.delay_dur + p.cue_dur + p.feedback_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        % target 1
        Screen('DrawDots',w,p.ppd*[1;-1].*p.wm_coords(tt,:).', p.wm_size*p.ppd, p.wm_color, p.center, 2);
        
        % fixation
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.go_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

        Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.go_color,p.center,2);
        %Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.bg_color,p.center,2); 
        
        %Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 

        Screen('Flip',w);
        
        if isnan(p.feedback_start(tt))
            p.feedback_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
 
    end
    
    % ITI (XDAT 6) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('Trial %i: Inter Trial Interval %f\n',...
                tt, p.itis(tt));
    if p.do_et == 1
        Eyelink('Message','TarX %s', num2str(0));
        Eyelink('Message','TarY %s', num2str(0));
        
        Eyelink('Message','xDAT %i',6);
    end

    % TODO: base this on sum of all events up til then? or do that
    % everywhere?
    while GetSecs < trial_start + p.warn_dur + p.targ_dur + p.delay_dur + p.cue_dur + p.feedback_dur + p.itis(tt)
        
        % draw aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.dim_amt*p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.dim_amt*p.fix_color,p.center,2); 
        
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        
        Screen('Flip',w);
        
        if isnan(p.iti_start(tt))
            p.iti_start(tt) = GetSecs;
            % save [note: in scanner, do this at beginning of ITI after first flip]
            save(p.filename,'p');
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');

            save(p.filename,'p');
            return;
        end
        
    end
      
end

% ------- wait for p.end_wait -------- %
end_tmp = GetSecs;

resp = 0;
while (GetSecs-end_tmp) < p.end_wait
    [resp, ~] = checkForResp([], p.esc_key);
    if resp == -1
        sca;
        fprintf('ESC pressesd during post-experiment wait time\n');
        ShowCursor;
        return;
    end
end
clear resp;
p.end_expt = GetSecs;

save(p.filename,'p');

% END OF EXPERIMENT - TEXT
if p.do_et == 1
    Eyelink('Message','xDAT %i',11);
end

Screen('FillRect',w,[0 0 0]);
draw_aperture();
txt = sprintf('End of run %i',p.run);
DrawFormattedText(w,txt,'center',p.center(2)-4*p.ppd,p.fix_color);
%Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);

Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.dim_amt*p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.dim_amt*p.fix_color,p.center,2); 

Screen('Flip',w);

resp = 0;
while resp == 0
    [resp, ~] = checkForResp(p.space, p.esc_key);
end
clear resp;

if p.do_et == 1
    Eyelink('StopRecording');
    Eyelink('ReceiveFile',[p.eyedatafile '.edf'],[p.eyedatafile '.edf']);
    
    p.eyedatafile_renamed = [p.filename(1:(end-3)) 'edf'];
    movefile([p.eyedatafile '.edf'],p.eyedatafile_renamed);
    
    Eyelink('ShutDown');
end

Screen('CloseAll');
ShowCursor;

return