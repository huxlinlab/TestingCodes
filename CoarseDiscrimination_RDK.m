%%%USE THIS CODE OF COARSE MAPPING/TRAINING%%%
clear all; close all;
oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
Screen('Preference', 'SkipSyncTests', 1);
KbName('UnifyKeyNames');
 
%-----Subject Settings------------------ 
H_ecc_stim                      = -17;      %  (degs, neg is left)
V_ecc_stim                      = -5;      % Vertical stimulus ecc (degs, neg is down)
  
initials                        = 'NBP20181212';
exp_type                        = 'Coarse Discrimination';                  % Newsome for abocve noise types, tania for direction ranges, both for combined, PTM for orientation
n_staircases                    = 1;                                        % How many staircases you want to run per condition. 2 should be fine
n_trials                        = 100;                                      % Number of trials per staircase
correct_trials                  = 0;     
V_ecc_stim                      = -V_ecc_stim; 
H_ecc_fix                       = 0;                                         % Horizontal fixation ecc (degs, neg is left)
V_ecc_fix                       = 0;                                        % Vertical fixation ecc (degs, neg is up)
  

%-----Dot Settings-- --------------------
stimulus_duration               = 500;                                      % ms
aperture_radius                 = 2.5;                                      % Degrees
ITI                             = 1.5;
dot_density                     = 3.5;                                      % Dots per square degree, 1.7 (Newsome&Pare '88)
dot_size                        = 14;                                       % diamteter, arcmin
dot_color                       = 0;                                        % Grayscale units
dot_speed                       = 5;                                        % deg/s
dot_lifetime                    = 250;                                      % in ms, for direction range only
angle_set                       = 1;                                        % Remember to set this 0 will be horizontal, 1 will be vertical
background                      = 128;                                      % Grayscale Units
linearize                       = 1;

%%%% subjective rating %%%%
lowConf = 'a'; medConf  = 's'; highConf = 'd';
dataOri            = zeros(n_trials*n_staircases,1);
dataResp                        = zeros(n_trials*n_staircases,1); 
dataAccu                        = zeros(n_trials*n_staircases,1);
dataConf                        = zeros(n_trials*n_staircases,1);

%-----Rig Settings----------------------
viewing_dist                    = 42;
screen_width                    = 40;                                 %in cm
resolution                      = [800 600];
theta                           = atand((screen_width/2)/viewing_dist);
scale_factor                    = theta*60/(resolution(1)/2);   
%scale_factor                    = 3.025;                                    % Arcmin/pixel ~2 for 76 cm, 3.49 at 42cm
frame_rate                      = 120;                                   % Screen frame rate (hz)
ISI                             = (1/frame_rate)*1000;                                   % ms, for best results, round to multiples of frame rate
  
%-----Eyetracking----------------------
ET                              = 1;                                    %Flag to run with or without eyetracking. 0 no, 1 yes
pref_eye                        = 2;                                    %Which eye to pick, 0 left, 1 right, 2 both 
dummymode                       = 0;                                    %Flag to use the mouse as false eyetracking, turn off if eyetracking!
fix_box_size                    = 4;                                    %Size of fixation box in degrees, will force a break if fixation leaves this window
asymmetric_fix_box_size         = [0 2 2 3 2]; %Use a rectangle instead of square for fixation, 1st val is yes or no, 2-5 are to define pos and neg height and width
asym_adder                      = 0; %adds an extra degree to a specific dimension of fixation box, will be removed soon for better implimentation
edfFile                         = initials;
ET_res                          = [1280 1024];                             %Resolution that the ET screen runs at (should be 4:3, 1280x1024 is native to monitor
first_run                       = 1; %just a  marker of whether we've run yet or not
get_bg                          = 0; %tells the ET program to grab a screenshot of the stimuls on the first run and transfer it to the ET computer for reference
show_fix_box                    = 0;
broke_fixation_count            = 0;
Cal_ecc                         = 8; %in degrees, how far out do we want the calibration targets? Only if using custom
use_small_cal                   = 1; %Use smaller, inverted targets, may get better fixations

%-----Housekeeping----------------------
% Scale things based on viewing distance, and convert other stuff to
% the units PsychToolbox wants...
h_ecc_orig                      = H_ecc_stim;
v_ecc_orig                      = -V_ecc_stim; 
results                         = zeros(n_trials*n_staircases,5);
ndots                           = round(aperture_radius^2*pi*dot_density);
aperture_radius                 = 60*aperture_radius/scale_factor;
stimulus_radius                 = aperture_radius;
H_ecc_stim                      = H_ecc_stim*60/scale_factor;
H_ecc_fix                       = H_ecc_fix*60/scale_factor;
V_ecc_stim                      = V_ecc_stim*60/scale_factor;
V_ecc_fix                       = V_ecc_fix*60/scale_factor;
dot_step                        = dot_speed*60/scale_factor/(1000/ISI);
dot_size                        = dot_size/scale_factor;
mv_length                       = round(stimulus_duration/(1000/frame_rate));
age{1,mv_length}                = [];
age(1:mv_length)                = {zeros(1, ndots)};
positions{1,mv_length}          = [];
postions(1:mv_length)           = {zeros(1, ndots)};

% Conversions for direction ranges
Radius                          = aperture_radius;
Velocity                        = dot_step;
lifetime                        = dot_lifetime/(1000/frame_rate);

%-----Randomize Trial Order-------------
total_trials                    = n_trials*n_staircases;
perm                            = randperm(total_trials);
perm                            = mod(perm,n_staircases)+1;

%-----Spatial Envelope------------------
bps                             = (stimulus_radius)*2+1;

%-----Open Screens----------------------
screens                         = Screen('Screens');
screenNumber                    = max(screens);
[w, rect]                       = Screen('OpenWindow',screenNumber,0,[],[],2);
Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
screen_rect                     = Screen('Rect',w);
Screen('FillRect',w, background);
Screen('Flip', w);
Screen('FillRect',w, background);
Screen('TextSize',w,20);
Screen('TextFont',w,'Charcoal');

%-----Stimulus Rectangles---------------
movie_rect                      = [0,0,bps,bps];
scr_left_middle                 = fix(screen_rect(3)/2)-round(bps/2);
scr_top                         = fix(screen_rect(4)/2)-round(bps/2);
screen_rect_middle              = movie_rect + [scr_left_middle, scr_top, scr_left_middle, scr_top];
screen_patch                    = screen_rect_middle+[H_ecc_stim,V_ecc_stim,H_ecc_stim,V_ecc_stim];
sr_hor                          = round(screen_rect(3)/2);
sr_ver                          = round(screen_rect(4)/2);
fix_hor                         = sr_hor+H_ecc_fix;
fix_ver                         = sr_ver+V_ecc_fix;
stim_hor                        = sr_hor+H_ecc_stim;
stim_ver                        = sr_ver+V_ecc_stim;
fix_rect                        = SetRect(0, 0, 4*scale_factor, 4*scale_factor); %set fix size here
fix_rect                        = CenterRectOnPoint(fix_rect,fix_hor,fix_ver);
noise                           = [0 40 80 120 160 200 240 280 320 355];
stair1                          = 1;
stair2                          = 3;
stair3                          = 5;
staircount1                     = 0;
staircount2                     = 0;
staircount3                     = 0;
stair_array                     = zeros(1,n_trials*n_staircases);
for x = 1:n_staircases
    for i = 1:n_trials
        stair_array(i+n_trials*(x-1)) = x;
    end
end
stair_array = Shuffle(stair_array);

if linearize
%     if projector
%         screen_clut = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
%         screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
%     else
        addpath('~/Desktop/Matlab Maintenance');
        load 'hpp12130cal.mat'
        Screen('LoadNormalizedGammaTable',screenNumber,gTmp);
%     end
end

if ET       
    fix_box_size = fix_box_size*60/scale_factor;
    asymmetric_fix_box_size(2:5)= asymmetric_fix_box_size(2:5)*60/scale_factor;
    if asymmetric_fix_box_size(1)==0
        fix_box = SetRect(-.5*fix_box_size,-0.5*fix_box_size,.5*fix_box_size,0.5*fix_box_size);
        fix_box = CenterRectOnPoint(fix_box,fix_hor,fix_ver);
        fix_box(4)=fix_box(4)+asym_adder*60/scale_factor;
    else
        fix_box = SetRect(-.5*fix_box_size,-0.5*fix_box_size,.5*fix_box_size,0.5*fix_box_size);
        for i=1:4
            changeit = zeros(1,4);
            if fix_box(i)~=asymmetric_fix_box_size(i+1) || fix_box(i)~=asymmetric_fix_box_size(i+1)
                changeit(i)=1;
            end
        end
        fix_box = CenterRectOnPoint(fix_box,fix_hor,fix_ver);
    end
    %Initiates Eyelink defaults with details about the graphical environment

    el=EyelinkInitDefaults(w);

    if use_small_cal
%         el.backgroundcolour = BlackIndex(el.window);
        el.backgroundcolour = WhiteIndex(el.window);
        
% el.msgfontcolour  = WhiteIndex(el.window);
    el.msgfontcolour  = BlackIndex(el.window);    
        el.imgtitlecolour = WhiteIndex(el.window);
                el.calibrationtargetcolour = BlackIndex(el.window);

%         el.calibrationtargetcolour = WhiteIndex(el.window);
        el.calibrationtargetsize= 1;
        el.calibrationtargetwidth=0.5;
        % call this function for changes to the calibration structure to take
        % affect
        EyelinkUpdateDefaults(el);
    end


    %Init Eyelink connection, exit if no link exists
    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted. \n');
    end

    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.

    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end


    % open file to record data to
    i = Eyelink('Openfile',initials(1:8) );
    if i~=0
        fprintf('Cannot create EDF file ''%s'' ', edffilename);

        return;
    end;

    % make sure we're still connected.
    if Eyelink('IsConnected')~=1 && ~dummymode
        dummymode=1;
        return;
    end;

    Eyelink('command', 'add_file_preamble_text ''Recorded by EyelinkToolbox for Huxlin lab exp''');

    %Set up mapping of gaze position from tracker to screen pixel positions for
    %fixation

    Eyelink('command', 'screen_pixel_coords = %ld %ld %ld %ld', 0, 0, screen_rect(3)-1, screen_rect(4)-1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, screen_rect(3)-1, screen_rect(4)-1);


    % set calibration type.
    Eyelink('command', 'calibration_type = HV9');
    % you must send this command with value NO for custom calibration
    % you must also reset it to YES for subsequent experiments
    Eyelink('command', 'generate_default_targets = YES');
    Eyelink('command', 'generate_default_targets = YES');
    shc=rect(3)/2;
    svc=rect(4)/2;
    % STEP 5.1 modify calibration and validation target locations
    Eyelink('command','calibration_samples = 5');
    Eyelink('command','calibration_sequence = 0,1,2,3,4,5');
    Eyelink('command','calibration_targets = %d %d,%d %d, %d %d, %d %d, %d %d,' , shc,svc, shc,svc+Cal_ecc, shc-Cal_ecc,svc, shc+Cal_ecc,svc, shc,svc-Cal_ecc);
    Eyelink('command','validation_samples = 5');
    Eyelink('command','validation_sequence = 0,1,2,3,4,5');
    Cal_ecc=Cal_ecc-100;

    Eyelink('command','validation_targets = %d %d,%d %d, %d %d, %d %d, %d %d,', shc,svc, shc,svc+Cal_ecc, shc-Cal_ecc,svc, shc+Cal_ecc,svc, shc,svc-Cal_ecc);         
    %set EDF file contents
    [v, vs] = Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.]n', vs);
    vsn = regexp(vs,'\d','match');
    % remote mode possible add HTARGET ( head target)
    if v == 3 && str2double(vsn{1}) == 4 % if EL 1000 and tracker version 4.xx
        % remote mode possible add HTARGET ( head target)
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT,HTARGET');
        % set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT,HTARGET');
    else
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT');
        % set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    end

    if ~dummymode
        % Hide the mouse cursor and setup the eye calibration window
        Screen('HideCursorHelper', w);
    end

    WaitSecs(0.5);
    Eyelink('Command', 'set_idle_mode');
    Eyelink('command','clear_screen 1');

    WaitSecs(0.5);
    Beeper;

    %Let's get started! camera setup mode, calibration and validation
    EyelinkDoTrackerSetup(el);
    WaitSecs(0.1)
    Beeper;
   
end

Screen('DrawText',w,strcat('Dots:   ',exp_type),100,100,0);

if angle_set == 1
   Screen('DrawText',w,'Use the LEFT/RIGHT arrows to respond',100,130,0);
else
   Screen('DrawText',w,'Use the UP/DOWN arrows to respond',100,130,0);
end

Screen('DrawText',w,[int2str(total_trials),'  trials - press any key to start'],100,160,0);
Screen('Flip',w);
KbWait(-1,2);
validKey = 0;

trial = 1;
TSTART = tic;
while trial<total_trials+1
    aa = GetSecs;
    Screen('FillRect',w, background);
    Screen('Flip', w);
    Screen('FillOval',w,0,fix_rect);
    %Screen('DrawText',w,[num2str(trial) '/' num2str(total_trials)], resolution(1)/2, resolution(2)*3/4, 0);
    Screen('Flip',w);    
   
    if ET
        Eyelink('Message', 'Fixation cross drawn');
    end
    
    if staircount1 == 3
        stair1 = stair1 + 1;
        if stair1 == 11
            stair1 = 10;
        end
        staircount1 = 0;
    end
    if staircount2 == 3
        stair2 = stair2 + 1;
        if stair2 == 11
            stair2 = 10;
        end
        staircount2 = 0;
    end
    if staircount3 == 3
        stair3 = stair3 + 1;
        if stair3 == 11
            stair3 = 10;
        end
        staircount3 = 0;
    end
    which_stair = stair_array(trial);

    if ET
        %Send start trial message
        Eyelink('Message', 'TRIALID %d', trial');
        Eyelink('command', 'record_status_message "TRIAL %d/%d"', trial,total_trials');
        %Clear tracker display and draw box at center
        Eyelink('Command', 'set_idle_mode');
        Eyelink('Command', 'clear_screen %d', 1)
        WaitSecs(0.55)
        Eyelink('StartRecording');
        eye_used = pref_eye; 
    end 

    direction = ceil(2*rand);
    if direction==1
        correct = 'LeftArrow';
        incorrect = 'RightArrow';
    else
        correct = 'RightArrow';
        incorrect = 'LeftArrow';
    end
    
    if which_stair == 1
        dir_rangeP = noise(stair1);
    elseif which_stair == 2
        dir_rangeP = noise(stair2);
    elseif which_stair == 3
        dir_rangeP = noise(stair3);
    end
    
    dir_range = dir_rangeP;
    angle=direction*180;           

    %-----Position
    positions{1}(:,1) = (rand(ndots,1)-.5)*bps;
    positions{1}(:,2) = (rand(ndots,1)-.5)*bps;

    for i=1:ndots
        while sqrt(positions{1}(i,1)^2+positions{1}(i,2)^2)>stimulus_radius
            positions{1}(i,:) = [ceil((rand-0.5)*bps),ceil((rand-0.5)*bps)];
        end
    end

    vectors = pi*(angle+(rand(ndots,1)-.5)*dir_range)/180;

    %-----Lifetime
    age{1} = ceil(lifetime*rand(ndots,1));

    for j=2:mv_length
        %----------Update Dots----------
        for i=1:(ndots);
            %-----Move to new positions, wrap if necessary
            positions{j}(i,1) = positions{j-1}(i,1)+Velocity*cos(vectors(i));
            positions{j}(i,2) = positions{j-1}(i,2)+Velocity*sin(vectors(i));

            %-----Age dots
            age{j}(i) = age{j-1}(i)+1;
            %-----Kill and respawn dead dots
            if age{j}(i)>lifetime
                positions{j}(i,1) = (rand-.5)*bps;
                positions{j}(i,2) = (rand-.5)*bps;
                age{j}(i) = 1;
            end

            while sqrt(positions{j}(i,1)^2+positions{j}(i,2)^2)>stimulus_radius
                positions{j}(i,:) = [ceil((rand-0.5)*bps),ceil((rand-0.5)*bps)];
            end
        end
    end

    %Finish the ITI
    WaitSecs(ITI-(GetSecs-aa));
    Screen('FillOval',w,0,fix_rect);
    %Screen('DrawText',w,[num2str(trial) '/' num2str(total_trials)], resolution(1)/2, resolution(2)*3/4, 0);
    Screen('Flip',w);

    priorityLevel=MaxPriority(w);
    Priority(priorityLevel);

    if ET
        Eyelink('Message', 'ITI Finished');
    end
    % Wait for fixation
    if ET
        Eyelink('StartRecording')
        Eyelink('Message', 'Waiting for Fixation');
        quitflag = 0;
        while ~quitflag
            if Eyelink( 'NewFloatSampleAvailable') > 0 % in dummy mode use mousecoordinates
                sample = Eyelink('NewestFloatSample');
            else
                [x,y,button] = GetMouse(w);
                sample.gx=x;
                sample.gy=y;
            end
            if eye_used<2
                if IsInRect(sample.gx(eye_used+1),sample.gy(eye_used+1),fix_box)
                    sample.gx
                    sample.gy
                    1;
                    quitflag = 1;
                end
            else
                if IsInRect(sample.gx(1),sample.gy(1),fix_box) && IsInRect(sample.gx(2),sample.gy(2),fix_box)
                    quitflag = 1;
                end
            end
            WaitSecs(.25);
        end
        Eyelink('Message', 'Successful Fixation');
    end 
    fixated = 1;
    if ET
        Eyelink('Message', 'Start Stimulus');
    end
    Beeper(1000); 
    % Play the movie
    i=1;
    while i <= mv_length
        Screen('FillRect',w, background);
        if exist('positions_n','var')
            Screen('DrawDots', w, transpose(positions_n{i}),dot_size,dot_color,[stim_hor stim_ver],2);
        end
        Screen(w,'DrawDots',transpose(positions{i}),dot_size,dot_color,[stim_hor stim_ver],2);

        % Draw the black fixation cross
        Screen('FillOval',w,0,fix_rect);
        %Screen('DrawText',w,[num2str(trial) '/' num2str(total_trials)], resolution(1)/2, resolution(2)*3/4, 0);
        Screen('Flip',w);
        i=i+1;
        if i>mv_length
            break
        end
        if ET
            if Eyelink( 'NewFloatSampleAvailable') > 0 % in dummy mode use mousecoordinates
                sample = Eyelink('NewestFloatSample');
            else
                [x,y,button] = GetMouse(w);
                sample.gx=x;
                sample.gy=y;
            end
            if eye_used<2
                          if length(sample.gx) ==1
                    if IsInRect(sample.gx, sample.gy,fix_box) == 0
                        Eyelink('Message', 'Fixation Broken');
                        Beeper;
                        Beeper;
                        Beeper;
                        fixated=0;
                        break;
                    end
                else
                    if IsInRect(sample.gx(eye_used+1),sample.gy(eye_used+1),fix_box) == 0
                        Eyelink('Message', 'Fixation Broken');
                        Beeper;
                        Beeper;
                        Beeper;
                        fixated=0;
                        break;
                    end
                end
            else
                if IsInRect(sample.gx(1),sample.gy(1),fix_box) ==0 || IsInRect(sample.gx(2),sample.gy(2),fix_box) == 0
                    Eyelink('Message', 'Fixation Broken');
                    Beeper;
                    Beeper;
                    Beeper;
                    fixated=0;
                    break;
                end
            end
        end
    end    
    Screen('FillOval',w,0,fix_rect);
    %Screen('DrawText',w,[num2str(trial) '/' num2str(total_trials)], resolution(1)/2, resolution(2)*3/4, 0);
    Screen('Flip',w);
    Priority(0);

    % Get the response
    if ET && fixated==1
        Eyelink('Message', 'Successfuly Fixated Trial');
    end
    if fixated==0
        WaitSecs(0.1);
        Eyelink('StopRecording');
        broke_fixation_count = broke_fixation_count+1; 
        if broke_fixation_count > 9
            EyelinkDoTrackerSetup(el);
            broke_fixation_count = 0;
        end
    else
        validKey = 0;
        WaitSecs(0.1);
        tic
        if ET
            Eyelink('StopRecording');
        end
        while ~validKey
            [secs, keyCode, deltaSecs] = KbWait(-1);
            if keyCode(KbName(incorrect))
                validKey = 1;
                keyResp  = KbName(incorrect); 
                rs = 0;
                if which_stair == 1
                    results(trial,2) = noise(stair1);
                    stair1 = stair1 - 1;
                    staircount1 = 0;
                    if stair1 == 0
                        stair1 = 1;
                    end
                elseif which_stair == 2
                    results(trial,2) = noise(stair2);
                    stair2 = stair2 - 1;
                    staircount2 = 0;
                    if stair2 == 0
                        stair2 = 1;
                    end
                elseif which_stair == 3
                    results(trial,2) = noise(stair3);
                    stair3 = stair3 - 1;
                    staircount3 = 0;
                    if stair3 == 0
                        stair3 = 1;
                    end
                end
            elseif keyCode(KbName(correct))
                rs = 1;
                keyResp  = KbName(correct);  
                correct_trials = correct_trials+1;
                validKey = 1;
                results(trial,5) = 1;
                if which_stair == 1
                    staircount1 = staircount1 + 1;
                    results(trial,2) = noise(stair1);
                elseif which_stair == 2
                    staircount2 = staircount2 + 1;
                    results(trial,2) = noise(stair2);
                elseif which_stair == 3
                    staircount3 = staircount3 + 1;
                    results(trial,2) = noise(stair3);
                end
            end
        end

        %%% CONFIDENCE
%         FlushEvents('keyDown');
%         confKey = 0;
%         WaitSecs(0.1);
%         while ~confKey
%             [secs, keyCode, deltaSecs] = KbWait;
%             if keyCode(KbName(lowConf))
%                 cf=1;
%                 confKey = 1;
%             elseif keyCode(KbName(medConf))
%                 cf=2;
%                 confKey = 1;
%             elseif keyCode(KbName(highConf))
%                 cf=3;
%                 confKey = 1;
%             end
%         end
%         dataOri(trial)  = KbName(correct);
%         dataResp(trial) = keyResp;
%         dataAccu(trial) = rs;
%         dataConf(trial) = cf;
        if rs==0
            Beeper(600,.5,.1);Beeper(400,.5,.1);
        elseif rs==1
            Beeper(800,.5,.1);Beeper(1000,.5,.1);
        end
        reaction_time = toc;
        results(trial,1) = trial;
        results(trial,3) = direction;
        results(trial,4) = reaction_time;
        trial = trial+1;
        Screen('FillOval',w,0,fix_rect);
        %Screen('DrawText',w,[num2str(trial) '/' num2str(total_trials)], resolution(1)/2, resolution(2)*3/4, 0);
        Screen('Flip',w);
    end
end   
length_of_session = toc(TSTART);
Screen('CloseAll');
Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);

%mean(dataConf)
noise(stair1)
correct_trials/n_trials

% Save all data.
cd('~/Desktop/Subjects')
folder = exist(initials,'file');
if folder ~= 7
    mkdir(initials)
end 
cd(initials);
time = date;
save(strcat(initials,'_coarse_discrim_subj_',int2str(h_ecc_orig),'_',int2str(v_ecc_orig),'_',time));
cd('~/Desktop/Updated Matlab Codes')
