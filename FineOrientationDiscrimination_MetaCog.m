%% Third draft Fine Orientation Discrimination task with Meta-Cognitive assessment 03/21/24
% Hard coded to 3 staircases
% Required Subfunctions:
% presentScore - displays trial count, current score, fixation point
% dujeEnvelope - Duje's version NOT matlab default
% Eyetracking suite: startStimulusEyelink, setupEyelink,
% duringStimulusEyelink
% genImagesDT - creates mask
% checkScreenSettings

% Suppresses warnings native to Psychtoolbox, shuffle trials, other housekeeping
oldVisualDebugLevel             = Screen('Preference', 'VisualDebugLevel', 3);
oldSupressAllWarnings           = Screen('Preference', 'SuppressAllWarnings', 1);

Screen('Preference', 'SkipSyncTests', 1);
KbName('UnifyKeyNames');
rng('shuffle')                                          % Shuffles trials to make each session unique

%% Identify starting location
initials                        = 'CB048092625';        % Identify patient for saving purposes
1
H_ecc_stim                      = -4;                    % X Location in degrees of target stim
V_ecc_stim                      = 5;  
% Y Location in degrees of target stim

ET                              = 0;                    % Flag to run with or without eyetracking. 0 no, 1 yes

% Shift fixation point if necessary (rarely used)
H_ecc_fix                       = 0;                                                       
V_ecc_fix                       = 0;

%% Initialize variables
% Monitor information
viewing_dist                    = 42;                                       % in cm
screen_width                    = 54.3;                                     % in cm
resolution                      = [1920 1080];
theta                           = atand((screen_width/2)/viewing_dist);
scale_factor                    = theta*60/(resolution(1)/2);               % in arcmin/px
frame_rate                      = 120;                                       % Screen frame rate (hz)
fix_shift                       = (H_ecc_fix*60)/scale_factor;
windowRect                      = [];

% Beautification
font                            = 'Arial';
fontSize                        = 40;
fix_size                        = 8;
background                      = 128;                                      % greyscale units

% Open Screens to be used in session
screens                         = Screen('Screens');
screenNumber                    = max(screens);                             % always display to secondary monitor if available
[w, rect]                       = Screen('OpenWindow',screenNumber,0,windowRect,[],2);
screen_rect                     = Screen('Rect',w);

% Check that monitor resolution and frame rate are accuate. DOES NOT WORK
% FOR SCREEN WIDTH AT THIS TIME
% checkScreenSettings(resolution, frame_rate, screen_width, screenNumber);

% Linearize display with pre-generated, monitor specific Gamma Table
addpath('~/Desktop/Matlab Maintenance');
load 'Alienware_room2_082922_cubic.mat';  % CHANGE TO COLOR TABLES!!!
Screen('LoadNormalizedGammaTable', screenNumber, gTmp);

Screen(w, 'BlendFunction', GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('FillRect', w, background);
Screen('Flip', w);
Screen('TextSize', w, fontSize/2);
Screen('TextFont', w, font);

% Other task parameters
delayMask                       = .25;                  % Indicated length of time mask should be delayed for in Seconds (.25 is standard)
n_trials                        = 90;                   % Trials per staircase
bonusTrials                     = 30;                   % # of additional, non-staircase trials. Should be ~10% of total trial count
ITI                             = 1;                    % Time between trials in seconds

% Set up presentation angles and staircases
baseAngle                       = [53.1 33.2 20.75 12.97 8.1 5.1 3.2 2.0 1.2 0.8 0.5 0.1];                              % Base orientation angle off of vertical
trialStructure                  = [zeros(bonusTrials,1);ones(n_trials,1); ones(n_trials,1)*2; ones(n_trials,1)*3];      % Add "catch" trials to trial structure
trialStructure                  = trialStructure(randperm(length(trialStructure)));                                     % Shuffle staircase presentation order

stairStep1                      = 1;                    % Starting difficulty step for 4:1 staircase
stairStep2                      = 4;                    % Starting difficulty step for 3:1 staircase
stairStep3                      = 9;                    % Starting difficulty step for 2:1 staircase

% Stim parameters
onset                           = 0.25;                 % Ramp onset - 250ms
offset                          = 0.25;                 % Ramp offset - 250ms
duration                        = 0;                    % Duration of stimulus in MS (set to 0 for gabors with onset/offset)
SF                              = 1;                    % c/deg
stimsize                        = 2.5;                  % stim radius in deg
contrast                        = 100;                  % percent contrast 
spatial_envelope                = 1;                    % 0 = disk, 1 = Gabor patch, 2 = raised cosine envelope
which_envelope                  = 1;                    % No idea what this does

results                         = zeros(size(trialStructure, 1) ,7);
stairCount1                     = 0;
stairCount2                     = 0;
stairCount3                     = 0;

Env.x1                          = zeros(139, 139);
Env.y1                          = zeros(139, 139);
Env.circle1                     = zeros(139, 139);

% Mask parameters
orientation                     = 30;                   % irrelevant as you have coherence of 0
Coh                             = 0;
sizeW                           = 151;                  % Width of the stimulus in pixels. MUST BE ODD
SF_cpp                          = 0.001;                % Mean spatial frequency of images in cycles per pixel.
std_SF_cpp                      = 0.01;                 % Std deviation of spatial frequency in cycles per pixel.
envelopesize                    = floor(sizeW/2);

runningScore                    = 0;
trialCount                      = 1;

% Setups where the stimulus will appear on the screen
h_ecc_orig                      = H_ecc_stim;           % Record stimulus X location
v_ecc_orig                      = V_ecc_stim;           % Record stimulus Y location
V_ecc_stim                      = -V_ecc_stim;          % Invert the Y location

stimulus_radius                 = 60*stimsize/scale_factor; % stimulus size adjusted for scale factor
motion_duration                 = 200;
time_sigmaP                     = motion_duration/1000*frame_rate;
Gaussian_stdev                  = round(stimulus_radius/1.5);
f                               = (SF*scale_factor/60)*2*pi;

[x,y]                           = meshgrid(-envelopesize:envelopesize,-envelopesize:envelopesize);
circle                          = envelopesize^2-(x.^2+y.^2);
circle                          = double(circle>0);
R                               = (sqrt((x).^2 + (y).^2) + eps).*circle;
R                               = R/max(max(R));
cos2D                           = (cos(R*pi)+1)/2;
circleCOS                       = (cos2D.*circle);

% Adjust stim and fixation locations for scale factor
H_ecc_stim                      = H_ecc_stim*60/scale_factor;
H_ecc_fix                       = H_ecc_fix*60/scale_factor;
V_ecc_stim                      = V_ecc_stim*60/scale_factor;
V_ecc_fix                       = V_ecc_fix*60/scale_factor;

% make the spatial envelope (dunno what any of this does)
ii                              = 1;
stimulus_radius                 = round(stimulus_radius);
[x,y]                           = meshgrid(-stimulus_radius(ii):stimulus_radius(ii),-stimulus_radius(ii):stimulus_radius(ii));
bps                             = (stimulus_radius(ii))*2+1;circle=((stimulus_radius(ii))^2-(x.^2+y.^2));
for i=1:bps
    for j =1:bps
        if circle(i,j) < 0; circle(i,j) = 0;
        else
            circle(i,j) = 1;
        end
    end
end
if spatial_envelope == 1
    circle = (exp(-(((x)/(sqrt(2)*Gaussian_stdev(ii)/2)).^2)-((y/(sqrt(2)*Gaussian_stdev(ii)/2)).^2)).*circle);
elseif spatial_envelope == 2
    R = (sqrt(x.^2 + y.^2) + eps).*circle;R = R/max(max(R));
    cos2D = (cos(R*pi)+1)/2;circle = (cos2D.*circle);
end
Env(ii).x1 = x; Env(ii).y1 = y; Env(ii).circle1 = circle;

% make temporal envelope
n_onset                         = onset*frame_rate;
n_offset                        = offset*frame_rate;
onset_env                       = sind((0:n_onset)*90/n_onset);
offset_env                      = sind(fliplr(0:n_offset)*90/n_offset);
tempenv                         = horzcat(onset_env,ones(1,duration*frame_rate),offset_env);

% Prep screen presentation
movie_rect                      = [0,0,bps,bps];
scr_left_middle                 = fix(screen_rect(3)/2)-round(bps/2);
scr_top                         = fix(screen_rect(4)/2)-round(bps/2);
screen_rect_middle              = movie_rect + [scr_left_middle, scr_top, scr_left_middle, scr_top];
screen_patch                    = screen_rect_middle + [H_ecc_stim, V_ecc_stim, H_ecc_stim, V_ecc_stim];
screen_patch3                   = screen_rect_middle + [-50,-50,50,50];

sr_hor                          = round(screen_rect(3)/2);
sr_ver                          = round(screen_rect(4)/2);
fix_hor                         = sr_hor + H_ecc_fix;
fix_ver                         = sr_ver + V_ecc_fix;

mv_length                       = length(onset_env) + length(offset_env) + duration*frame_rate;
movie{1,mv_length}              = [];
movie(1:mv_length)              = {zeros(139, 139)};
movie_play{1,mv_length}         = [];

fix_rect                        = SetRect(0, 0, fix_size*scale_factor, fix_size*scale_factor);
fix_rect                        = CenterRectOnPoint(fix_rect, fix_hor, fix_ver);

% Eyetracking setup
if ET
    pref_eye                        = 1; %Which eye to pick, 0 left, 1 right, 2 both 
    broke_fixation_count            = 0;
    [fix_box, el] = setupEyelink(initials, scale_factor, fix_hor, fix_ver, w, screen_rect, rect);
    fix_box = round(fix_box);
end

% Adjust counter and score position based on stimulus position
counterSide = -.8;

% Open landing page and wait for subject to start the task
Screen('FillRect',w, background);
Screen('DrawText',w,'Aim for the High Score!', resolution(1)/2-stimulus_radius*1.5, resolution(2)/2-stimulus_radius*3, 0);
Screen('DrawText',w,'Confident + Correct = +50 points', resolution(1)/2-stimulus_radius*2, resolution(2)/2-stimulus_radius*2, 0);
Screen('DrawText',w,'Not Confident + Correct = +20 points', resolution(1)/2-stimulus_radius*2, resolution(2)/2-stimulus_radius*1.5, 0);
Screen('DrawText',w,'Not confident + Incorrect = 0 points', resolution(1)/2-stimulus_radius*2, resolution(2)/2-stimulus_radius, 0);
Screen('DrawText',w,'Confident + Incorrect = -50 points', resolution(1)/2-stimulus_radius*2, resolution(2)/2-stimulus_radius*.5, 0);
Screen('DrawText',w,'Response keys: 1 if confident left, 2 if uncertain left, 9 if uncertain right, 0 if confident right', ...
    resolution(1)/2-stimulus_radius*5.625, resolution(2)/2+stimulus_radius, 0);
Screen('DrawText',w,'Press any key to start', resolution(1)/2-stimulus_radius*1.375, resolution(2)/2+stimulus_radius*1.5, 0);
Screen('Flip', w);
[secs, keyCode, deltaSecs] = KbWait(-1);

%% Begin task
while trialCount < size(trialStructure,1)+1
    fixated = 1; % Default fixated flag to allow for program to run without eyetracking, is overwritten when tracking is on

    % Draws fixation spot and presents score info
    presentScore(w,background, fontSize, runningScore, resolution, stimulus_radius, counterSide, trialCount, trialStructure, fix_rect);

    % Select Staircase and current difficulty
    % If statements limit range to stay within staircase bounds
    if trialStructure(trialCount) == 0
        if stairStep1 < 5
            angle_deviationP = baseAngle(1);
        else
            angle_deviationP = baseAngle(stairStep1-3);
        end
    elseif trialStructure(trialCount) == 1
        if stairCount1 == 4
            stairStep1 = stairStep1 + 1;
            stairCount1 = 0;
            if stairStep1 > numel(baseAngle)
                stairStep1 = numel(baseAngle);
            end
        end
        angle_deviationP = baseAngle(stairStep1);
    elseif trialStructure(trialCount) == 2
        if stairCount2 == 3
            stairStep2 = stairStep2 + 1;
            stairCount2 = 0;
            if stairStep2 > numel(baseAngle)
                stairStep2 = numel(baseAngle);
            end
        end
        angle_deviationP = baseAngle(stairStep2);
    elseif trialStructure(trialCount) == 3
        if stairCount3 == 2
            stairStep3 = stairStep3 + 1;
            stairCount3 = 0;
            if stairStep3 > numel(baseAngle)
                stairStep3 = numel(baseAngle);
            end
        end
        angle_deviationP = baseAngle(stairStep3);
    end
    
    % Prep eyetracker for stimulus onset
    if ET
        Eyelink('Message', 'Fixation cross drawn');
        Eyelink('Command', 'draw_cross %d %d 15', fix_hor, fix_ver);
        Eyelink('Command', 'draw_box %d %d %d %d 15', fix_box(1), fix_box(2), fix_box(3), fix_box(4));
    end
    
    % Randomize tilt left or right, assign appropriate response key values
    if CoinFlip(1,.5)
        angle_deviationP = angle_deviationP*-1;
        correctHC = '1!';
        correctLC = '2@';
        incorrectHC = '0)';
        incorrectLC = '9(';
    else
        correctHC = '0)';
        correctLC = '9(';
        incorrectHC = '1!';
        incorrectLC = '2@';
    end
    
    %% Play Sample and Target
    priorityLevel = MaxPriority(w);
    Priority(priorityLevel);

    amplitude = background*contrast/100;

    x = Env.x1;y = Env.y1;circle = Env.circle1;
    bps = max(size(circle));
    [time_gauss] = dujeEnvelope(((time_sigmaP)/frame_rate),frame_rate,which_envelope,amplitude);

    angle = deg2rad(angle_deviationP);
    a = cos(angle)*f; 
    b = sin(angle)*f;
    
    if ET
        [fixated] = startStimulusEyelink(trialCount, size(trialStructure,1), pref_eye, w, fix_box, fix_hor, fix_ver);
    end

    Beeper(1000) % Sound to indicate stimulus onset

    %% For Stimulus
    i = 1;
    while i <= mv_length
        ramp_amp = amplitude * tempenv(i); % adjust stimulus amplitude by onset/offset to allow for fading
        movie{i} = round((sin(a*x+b*y) .* circle * ramp_amp) + background);
        movie_play{i} = Screen('MakeTexture', w, movie{i});
        Screen('FillRect', w, background); % Background must be filled BEFORE playing stimulus video
        Screen('DrawTexture', w, movie_play{i}, movie_rect,screen_patch); % Plays the stimulus video

        Screen('TextSize', w, fontSize);
        Screen('DrawText', w, num2str(runningScore), resolution(1)/2-(size(num2str(runningScore),2)*fontSize*.25), resolution(2)/2-stimulus_radius/2*counterSide, 0);
        Screen('TextSize', w, fontSize/2);
        Screen('DrawText', w, [num2str(trialCount), '/', num2str(size(trialStructure,1))], resolution(1)/2-(fontSize/2), resolution(2)/2-stimulus_radius*1.5*counterSide, 0);
        Screen('FillOval', w, 0, fix_rect);
        Screen('Flip', w);

        i=i+1;
        if i>mv_length
            break
        end
        if ET % Enforces fixation while stimulus is playing       
            Eyelink('Command', 'draw_cross %d %d 15', fix_hor, fix_ver);
            Eyelink('Command', 'draw_box %d %d %d %d 15', fix_box(1), fix_box(2), fix_box(3), fix_box(4));
            [breakFlag] = duringStimulusEyelink(w, fix_box, pref_eye);
            if breakFlag
                fixated = 0;
                break;
            end
        end
    end
    
    presentScore(w, background, fontSize, runningScore, resolution, stimulus_radius, counterSide, trialCount, trialStructure, fix_rect);
    
    if ET && fixated==1
        Eyelink('Message', 'Successfuly Fixated Trial');
    end
    if fixated == 0
        if ET
            Eyelink('StopRecording');
        end
        broke_fixation_count=broke_fixation_count+1;
        if broke_fixation_count > 9
            EyelinkDoTrackerSetup(el);
            broke_fixation_count=0;
        end
    elseif fixated == 1 || ET == 0
        % Draw Mask
        im = uint8(genImagesDT(sizeW, SF_cpp, std_SF_cpp, orientation, Coh,circleCOS) * contrast + background);
        mask_play = Screen('MakeTexture',w,im); 
        Screen('FillRect',w, background);
        Screen('DrawTexture', w, mask_play,movie_rect,screen_patch);
        Screen('TextSize',w, fontSize);
        Screen('DrawText',w, num2str(runningScore), resolution(1)/2-(size(num2str(runningScore),2)*fontSize*.25), resolution(2)/2-stimulus_radius/2*counterSide, 0);
        Screen('TextSize',w, fontSize/2);
        Screen('DrawText',w, [num2str(trialCount), '/', num2str(size(trialStructure,1))], resolution(1)/2-(fontSize/2), resolution(2)/2-stimulus_radius*1.5*counterSide, 0);
        Screen('FillOval',w,0,fix_rect);
        Screen('Flip',w);

        WaitSecs(delayMask) % Mask display time in Seconds

        presentScore(w, background, fontSize, runningScore, resolution, stimulus_radius, counterSide, trialCount, trialStructure, fix_rect);

        %% End Presentation, collect response
        tic
        rs = 0;
        validKey = 0;
        while validKey == 0
            [secs, keyCode, deltaSecs] = KbWait([]);
            % Determine if the correct key was selected
            if keyCode(KbName('1!')) || keyCode(KbName('0)')) || keyCode(KbName('2@')) || keyCode(KbName('9('))
                responseTime = toc;    % Record response time
                validKey = 1;
                if keyCode(KbName(correctHC))
                    rs = 1;
                    confidence = 1;
                    runningScore = runningScore + 50;
                elseif keyCode(KbName(correctLC))
                    rs = 1;
                    confidence = 0;
                    runningScore = runningScore + 20;
                elseif keyCode(KbName(incorrectLC))
                    confidence = 0;
                elseif keyCode(KbName(incorrectHC))
                    confidence = 1;
                    runningScore = runningScore - 50;
                end
            elseif keyCode(KbName('Escape'))
                 EyelinkDoTrackerSetup(el)
            end
        end

        % Play correct/incorrect response sounds
        if rs == 1
            if trialStructure(trialCount) == 1
                stairCount1 = stairCount1 + 1;
            elseif trialStructure(trialCount) == 2
                stairCount2 = stairCount2 + 1;
            elseif trialStructure(trialCount) == 3
                stairCount3 = stairCount3 + 1;
            end
            Beeper(1200)
        else
            if trialStructure(trialCount) == 1
                if stairStep1 > 1
                    stairStep1 = stairStep1 - 1;
                end
                stairCount1 = 0;
            elseif trialStructure(trialCount) == 2
                if stairStep2 > 1
                    stairStep2 = stairStep2 - 1;
                end
                stairCount2 = 0;
            elseif trialStructure(trialCount) == 3
                if stairStep3 > 1
                    stairStep3 = stairStep3 - 1;
                end
                stairCount3 = 0;
            end
            Beeper(800)
        end
        keyCode = zeros(1,256);                             % Clear keycode

        % Record trial parameters to results structure
        results(trialCount,1) = trialCount;                 % Trial Number
        results(trialCount,2) = abs(angle_deviationP);      % Difficulty level
        results(trialCount,3) = responseTime;               % Response time
        results(trialCount,4) = angle_deviationP;           % Record base angle
        results(trialCount,5) = rs;                         % Correct/Incorrect
        results(trialCount,6) = confidence;                 % High/Low
        results(trialCount,7) = trialStructure(trialCount); % Which staircase

        presentScore(w, background, fontSize, runningScore, resolution, stimulus_radius, counterSide, trialCount, trialStructure, fix_rect);

        trialCount = trialCount + 1;                        % Increment trial count
        WaitSecs(ITI)                                       % time between trials
    end
end

%% End program - close screen and save output
Screen('CloseAll');
cd('~/Desktop/Subjects')
folder = exist(initials,'file');
if folder ~= 7
    mkdir(initials),
end
cd(initials); 
time = string(datetime('today'));
save(strcat(initials,'_OriMetaCog_',int2str(h_ecc_orig),'_',int2str(v_ecc_orig),'_',time));
