%% First draft feature matching program 8/11/21

%% Program Setup
% Suppresses warnings native to Psychtoolbox, shuffle trials, other housekeeping
oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
Screen('Preference', 'SkipSyncTests', 1);
KbName('UnifyKeyNames');
rng('shuffle')                                                              % Shuffles trials to make each session unique

n_trials                        = 315;                                      % Currently hardcoded. Allows for 105 matching trials and 210 trials of equally mixed faster/slower. 
                                                                            % Designed around running 70 trials for each of three base speeds

H_ecc_stim                      = -5;                                       % X Location in degrees of target stim
V_ecc_stim                      = 5;                                        % Y Location in degrees of target stim

testingPattern                  = 'mirror';                                 % 'mirror'/'same'/'diagonal' -- Indicate where you want the sample to be placed - Mirrored location, same hemifield, or diagonal

H_ecc_fix                       = 0;                                        
V_ecc_fix                       = 0;

% Monitor information
frame_rate                      = 60; 
resolution                      = [1280 1024]; 
screen_width                    = 34;                                       % cm
viewing_dist                    = 42;                                       % cm
fig_width                       = resolution(1)/4;
fig_height                      = resolution(2)/4;
theta                           = atand((screen_width/2)/viewing_dist);
scale_factor                    = theta*60/(resolution(1)/2);               % Arcmin/pixel
ISI                             = (1/frame_rate)*1000; 
windowRect                      = [];

% Beautification
font                            = 'Arial';
fontSize                        = 20;
fix_size                        = 8;
cue_color                       = [255 255 255];
background                      = 128;

%%% Open Screens to be used in session
screens                         = Screen('Screens');
screenNumber                    = max(screens);
[w, rect]                       = Screen('OpenWindow',screenNumber,0,windowRect,[],2);
screen_rect                     = Screen('Rect',w);

Screen(w, 'BlendFunction', GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('FillRect',w, background);
Screen('Flip', w);
Screen('FillRect',w, background);
Screen('TextSize',w, fontSize);
Screen('TextFont',w, font);

% Set up staircases and base speeds
stair1                          = 1;                                            % Used to track position in the staircase
staircount1                     = 0;
baseAngle                       = [45 135];                                    % Base speeds to be tested. Needs at least two values, but can be duplicates
difficultyRange                 = [45 30 15 10 5 2.5 1];                     % Difficulty staircase, represented as percentages of the base speed

LRS                             = [ones(105,1)*-1;ones(105,1);zeros(105,1)];    % Generate list of faster/slower/same (hardcoded counts
LRS                             = Shuffle(LRS);                                 % Shuffle array to randomize
maxJitter                       = 30;
jitterAdjust                    = 15;

% Dot Settings
stimulus_duration               = 500;                                          % ms
aperture_radius                 = 2.5;                                          % stimulus radius in degrees
dot_density                     = 3.5;                                          % number of dots per square deg of stimulus
initial_dot_size                = 14;                                           % diamteter, arcmin
dot_color                       = 0;                                            % Grayscale units
dot_lifetime                    = 200;                                          % in ms, for direction range only
dot_size                        = floor(initial_dot_size/scale_factor);         % size of each dot adjusted for scale factor

% Setups where the stimulus will appear on the screen
h_ecc_orig                      = H_ecc_stim;                                   % Record stimulus X location
v_ecc_orig                      = V_ecc_stim;                                   % Record stimulus Y location
V_ecc_stim                      = -V_ecc_stim;                                  % Invert the Y location to fit the screen settings
ndots                           = round(aperture_radius^2*pi*dot_density);      % Total number of dots in the stimulus
stimulus_radius                 = 60*aperture_radius/scale_factor;              % stimulus size adjusted for scale factor
lifetime                        = dot_lifetime/(1000/frame_rate);
bps                             = (stimulus_radius)*2+1;                        % no idea what this is
dot_speed                       = 10;
Velocity                        = dot_speed*60/scale_factor/(1000/ISI);         

% Adjust stim and fixation locations for scale factor
H_ecc_stim                      = H_ecc_stim*60/scale_factor;
H_ecc_fix                       = H_ecc_fix*60/scale_factor;
V_ecc_stim                      = V_ecc_stim*60/scale_factor;
V_ecc_fix                       = V_ecc_fix*60/scale_factor;

% Pre allocate various variables for dot generation
mv_length                       = round(stimulus_duration/(1000/frame_rate));
age{1,mv_length}                = [];
age(1:mv_length)                = {zeros(1, ndots)};
positions{1,mv_length}          = [];
postions(1:mv_length)           = {zeros(1, ndots)};
positions2{1,mv_length}         = [];
postions2(1:mv_length)          = {zeros(1, ndots)};
results                         = zeros(n_trials,5);

% Prep screen presentation
movie_rect                      = [0,0,bps,bps];
scr_left_middle                 = fix(screen_rect(3)/2)-round(bps/2);
scr_top                         = fix(screen_rect(4)/2)-round(bps/2);
screen_rect_middle              = movie_rect + [scr_left_middle, scr_top, scr_left_middle, scr_top];
screen_patch                    = screen_rect_middle + [H_ecc_stim, V_ecc_stim, H_ecc_stim, V_ecc_stim];

sr_hor                          = round(screen_rect(3)/2);
sr_ver                          = round(screen_rect(4)/2);
fix_hor                         = sr_hor + H_ecc_fix;
fix_ver                         = sr_ver + V_ecc_fix;
stim_hor                        = sr_hor + H_ecc_stim;
stim_ver                        = sr_ver + V_ecc_stim;

switch testingPattern
    case {'mirror'}
        stim_hor2                       = sr_hor - H_ecc_stim;
        stim_ver2                       = sr_ver + V_ecc_stim;
    case {'same'}
        stim_hor2                       = sr_hor + H_ecc_stim;
        stim_ver2                       = sr_ver - V_ecc_stim;
    case {'diagonal'}
        stim_hor2                       = sr_hor - H_ecc_stim;
        stim_ver2                       = sr_ver - V_ecc_stim;
end

fix_rect                        = SetRect(0, 0, fix_size*scale_factor, fix_size*scale_factor);
fix_rect                        = CenterRectOnPoint(fix_rect, fix_hor, fix_ver);

for trialCount = 1:n_trials    
    validKey = 0;
    
    % Adjust difficulty based on correct trial count
    if staircount1 == 2
        stair1 = stair1 + 1;
        if stair1 > length(difficultyRange)
            stair1 = length(difficultyRange);
        end
        staircount1 = 0;
    end
    
    % Selects current difficulty as a percent and modifies it by LRS
    difficultyMagnitude = 1 + (difficultyRange(stair1)/100) * LRS(trialCount);
        
    % Generate velocity values for the two stimuli
    Velocity = dot_speed*60/scale_factor/(1000/ISI);    % Sample speed
    
    % Determine what keystroke is "correct" based on LRS
    if LRS(trialCount) == 1 
        correct = 'RightArrow';
    elseif LRS(trialCount) == -1
        correct = 'LeftArrow';
    elseif LRS(trialCount) == 0
        correct = 'UpArrow';
    end
    
    while ~validKey
                
        bb = GetSecs;
        % Generate Sample Dots
        angle = randsample(baseAngle,1)+(randi(maxJitter)-jitterAdjust);
        [positions] = GenerateDots(ndots, bps, stimulus_radius, angle, Velocity, lifetime, mv_length);
        
        % Generate Target Dots
        angle2 = angle+difficultyRange(stair1)*LRS(trialCount);
        [positions2] = GenerateDots(ndots, bps, stimulus_radius, angle2, Velocity, lifetime, mv_length);
        
        %% Play Sample and Target
        % Finish the ITI
        WaitSecs(.5-(GetSecs-bb));
        priorityLevel = MaxPriority(w);
        Priority(priorityLevel);

        Beeper(1000)
        WaitSecs(0.05);

        % Play the movies
        i = 1;
        while i <= mv_length
            Screen('FillRect',w, background);
            Screen(w,'DrawDots',transpose(positions{i}), dot_size, dot_color, [stim_hor2 stim_ver2],2);
            Screen(w,'DrawDots',transpose(positions2{i}), dot_size, dot_color, [stim_hor stim_ver],2);
            Screen('FillOval',w,0,fix_rect);
            Screen('Flip',w);
            i=i+1;
%             [secs, keyCode, deltaSecs] = KbWait;
%             WaitSecs(0.05);
            if i>mv_length
                break
            end
        end
        
        tic
        Screen('Flip',w);
        rs = 0;
        [secs, keyCode, deltaSecs] = KbWait;
        responseTime = toc;
        
        % Determine if the correct key was selected
        if keyCode(KbName('LeftArrow'))
            validKey = 1;
            if keyCode(KbName(correct))
                rs = 1;
            end
        elseif keyCode(KbName('UpArrow'))
            validKey = 1;
            if keyCode(KbName(correct))
                rs = 1;
                staircount1 = staircount1 - 1;
            end
        elseif keyCode(KbName('RightArrow'))
            validKey = 1;
            if keyCode(KbName(correct))
                rs = 1;
            end
        end
        results(trialCount,2) = difficultyRange(stair1);
        
        % Play correct/incorrect response sounds, adjust staircase based on
        % performance
        if rs == 1
            Beeper(1200)
            staircount1 = staircount1 + 1;
        else
            Beeper(800)
            stair1 = stair1 - 1;
            staircount1 = 0;
            if stair1 == 0
                stair1 = 1;
            end
        end
    end
    % Record trial parameters to results structure
    results(trialCount,1) = trialCount;
    results(trialCount,3) = LRS(trialCount);
    results(trialCount,4) = rs;
    results(trialCount,5) = responseTime;
    Screen('FillOval',w,0,fix_rect);
    Screen('Flip',w);
end
Screen('CloseAll');
