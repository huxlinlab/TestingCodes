%% Eyelink setup     
function [fix_box, el] = setupEyelink(initials, scale_factor, fix_hor, fix_ver ,w , screen_rect, rect)
dummymode                       = 0;                                     % Flag to use the mouse as false eyetracking, turn off if eyetracking!
fix_box_size                    = 5;                                     % Size of fixation box in degrees, will force a break if fixation leaves this window
edfFile                         = initials;
ET_res                          = [1280 1024];                           % Resolution that the ET screen runs at (should be 4:3, 1280x1024 is native to monitor
broke_fixation_count            = 0;
Cal_ecc                         = 8;                                     % in degrees, how far out do we want the calibration targets? Only if using custom
use_small_cal                   = 1;                                     % Use smaller, inverted targets, may get better fixations

fix_box_size = fix_box_size*60/scale_factor;
fix_box = SetRect(-.5*fix_box_size,-0.5*fix_box_size,.5*fix_box_size,0.5*fix_box_size);
fix_box = CenterRectOnPoint(fix_box,fix_hor,fix_ver);
%Initiates Eyelink defaults with details about the graphical environment

el=EyelinkInitDefaults(w);

if use_small_cal
    el.backgroundcolour = BlackIndex(el.window);
    el.msgfontcolour = WhiteIndex(el.window);
    el.imgtitlecolour = WhiteIndex(el.window);

    el.calibrationtargetcolour = WhiteIndex(el.window);
    el.calibrationtargetsize = 1;
    el.calibrationtargetwidth = 0.5;
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
if size(initials,2) > 8
    i = Eyelink('Openfile',initials(1:8));
else
    i = Eyelink('Openfile',initials(1:size(initials,2)));
end
if i~=0
    fprintf('Cannot create EDF file ''%s'' ', edffilename);

    return;
end

% make sure we're still connected.
if Eyelink('IsConnected')~=1 && ~dummymode
    dummymode=1;
    return;
end

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
Eyelink('command', 'set_idle_mode');
Eyelink('command', 'clear_screen 1');

WaitSecs(0.5);
Beeper;

%Let's get started! camera setup mode, calibration and validation
EyelinkDoTrackerSetup(el);
WaitSecs(0.1)
Beeper;