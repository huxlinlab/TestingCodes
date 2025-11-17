%% Start Stimulus  
function [fixated] = startStimulusEyelink(trialCount, n_trials, pref_eye, w, fix_box, fix_hor, fix_ver)
%Send start trial message
Eyelink('Message', 'TRIALID %d', trialCount');
Eyelink('command', 'record_status_message "TRIAL %d/%d"', trialCount, n_trials');
%Clear tracker display and draw box at center
Eyelink('Command', 'set_idle_mode');
Eyelink('Command', 'clear_screen %d', 1)
WaitSecs(0.55)
Eyelink('StartRecording');

Eyelink('Message', 'ITI Finished');

% Wait for fixation
Eyelink('StartRecording')
Eyelink('Message', 'Waiting for Fixation');
quitflag = 0;
while ~quitflag
  %  Eyelink('Command', 'draw_cross %d %d 15', fix_hor, fix_ver);
   % Eyelink('Command', 'draw_box %d %d %d %d 15', fix_box(1), fix_box(2), fix_box(3), fix_box(4));
    if Eyelink( 'NewFloatSampleAvailable') > 0 % in dummy mode use mousecoordinates
        sample = Eyelink('NewestFloatSample');
    else
        [x,y] = GetMouse(w);
        sample.gx=x;
        sample.gy=y;
    end
    if pref_eye<2
        if IsInRect(sample.gx(pref_eye+1),sample.gy(pref_eye+1),fix_box)
            sample.gx
            sample.gy, 
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

fixated = 1;

Eyelink('Message', 'Start Stimulus');
