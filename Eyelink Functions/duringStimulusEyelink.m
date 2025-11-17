%% During stimulus            
function [breakFlag] = duringStimulusEyelink(w, fix_box, pref_eye)
breakFlag = 0;
if Eyelink( 'NewFloatSampleAvailable') > 0 % in dummy mode use mousecoordinates
    sample = Eyelink('NewestFloatSample');
else
    [x,y] = GetMouse(w);
    sample.gx=x;
    sample.gy=y;
end
if pref_eye<2
    if length(sample.gx) ==1
        if IsInRect(sample.gx, sample.gy,fix_box) == 0
            Eyelink('Message', 'Fixation Broken');
            Beeper;
            Beeper;
            Beeper;
            breakFlag = 1;
        end
    else
        if IsInRect(sample.gx(pref_eye+1),sample.gy(pref_eye+1),fix_box) == 0
            Eyelink('Message', 'Fixation Broken');
            Beeper;
            Beeper;
            Beeper;
            breakFlag = 1;
        end
    end
else
    if IsInRect(sample.gx(1),sample.gy(1),fix_box) == 0 || IsInRect(sample.gx(2),sample.gy(2),fix_box) == 0
        Eyelink('Message', 'Fixation Broken');
        Beeper;
        Beeper;
        Beeper;
        breakFlag = 1;
    end
end
