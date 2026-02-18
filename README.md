# TestingCodes
This repository currently contains three primary testing programs:

CoarseDiscrimination_RDK.m

FineDiscrimination_RDK.m

FineOrientationDiscrimination_MetaCog.m

All programs were built in Matlab using Psychtoolbox and designed for use with Mac computers. Testing performed with Matlab 2023b, Psychtoolbox version 3.0.19, on an Intel-based Mac running 13.7.8 Ventura. Modification may be required for use with other systems.
Parameters such as monitor information, testing location, stimulus size, dot density, contrast, and timing variables should be adjusted per patient and testing rig.

Each of these testing codes may require a number of varying subfunctions to use properly. Most subfunctions are listed for each code, but may additionally require the supplied Eyelink Functions for eyetracking, as well as additional Matlab toolboxes.

A brief description of each testing code follows:
CoarseDiscrimination_RDK - This provides a left vs right discrimination task. The stimulus is a random dot kinetogram. Difficulty is adjusted on a 3:1 staircase by increasing the direction range each individual dot can move by 40deg steps. The base direction of left vs right does not change.

FineDiscrimination_RDK - Utilizes the same dot generation program as CoarseDiscrimination_RDK. Provides an Up vs Down discrimination task. Individual dots always move in the provided direction, and difficulty is adjusted with a 3:1 staircase by altering the angle of motion in semi-log steps from 85deg (near vertical) to 0.1deg (near hortizontal). Pre-cues can be turned on or off to add a Feature-Based attention cue prior to stimulus onset. Base direction (left or right) is randomized per trial, with Up or Down the directional axis of interest.

FineOrientationDiscrimination_MetaCog - Provides an oriented Gabor (raised Cosine) tilted left or right from vertical. Contrast is fixed, and a post-mask is provided after each trial. Difficulty is adjusted on a 3:1 staircase to alter the angle of orientation closer to vertical. Confidence reports can be activated to allow patients to rate their confidence for a given trial response, allowing for assessment of cognitive and meta-cognitive abilities.

All programs are works in progress and will potentially be updated in the future. Please check for updated programs to ensure use of the most reliable and robust scripts. Contact the Huxlin Lab for any additional inquiries or troubleshooting.
