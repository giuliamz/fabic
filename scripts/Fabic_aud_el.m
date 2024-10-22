clear all
clc

% Request parameters via dialog box

prompt = {'Subject ID [MA#]', 'Subject nr', 'Task Name [Pract; Exp]', 'Run number'};
title = 'Input';
lines = [1 30];
default = {'TEST', '1', 'Exp', '1'};
params = inputdlg(prompt, title, lines, default);
if ~isempty(params) % ok is pressed
    [S.subj_ID, S.subj_n, S.taskName, S.run] = params{[1 2 3 4]};   
else % cancel is pressed
    fprintf('Run has been aborted...\n')
    return
end

%Ttab is an excell file we built with info about trial number, trial condition, response condition, stimuli position 
%Ctab is a cell version of TTab
Tablename = append('trial_list\', 'trial_list_r', S.run, '_s', S.subj_n, '.xlsx');
Rtab=readtable(Tablename);
Rtab=table2cell(Rtab);

% Create folder data and subfolder with subj name
savePath = fullfile(pwd, 'data', S.subj_ID);
if ~isfolder(savePath) % create new folder if needed
    mkdir(savePath)
end
%% stimuli parameters    
% Timing of stimuli
S.timing.fixation = 3000;                   %fixation duration in milliseconds

S.fixpres = 500; %fixation time pre stimulus onset
S.restime = 1500; %response time

S.stim.eccentricity = 9;                    %degrees of visual angle

S.stim.aud_fsample = 48000;                 %sampling frequency in Hz


S.TrialList = Rtab;
S.TrialTot = size(Rtab,1);
S.order.comm = cell2mat(Rtab(:,1));
S.order.response = cell2mat(Rtab(:,2));
S.order.aud_loc = cell2mat(Rtab(:,3));
S.order.vis_loc = cell2mat(Rtab(:,4));
S.vid.name = Rtab (:,5);
S.aud.name = Rtab(:,6);
S.intro = 1000;


S.mon.dist    = 48;            %in centimeters
S.mon.width   = 53;
S.mon.height  = 30;
S.mon.res = [1920 1080];
S.mon.vis_hor_offset    = 0;   %in pixels 
S.mon.vis_ver_offset    = 0;   %in pixels
%S.visangle = 275; % 9 degrees in pixels
%S.visangle = 559.66; % 9 degrees dist 1000 mm
S.visangle = 447.73; %9 degress dist 800 mm


%define stimuli position
c_rect = [0 0 960 540]; % rectangle for scale

KbName('UnifyKeyNames');
%define keys for input
P.key1=KbName('LeftArrow');
P.key2=KbName('DownArrow');
P.key3=KbName('RightArrow');
P.quitKey=KbName('ESCAPE');
%P.quitkey=KbName('ESCAPE');



   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Prepare Audio Device %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % Open sound driver for high timing precision, stereo playback with a certain given sampling frequency: 
    % This has to be called with input argument 1 in order to achieve really low latency!
    InitializePsychSound(1);
    % Set the audio device number and driver's mode (1 = playback only (no
    % audio capture)) and then open the device.
    % Find the correct audiodriver
    audioDevices = PsychPortAudio('GetDevices');
    audioDeviceNames = {audioDevices.DeviceName}';     %The names of all the sound drivers (including the virtual drivers)
    audioDeviceIndices = {audioDevices.DeviceIndex}';  %The indices of all sound drivers
    P.audioDeviceChannels = 2;                           %The number of available channels
    latBias = 0;
    % reqLatencyClass level 2 means: Take full control over the audio device, even if this causes other sound applications to fail or shutdown.
    reqLatencyClass = 2;
    %stimDeviceNr = audioDeviceIndices{strcmp('Speakers (Realtek High Definition Audio)',audioDeviceNames)}; %windows ufficio
    %stimDeviceNr = audioDeviceIndices{strcmp('Built-in Output',audioDeviceNames)}; %mac ufficio
    stimDeviceNr = audioDeviceIndices{strcmp('Speakers (Creative SB X-Fi)',audioDeviceNames)};
    %stimDeviceNr = audioDeviceIndices{strcmp('Speakers (High Definition Audio Device)',audioDeviceNames)};
    % Open the device
    P.pahandle_stim = PsychPortAudio('Open', stimDeviceNr, [], reqLatencyClass, 48000, P.audioDeviceChannels);
    % Set the visuo-auditory delay compensation
    PsychPortAudio('LatencyBias', P.pahandle_stim, latBias);
    % Set the volume of the master-device fully open 
    PsychPortAudio('Volume', P.pahandle_stim, 2);

     
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize Screen %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%

    % Get the list of Screens and choose the one with the highest screen number. Choosing the display with the highest dislay number is a best guess about where you want the stimulus displayed.
    % Screen 0 is, by definition, the display with the menu bar. Often when two monitors are connected the one without the menu bar is used as the stimulus display.  
    P.screenNumber = max(Screen('Screens'));
    PsychImaging('PrepareConfiguration');

    %PsychDebugWindowConfiguration;

    % Colours
    P.white = WhiteIndex(P.screenNumber);
    P.black = BlackIndex(P.screenNumber);
    %P.grey = round((P.white + P.black) * (1 - S.stim.vis_contr));

    % Open window. 
    [P.win,P.winRect]=PsychImaging('OpenWindow',P.screenNumber,P.black);                
    % Hide the cursor
    HideCursor;
    % Switch to realtime-priority to reduce timing jitter and interruptions caused by other applications and the operating system itself:
    Priority(MaxPriority(P.win));
    % Prevent keyboard presses to be executed within MATLAB (we don't want to edit the scripts!). This cannot be used if KbQueue commands are required.            
    %ListenChar(2);
    % Enable antialiasing blending function (so we can use line smoothing in the DrawLines command)
    Screen('BlendFunction', P.win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    % Text Specifics
    Screen('TextFont',P.win,'Arial');
    Screen('TextSize',P.win, 70);
    Screen('TextStyle',P.win,0);                %0=normal, 1=bold, 2=italic, 4=underline
    Screen('TextColor',P.win,P.white);

    % Measure monitor refresh interval.
    % Flip three times (bug). Otherwise screen will start flickering
    Screen('Flip',P.win,[],P.black);    Screen('Flip',P.win,[],P.black);     Screen('Flip',P.win,[],P.black);
    % This will trigger a calibration loop of minimum 100 valid samples and return the estimated inter-flip-interval in 'P.ifi': interflip interval (in seconds).
    % We require an accuracy of 0.5 ms == 0.0005 secs. If this level of accuracy can't be reached, we time out after 5 seconds.  
    [P.ifi,~,~] = Screen('GetFlipInterval',P.win,100,0.0005,5);
    % Calculate the faxation time in flips      
     P.fixpres = round(((S.fixpres)/1000)/P.ifi);      
     % Calculate the Post-Stimulus faxation and response time in flips [The time between stimulus and response cue] 
     P.restime = round(((S.restime)/1000)/P.ifi);
  
    %P.attcue_in_flips = round((S.timing.attcue/1000)/P.ifi);
    P.fixation_in_flips = round((S.timing.fixation/1000)/P.ifi);
    P.intro_in_flips = round((S.intro/1000)/P.ifi);

     % Define window center, width and height (in pixels).
     [P.winCenter_x,P.winCenter_y] = RectCenter(P.winRect);                   
     % Calculate the size of one pixel (in meters).
     minXpix = P.winRect(1)+1; 
     minYpix = P.winRect(2)+1;
     maxXpix = P.winRect(3)-1;
     maxYpix = P.winRect(4)-1;
     winWidth = maxXpix-minXpix+1;                                            %width of screen in nr of pixels
     winHeight = maxYpix-minYpix+1;                                           %height of screen in nr of pixels   
     P.width_of_1pix = S.mon.width/winWidth;                        %in meters
     P.height_of_1pix = S.mon.height/winHeight;                     %in meters
     % Express the distance to the screen in number of pixels (separate for width and height) 
     P.dist2Screen_nPixWidth = S.mon.dist/P.width_of_1pix;
     P.dist2Screen_nPixHeight = S.mon.dist/P.height_of_1pix;
     % We can now calculate the number of pixels from the screen_centre by: Nr_pixels = round(tand(visual_angle)*P.dist2Screen_nPixWidth);    (or use P.dist2Screen_nPixHeight instead)
     % Do we need an offset of the centre of the screen?
     P.winCenter_x = P.winCenter_x + S.mon.vis_hor_offset;
     P.winCenter_y = P.winCenter_y + S.mon.vis_ver_offset;
     clear minXpix minYpix maxXpix maxYpix winWidth winHeight
     
for iTrial = 1:S.TrialTot
    % auditory locations
    if S.order.aud_loc(iTrial) == 1 %left position
        S.trials.aud_locations(iTrial) = 171;
    elseif S.order.aud_loc(iTrial) == 2 %central position
        S.trials.aud_locations(iTrial) = 180;
    elseif S.order.aud_loc(iTrial) == 3 %right position
        S.trials.aud_locations(iTrial) = 189;
    end
    if S.order.vis_loc(iTrial) == 1
        x_pos(iTrial)=P.winCenter_x-S.visangle;        
    elseif S.order.vis_loc(iTrial) == 2
        x_pos(iTrial)=P.winCenter_x;
    elseif S.order.vis_loc(iTrial) == 3
        x_pos(iTrial)=P.winCenter_x+S.visangle;
    end  
        %x_pos(iTrial)=P.winCenter_x;
    
end

    %%%%%%%%%%%%%%%%%%%%%
    %%% Prepare Cross %%%
    %%%%%%%%%%%%%%%%%%%%%
    % Drawing settings
    S.draw.fixx_size      = 0.7;  %horizontal and vertical visual angle (in degrees) for the fixation cross size (diameter, not radius)
    S.draw.fixx_LineWidth = 5;    %in pixels (integers only) for the fixation cross line-width

    fixxWidth_pix = round(tand(S.draw.fixx_size/2)*P.dist2Screen_nPixWidth);            %Half the size in pixels
    fixxHeight_pix = round(tand(S.draw.fixx_size/2)*P.dist2Screen_nPixHeight);
    % Coords for DrawLines always in order: 1st row = x_coords, 2nd row = y_coords: 1st column pair is start_coord, 2nd column is stop_coord).
    P.fixxCoords = [ [-fixxWidth_pix fixxWidth_pix; 0 0] [0 0; -fixxHeight_pix fixxHeight_pix] ];
    P.cross_colour = [128 128 128];
    P.text_colour = [128 128 128];


%% video frames spatialization
nTrial=0;
numTrials = S.TrialTot; %el trial num
z=0; %early responses counter

if strcmp(S.taskName, 'Pract')
    text=['Utilizza le frecce a sx, in basso e a dx per indicare \n' ...
    'la posizione dell audio se compare A, \n' ...
    'indica la posizione del video se compare V. \n' ...
    'Premi un tasto quando sei pront*'];
    DrawFormattedText(P.win, text, 'center', 'center', P.text_colour);
    Screen('Flip', P.win, []);
                        
    KbWait
elseif strcmp(S.taskName, 'Exp')
    text='Premi un tasto quando sei pront*';
    DrawFormattedText(P.win, text, 'center', 'center', P.text_colour);
    Screen('Flip', P.win, []);
                        
    KbWait
end

for iTrial = 1:S.TrialTot
    if iTrial == 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% EyeLink Calibration %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    dummymode=0;       % set to 1 to initialize in dummymode 
    fixWinSize = 100;
    infix=0;
    dotSize = 10;
    graceTime = GetSecs + 200/1000;
    
    fixationWindow = [-fixWinSize -fixWinSize fixWinSize fixWinSize];
    fixationWindow = CenterRect(fixationWindow, P.winRect);

    el=EyelinkInitDefaults(P.win);
    ListenChar(2);
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit(dummymode, 1)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end
    
    % for lower resolutions you might have to play around with these values
    % a little. If you would like to draw larger targets on lower res
    % settings please edit PsychEyelinkDispatchCallback.m and see comments
    % in the EyelinkDrawCalibrationTarget function
    el.calibrationtargetsize= 1;
    el.calibrationtargetwidth=0.5;
    % call this function for changes to the calibration structure to take
    % affect
    EyelinkUpdateDefaults(el);
    
    % make sure we're still connected.
    if Eyelink('IsConnected')~=1 && ~dummymode
        cleanup;
        return;
    end;
        
    % open file to record data to
    Filename = append(S.subj_ID, '_', S.run, '.edf');
    edfFile=Filename;
    Eyelink('Openfile', edfFile);
 
    
%     % start recording eye position
%     Eyelink('StartRecording');
%     % record a few samples before we actually start displaying
%     WaitSecs(0.1);
%     % mark zero-plot time in data file
%     Eyelink('Message', 'SYNCTIME');
%     stopkey=KbName('space');
%     eye_used = -1;

    Screen('FillRect', el.window, 0);
    Screen('TextFont', el.window, el.msgfont);
    Screen('TextSize', el.window, el.msgfontsize);
    [width, height]=Screen('WindowSize', el.window);
    message='Press space to stop.';
    Screen('DrawText', el.window, message, 200, height-el.msgfontsize-20, 128);
    Screen('Flip',  el.window, [], 1);
    
    % SET UP TRACKER CONFIGURATION
    % Setting the proper recording resolution, proper calibration type,
    % as well as the data file content;
    Eyelink('command', 'add_file_preamble_text ''Recorded by EyelinkToolbox demo-experiment''');
   
    % This command is crucial to map the gaze positions from the tracker to
    % screen pixel positions to determine fixation
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
    
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);
    % set calibration type.
    Eyelink('command', 'calibration_type = HV5');
    Eyelink('command', 'generate_default_targets = YES');
    % set parser (conservative saccade thresholds)
    Eyelink('command', 'saccade_velocity_threshold = 35');
    Eyelink('command', 'saccade_acceleration_threshold = 9500');
    % set EDF file contents
        % 5.1 retrieve tracker version and tracker software version
    [v,vs] = Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    vsn = regexp(vs,'\d','match');
    
        if v ==3 && str2double(vsn{1}) == 4 % if EL 1000 and tracker version 4.xx
        
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
    
    % calibration/drift correction target
    Eyelink('command', 'button_function 5 "accept_target_fixation"');
   
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);

    % do a final check of calibration using driftcorrection
    EyelinkDoDriftCorrection(el);
    
        %fixation pretrial
        Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, P.cross_colour, [P.winCenter_x,P.winCenter_y],0);
        [~,R.Onset] = Screen('Flip',P.win);

        Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, P.cross_colour, [P.winCenter_x,P.winCenter_y],0);
        Screen('Flip',P.win, R.Onset + P.intro_in_flips*P.ifi - P.ifi/2);
    end 
    
        Eyelink('Message', 'TRIALID %d', iTrial);
        % This supplies the title at the bottom of the eyetracker display
        Eyelink('command', 'record_status_message "TRIAL %d/%d"', iTrial,numTrials);
        Eyelink('Command', 'set_idle_mode');
        % clear tracker display and draw box at center
        Eyelink('Command', 'clear_screen %d', 0);
        % draw fixation and fixation window shapes on host PC
        Eyelink('command', 'draw_cross %d %d 15', width/2,height/2);
        Eyelink('command', 'draw_box %d %d %d %d 15', fixationWindow(1), fixationWindow(2), fixationWindow(3), fixationWindow(4));
   % Write !V IAREA message to EDF file: creates interest area around image in DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Interest Area Commands
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, fixationWindow(1), fixationWindow(2), fixationWindow(3), fixationWindow(4),'fixationArea');
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 1, 0, 0, width, height,'ScreenArea');
        Eyelink('Command', 'set_idle_mode');
        WaitSecs(0.05);
        Eyelink('StartRecording');
        % record a few samples before we actually start displaying
        WaitSecs(0.1);
        % mark zero-plot time in data file
        Eyelink('Message', 'SYNCTIME');
        stopkey=KbName('space');
        %eye_used = -1;
        eye_used = Eyelink('EyeAvailable'); % get eye that's tracked  
        % returns 0 (LEFT_EYE), 1 (RIGHT_EYE) or 2 (BINOCULAR) depending on what data is
        if eye_used == 2
            eye_used = 1; % use the right_eye data
        end
        % record a few samples before we actually start displaying
        % otherwise you may lose a few msec of data
        WaitSecs(0.1);
    
    R.early(iTrial)={NaN};

        
    % Draw background onto full screen and get first flip timestamp
    % before images/movie presentation %
    Screen('FillRect', P.win, P.black);
    [~, ~, lastEventTime] = Screen('Flip', P.win);
    R.Onsettrial = lastEventTime
   
    % Draw the fixation cross
    Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, P.cross_colour, [P.winCenter_x,P.winCenter_y],0);
    [~,R.Onsetfix] = Screen('Flip',P.win);
    t0=GetSecs;
    Mfixation = true;
    ti=GetSecs;
    while (ti < R.Onsetfix + S.fixpres/1000) %|| Mfixation  % loop till error or space bar is pressed
        ti=GetSecs;
        % Check recording status, stop display if error
        error=Eyelink('CheckRecording');
        if(error~=0)
            break;
        end
%         % check for keyboard press
%         [keyIsDown, secs, keyCode] = KbCheck;
%         % if spacebar was pressed stop display
%         if keyCode(stopkey)
%             break;
%         end
        % check for presence of a new sample update
        if Eyelink( 'NewFloatSampleAvailable') > 0
            % get the sample in the form of an event structure
            evt = Eyelink( 'NewestFloatSample');
            if eye_used ~= -1 % do we know which eye to use yet?
                % if we do, get current gaze position from sample
                x = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                y = evt.gy(eye_used+1);
                % do we have valid data and is the pupil visible?
                if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                    % if data is valid, draw green cross
                    if 1==IsInRect(x,y, fixationWindow )
                        Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, [128 128 128], [P.winCenter_x,P.winCenter_y],0);
                        Screen('Flip',P.win);
                         Mfixation = false;
                    else
                        Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, [128 128 128], [P.winCenter_x,P.winCenter_y],0);
                        Screen('Flip',P.win);
                         Mfixation = true;
                    end
                    
                else
                    % if data is invalid (e.g. during a blink), draw red
                    % cross
                    Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, [128 128 128], [P.winCenter_x,P.winCenter_y],0);
                    Screen('Flip',P.win);
                end
            else % if we don't, first find eye that's being tracked
                eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
                if eye_used == el.BINOCULAR; % if both eyes are tracked
                    eye_used = el.LEFT_EYE; % use left eye
                end
            end
        end % if sample available
        ti=GetSecs;
    end % main loop
    

        
    Filenamevid=char(S.vid.name(iTrial));
    video=struct;
    video.name=Filenamevid;
    video.frame.format= '.jpg';

       
%     if S.order.comm(iTrial) == 1
%     %count the number of images in a folder
%     folder=fullfile('C:\Users\HPZ22018292\Documents\MATLAB\Fabic\vid_stimuli_res\', Filenamevid);
%     a = dir(fullfile(folder, '*.jpg'));
% 
%     video.path=append("\vid_stimuli_res\", Filenamevid);
%     video.frame.rate= 30;
%     video.frame.number= length(a);
%     %createFramesTextureStructure
%     video.frames = struct;
%     
%     % fill it with the images the PTB way
%         for s = 1:video.frame.number
%                 video.frames(s).readImage = imread(fullfile(cd, ...
%                                                             video.path, ...
%                                                             strcat(video.name, ...
%                                                                    num2str(s), ...
%                                                                    video.frame.format)));
%                 video.frames(s).videoFrameTexture = Screen('MakeTexture', ...
%                                                            P.win, ...
%                                                            video.frames(s).readImage);
%          % calculate the duration of a single frame based on the provided desired frame rate %
%         video.frame.duration = 1 / video.frame.rate;
%         % Correct the stimulus duration (to ifi's)
%         S.timing.vis_duration = video.frame.duration * 1000 * video.frame.number;
%         P.stim_dur_in_flips(iTrial) = round((S.timing.vis_duration/1000)/P.ifi);  
%         end 
% 
%     elseif S.order.comm(iTrial) == 0
%         folder=fullfile('C:\Users\HPZ22018292\Documents\MATLAB\Fabic\vid_stimuli_res\voc');
%         video.path=append("\vid_stimuli_res\voc");
%     
% 
%         video.voc.readImage = imread(fullfile(cd, ...
%                                             video.path, ...
%                                             strcat(video.name, ...
%                                                    video.frame.format)));
%         video.voc.videoFrameTexture = Screen('MakeTexture', ...
%                                                        P.win, ...
%                                                        video.voc.readImage);
%         P.stim_dur_in_flips(iTrial) = round((500/1000)/P.ifi);
%     end    
%    
%     %playVideoFrames 
%     
%            
%         c_rect = CenterRectOnPoint(c_rect, x_pos(iTrial), 525); 

        %load audio stim
        Filenameaud=append('aud_hrtf\', S.aud.name(iTrial),'.wav');
        y=audioread(char(Filenameaud));
        PsychPortAudio('FillBuffer',P.pahandle_stim,y');
        
            
        % Predict V onset to schedule A onset
        R.OnsetAV = PredictVisualOnsetForTime(P.win, R.Onsetfix + P.fixpres*P.ifi - P.ifi/2);
        
        % Write message to EDF file to mark the start time of stimulus presentation.
        Eyelink('Message', 'STIM_ONSET');  
        
        % We schedule the sound to begin after the ISI. This is to provide the system enough time to prepare the sound device.
        PsychPortAudio('Start', P.pahandle_stim, [], R.OnsetAV);  

%         % frames presentation loop
%         if S.order.comm(iTrial) == 1
%             for f = 1:video.frame.number
%     
%                    if f==1
%                        
%                        Screen('DrawTexture', ...
%                                 P.win, ...
%                                 video.frames(f).videoFrameTexture, [], c_rect , 0);
%                        if x_pos(iTrial) ~= P.winCenter_x
%                             Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, P.cross_colour, [P.winCenter_x,P.winCenter_y],0);
%                        end
%                           [~, ~, lastEventTime] = Screen('Flip', ...
%                                                          P.win, ...
%                                                          R.OnsetAV);
%                           R.video_onset(iTrial)=lastEventTime;
%                     elseif f~=1
%                              
%                              Screen('DrawTexture', ...
%                                        P.win, ...
%                                        video.frames(f).videoFrameTexture, [], c_rect , 0);
%                        if x_pos(iTrial) ~= P.winCenter_x
%                             Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, P.cross_colour, [P.winCenter_x,P.winCenter_y],0);
%                        end
%                        
%                          [~, ~, lastEventTime] = Screen('Flip', ...
%                                                         P.win, ...
%                                                         lastEventTime+0.031);
%                          [KBpressed, KbTime, keyCode] = KbCheck(-1);
% 
%                    
%                         if KBpressed %&& KbTime <= lastEventTime && ~keyCode(P.quitKey)
%                             R.early(iTrial)={KbName(keyCode==1)};                        
%                         end
%                     end 
%         
%             end
%         elseif S.order.comm(iTrial) == 0  
% 
%             Screen('DrawTexture', ...
%                         P.win, ...
%                         video.voc.videoFrameTexture, [], c_rect , 0);
%                        if x_pos(iTrial) ~= P.winCenter_x
%                             Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, P.cross_colour, [P.winCenter_x,P.winCenter_y],0);
%                        end
%                        [~, ~, lastEventTime] = Screen('Flip', ...
%                                              P.win, ...
%                                              R.OnsetAV-0.016);
%               R.video_onset(iTrial)=lastEventTime;
%               WaitSecs(0.5);
%               [KBpressed, KbTime, keyCode] = KbCheck(-1);
%  
%                 if KBpressed && KbTime <= (R.video_onset(iTrial)+500) && ~keyCode(P.quitKey)
%                     R.early(iTrial)={KbName(keyCode==1)};                        
%                 end
% 
%         end

        % Write message to EDF file to mark time when blank screen is presented
        Eyelink('Message', 'BLANK_SCREEN');
        
        Eyelink('StopRecording');
        % Write TRIAL_RESULT message to EDF file: marks the end of a trial for DataViewer
        % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Defining the Start and End of a Trial
        Eyelink('Message', 'TRIAL_RESULT 0');

        %count early responses
        if strcmp(R.early(iTrial), 'DownArrow') || strcmp(R.early(iTrial), 'LeftArrow') || strcmp(R.early(iTrial), 'RightArrow')
            z=z+1;
        end
        
        status =PsychPortAudio('GetStatus', P.pahandle_stim);
        R.audio_onset(iTrial)=status.StartTime;
        R.audio_offset(iTrial)=status.EstimatedStopTime;
        R.audio_duration(iTrial)=(R.audio_offset(iTrial)-R.audio_onset(iTrial))*1000; %audio duration in ms

%             %remove visual data and start PSI (time between stimulus and
%             %response cue)
%             if S.order.comm(iTrial) == 1
%                 Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, P.cross_colour, [P.winCenter_x,P.winCenter_y],0);
%                 [~,R.OnsetPSI] = Screen('Flip',P.win); % R.OnsetAV + P.stim_dur_in_flips*P.ifi - P.ifi/2);
%             elseif S.order.comm(iTrial) == 0
%                 Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, P.cross_colour, [P.winCenter_x,P.winCenter_y],0);
%                 [~,R.OnsetPSI] = Screen('Flip',P.win,R.OnsetAV + P.stim_dur_in_flips(iTrial)*P.ifi - P.ifi/2); 
%             end
            
            %Stop sound playback
            PsychPortAudio('Stop', P.pahandle_stim, 1); % waitForEndOfPlayback = 1 --> but it should have ended by now

            %report cue 

            if S.order.response(iTrial) == 1 %aud requested
                correct_key = [P.key1, P.key2, P.key3];
                cue = 'A';
                elseif S.order.response(iTrial) == 0 %vis requested
                correct_key = [P.key1, P.key2, P.key3];
                cue = 'A';
            end

            %% Ask for a response
            R.timeout = 1500; %response time
            QuestionAnsweredBool = 1;
            
            Screen('TextSize', P.win, 70);

 %draw response cue
%             if S.order.comm(iTrial) == 1
                DrawFormattedText(P.win, cue, P.winCenter_x-fixxWidth_pix, P.winCenter_y+fixxHeight_pix, P.cross_colour);
                Screen('DrawingFinished', P.win);
                [~,R.Onsetres] = Screen('Flip',P.win);
%             elseif S.order.comm(iTrial) == 0
%                 DrawFormattedText(P.win, cue, P.winCenter_x-fixxWidth_pix, P.winCenter_y+fixxHeight_pix, P.cross_colour);
%                 Screen('DrawingFinished', P.win);
%                 [~,R.Onsetres] = Screen('Flip',P.win,R.OnsetAV + P.stim_dur_in_flips(iTrial)*P.ifi - P.ifi/2);
%             end            
          
        while QuestionAnsweredBool == 1   
                % Get keypress
                [KBpressed, KbTime, keyCode] = KbCheck(-1);

                if KBpressed && (KbTime-R.Onsetres)*1000 <= R.timeout && ~keyCode(P.quitKey)

                     % Check key press...
                    if  sum((keyCode==1))> 1 || ~any(correct_key == find(keyCode==1))    % incorrect key...
                        R.incorrectButtonClick = 1;
                        
                    elseif sum((keyCode==1))== 1 && any(correct_key == find(keyCode==1)) % correct key...
                        
                        RT = (KbTime-R.Onsetres)*1000; %in ms
                        % Record the response
                        R.response(iTrial) = {KbName(keyCode)};%str2double(KbName(keyCode==1));
                        R.reactiontime(iTrial) = RT;
                    
                    % Wail until the end of the response window (1500ms)
                    WaitSecs((R.timeout-RT)/1000); %WaitSecs works in seconds
                    
                    % Break from the while loop
                    QuestionAnsweredBool = 0;
                    end

                elseif (GetSecs-R.Onsetres)*1000 > R.timeout
        
                % No answer was given with the correct mouse button within the response temporal window
                RT = R.timeout; %in ms
                R.response(iTrial) = {NaN};
                R.reactiontime(iTrial) = RT;
                                                
                R.missedButtonClick = 1;
                
                % Break from the while loop
                QuestionAnsweredBool = 0;
                end 

            % Escape pressed?    
            if KBpressed && keyCode(P.quitKey)
                QuestionAnsweredBool = 1;                  %Break from the while loop
                Screen('CloseAll')
            end 
            
        end 

PsychPortAudio('Stop', P.pahandle_stim, [], [], [], R.OnsetAV+0.500);

            %% present 3000 ms fixation each 9 communicative trials
            if S.order.comm(iTrial) == 1
                nTrial=nTrial+1;
                if nTrial==9
                    % do a final check of calibration using driftcorrection
                    %EyelinkDoDriftCorrection(el);
                    
%                     Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, P.cross_colour, [P.winCenter_x,P.winCenter_y],0);
%                     [~,R.OnsetBlockFix] = Screen('Flip',P.win, R.Onsetres + P.restime*P.ifi - P.ifi/2);
%         
%                     Screen('DrawLines', P.win, P.fixxCoords, S.draw.fixx_LineWidth, P.cross_colour, [P.winCenter_x,P.winCenter_y],0);
%                     [~,R.OffsetBlockFix] = Screen('Flip',P.win, R.OnsetBlockFix + P.fixation_in_flips*P.ifi - P.ifi/2);

                    DrawFormattedText(P.win, 'P', P.winCenter_x-fixxWidth_pix, P.winCenter_y+fixxHeight_pix, P.cross_colour);
                    Screen('DrawingFinished', P.win);
                    [~,R.OnsetBlockFix] = Screen('Flip',P.win);
                
                    t0p = GetSecs;
                    
                    while GetSecs - t0p <= 2.985
                [KBpressed, KbTime, keyCode] = KbCheck(-1)

                if KBpressed && keyCode(P.quitKey)
                   EyelinkDoDriftCorrection(el);
                   
                end
                    end
        
                    
                    DrawFormattedText(P.win, 'P', P.winCenter_x-fixxWidth_pix, P.winCenter_y+fixxHeight_pix, P.cross_colour);
                    Screen('DrawingFinished', P.win);
                    [~,R.OffsetBlockFix] = Screen('Flip',P.win, R.OnsetBlockFix + P.fixation_in_flips*P.ifi - P.ifi/2);

                    nTrial=0;
                    

                    R.BlockFix_duration = R.OffsetBlockFix - R.OnsetBlockFix;

                     %early responses feedback (behavioural only)
                    if z>=4
                        myText = 'ATTENZIONE: Rispondi solo quando vedi la lettera';
                        boundsText = Screen(P.win,'TextBounds',myText); %Returns the smallest enclosing rect for the drawn text [L,T,R,B]
                        DrawFormattedText(P.win, myText, round(P.winCenter_x-boundsText(3)/2), round(P.winCenter_y-boundsText(4)/2), P.text_colour);
                        Screen('Flip', P.win, []);
                        
                        WaitSecs(4);
                        z=0;  
                    end
                            
                                       
                end
                        
            end

            Screen('FillRect', P.win, P.black);
            [~,R.Offsetres]=Screen('Flip', P.win);
            
    
            %Log the durations of the flips, visual stimuli, ISIs, and PSI
            %in milliseconds
            R.ISI_durations(iTrial) = (R.OnsetAV-R.Onsetfix)*1000;
            R.stim_durations(iTrial) = (R.Onsetres-R.OnsetAV-0.017)*1000;
            %R.duration(iTrial) = (R.OnsetPSI-R.Onsetfix)*1000;
            R.res_duration(iTrial) = (R.Offsetres-R.Onsetres)*1000;
            R.asynch(iTrial) = (R.video_onset(iTrial)-R.audio_onset(iTrial))*1000;
            R.trial_duration(iTrial) = (R.Offsetres-R.Onsetfix)*1000;
            

 Screen('Close') %closes open textures
end
 
    % finish up: stop recording eye-movements,
    % close graphics window, close data file and shut down tracker
%     Eyelink('StopRecording');
    Eyelink('CloseFile');
    
    
     % download data file
    try
        fprintf('Receiving data file ''%s''\n', edfFile );
        status=Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        end
    catch %#ok<*CTCH>
        fprintf('Problem receiving data file ''%s''\n', edfFile );
    end

    % Shutdown Eyelink:
    Eyelink('Shutdown');
    % Close window:
    sca;
    
 Screen('CloseAll');   
 
 
 Log={};
 Log(:, 1) = Rtab(:,1);
 Log(:, 2) = Rtab(:,2);
 Log(:, 3) = Rtab(:,3);
 %Log(:, 4) = Rtab(:,4);
 Log(:, 5) = R.response';
 Log(:, 6) = num2cell(R.reactiontime');
 Log(:, 7) = R.early;

 Loglabel = {'communicative', 'response_cue', 'aud_loc', 'response', ...
 'RT', 'early'};
 Log = [Loglabel; Log];
 Filename = append('C:\Users\HPZ22018292\Documents\MATLAB\Fabic\data\Log_audio_', S.subj_ID, '_', S.run, '.xlsx');
 writecell(Log, Filename);
 %Filename1 = append('\Res_', S.subj_ID, '_', S.run, '.xlsx');
 %writecell(R.response', 'Res4_subj2_r4.xlsx');

