% Function is fed with inputs from MoL_Script
function [cumul_release, output] = profile_1D_Sph(DrugSol,tSol,DrugCDSol,Drug_0,...
                    R0,r,Volume,Figure_1,Figure_2,Figure_3,...
                    Figure_4,Figure_5,Animation,UploadData)
%% Extract [Drug] Throughout Hydrogel

% Find average [Drug] across all spatial slices at time (t)
release_Drug = mean(DrugSol,2); 
% Find average [Drug-CD] across all spatial slices at time (t)
release_CD = mean(DrugCDSol,2);
% The [Drug total] = [Drug] + [Drug-CD] = unbound drug + bound drug
release = release_Drug + release_CD;
% Released [Drug] = initial [Drug] - final [Drug] (inside the hydrogel)
cumul_release = release(1)-release(:);
% Normalize released [Drug] to the initial [Drug] across all spatial slices
frac_release = cumul_release./Drug_0; 
% Calculate the percentage (%) of [Drug] remaining at the last time point
cumul_frac = 100.*(release(1)-release(end))./Drug_0;

thres = 0.9;  % Threshold to consider the release profile complete (90%)
% Default values for time of 90% release (tBurst) and equlibrium time
% (tStop)
tBurst = 0;    
tStop = 0; 
% Counters
count = 1;
count2 = 1;

for i = 2:size(frac_release,1)  % Iterate throughout all timepoints of the 
   % release profile
   % Calculate the time required to reach 90% drug release 
   if frac_release(i) >= thres % If the hydrogel reaches 90% drug release 
       tBurst(count) = tSol(i);
       count = count + 1;
   end 
    
   % Calculate the time required to reach equilibrium based on the
   % specified thresholds (tolerances)
   if abs(release(i)-release(i-1)/(tSol(i)-tSol(i-1))) <= 0.005 
       % The profile slope is ~ zero
       tStop(count2) = tSol(i-1);
       count2 = count2 + 1;
   end  
end 

tBurst = tBurst(1); % Time to Reach 90% Drug Release is the first value 
                    % of tBurst
tStop = tStop(1); % Equilibrium Time is the 1st value of tStop

% Cumulative Release Profile
if Figure_1 ~= 0    % If Figure 1 is selected by the user
    figure
    plot(tSol, cumul_release*10^3); % Plot the Cumulative Release Profile
    hold on
    % Plot a dashed red line to represent the Time of 90% Drug Release
    if tBurst ~= 0
        plot(tBurst*ones(1,5),linspace(0,max(release(:)),5)*10^3, '--r',...
        'LineWidth',2);
        hold on
    end 
    % Plot a dashed green line to represent the Initial [Drug]
    plot(linspace(0, tSol(end), 5),linspace(Drug_0, Drug_0, 5)*10^3,'--g',...
        'LineWidth',2);
    hold off
    title('Cumulative Release Profile')
    xlabel('Time (hr)')
    ylabel('Drug Release (mM)')
    legend('Model Output','90% Drug Release','Initial Drug Conc.',...
        'Location','Southeast')
    savefig('Fig1') % Automatically save figure in working folder
end

% Fractional Release Profile
if Figure_2 ~= 0  % If Figure 2 is selected by the user
    figure
    plot(tSol, frac_release);  % Plot the Fractional Release Profile
    hold on
    % Plot a dashed red line to represent the Time of 90% Drug Release
    if tBurst ~= 0
        plot(tBurst*ones(1,5), linspace(0,max(frac_release(:,1)),5), '--r', ...
        'LineWidth',2);
    end 
    hold off    
    title('Fractional Release Profile')
    xlabel('Time (hr)');
    ylabel('Fractional Drug Release (%)');
    legend('Model Output','90% Drug Release','Location','Southeast')
    savefig('Fig2') % Automatically save figure in working folder
end 

% Change Between Points in Release Profile
if Figure_4 ~= 0    % If Figure 4 is selected by the user
    figure
    % Plot the point-to-point difference of the Fractional Release Profile
    plot(tSol(2:end,:), diff(frac_release)); 
    hold on
    % Plot a dashed red line to represent the Time of 90% Drug Release
    if tBurst ~= 0
        plot(tBurst*ones(1,5), linspace(0,max(diff(frac_release(:,1))),5), ...
        '--r','LineWidth',2);
    end 
    hold off
    title('Change in Fractional Release Profile')
    xlabel('Time (hr)');
    ylabel('Change in Release (\Delta %)');
    legend('Model Output','90% Drug Release')
    savefig('Fig3') % Automatically save figure in working folder
end 

slope = zeros(size(tSol)); % Preallocate space
for i = 2:numel(tSol)   % Iterate throughout the time points
    % Calculate the slope point-by-point
    slope(i)=((release(i)-release(i-1)).*Volume)/(tSol(i)-tSol(i-1));
end 

% Change in Slope
if Figure_3 ~= 0    % If Figure 3 is selected by the user
    figure
    plot(tSol,slope);   % Plot the point-by-point slope vs. time
    hold on
    % Plot a dashed red line to represent the Time of 90% Drug Release
    if tBurst ~= 0
        plot(tBurst*ones(1,5),linspace(0,max(slope),5), '--r','LineWidth',...
            2);
    end 
    hold off
    title('Change in Release Rate')
    xlabel('Time (sec)')
    ylabel('Slope (\Delta % / hr)')
    legend('Model Output','90% Drug Release')
    savefig('Fig4') % Automatically save figure in working folder
end

% Directional Release Profile (Heatmap)
if Figure_5 ~= 0    % If Figure 5 is selected by the user
    Fig5 = figure;
    xaxis = [0 tSol(end)];   % Time is the x-axis
    yaxis = linspace(0,r*100,R0+1);  % Space (r-dir) is the y-axis
    imagesc(xaxis,yaxis,DrugSol');  % Create heatmap of [Drug] in hydrogel
    xlabel('time span (hr)')
    ylabel('axial position (mm)')
    ticks = 6;      % Define number of xticks
    % Specify xtick labels
    time_labels = tSol(ceil(linspace(1,numel(tSol),ticks)));
    for i = 2:ticks % Iterate through xticks
        % Round the xtick values to 1 SigFig
        if ceil(abs(log10(time_labels(i)))) < 2 || ceil(log10(time_labels(i))) >= 0
           time_labels(i) = round(time_labels(i),1);
        else 
            time_labels(i) = round(time_labels(i),ceil(abs(log10(time_labels(i)))));
        end
    end 
    xticklabels(time_labels);
    colormap jet
    c = colorbar;
    c.Label.String = 'Molar Drug Conc. (mM)';
    set(gca,'YDir','normal')
    savefig('Fig5')  % Automatically save figure in working folder
    hold on
    % Automatically close the figure after 5 seconds of show
    % Note: This is done to prevent the figure from being overwritten by 
    % the animation 
    jFrame = get(handle(Fig5), 'JavaFrame');
    pause(5)
    jFrame.setMinimized(1);
end 

%% Upload Data

file = 0; % Default value
if UploadData ~=0   % If the user wants to upload experimental data
     % Filter for MALTAB, Excel, txt, or CSV files
     filter = {'*.xlsx';'*.m';'*.mat';'*.txt';'.*csv';'*.*'};
     % Prompt user to select the data file
     [file,pathname] = uigetfile(filter,'Select a data file');
     if file ~= 0   % If the filename is NOT empty
        data = importdata(strcat(pathname, file));  % Import data file
        t = data(:,1);  % Collect time data (hr) from the 1st column
        y = data(:,2);  % Collect release data (%) from the 2nd column

        % Preallocate space
        yCalc = zeros(size(y));
        tCalc = zeros(size(y));

        for i = 1:numel(t)  % Iterate through all time points in the data 
            % Find the closest time point in the model output to the
            % time point in the experimental data
            t_vec = abs(tSol-(t(i)));
            j = find(t_vec == min(t_vec));
            yCalc(i) = frac_release(j); % Find the modeled release that
            % correlates this matched time value
            tCalc(i) = tSol(j); % Return the matched modeled time value
        end 

        % Fractional Release Profiles (Model Output & Experimental Data)
        figure
        plot(tCalc,yCalc,'b',t,y,'r*')
        legend('Model Output','Experimental Data','Location','SouthEast')
        xlabel('Time (hr)')
        ylabel('Fractional Drug Release (%)')
        RMSE = sqrt((mean(y-yCalc).^2));
        Rsquared = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2);
        savefig('Fig6') % Automatically save figure in working folder
     end
end 
%% Save txt output

if tBurst == 0     % If 90% Drug Release is NOT reached
    tBurst_str = '90% Drug Release Not Reached';    
else 
    tBurst_str = strjoin({'Time for 90% Drug Release = ',num2str(tBurst),...
        ' hr'});
end 

if tStop == 0 || tStop < tBurst    % If Equilibrium is NOT reached
    tStop_str = 'Equilibrium Not Reached';
else 
    tStop_str = strjoin({'Time to Reach Equilibrium = ',num2str(tStop),' hr'});
end

RMSE_str = '';     % Default value of RMSE (no uploaded data)
Rsquared_str = '';     % Default value of R-squared (no uploaded data)
if UploadData ~= 0 & file ~= 0     % If there is uploaded data
    RMSE_str = strjoin({'RMSE: ',num2str(RMSE)});
    Rsquared_str = strjoin({'R-squared: ',num2str(Rsquared)});
end 

% Output text as a string array
output = {strjoin({'Total Molar Release = ', ...
    num2str((release(1)-release(end))*10^3),'mM '}),...
    strjoin({'Total Amount of Drug Released = ',...
    num2str((release(1)-release(end))*Volume*10^3),'mmoles'}),...
    strjoin({'Percentage of Loaded Drug Released = ',...
    num2str(cumul_frac),'%'}),tBurst_str,tStop_str,RMSE_str, Rsquared_str};
%% 3D Animation

if Animation ~=0   % If the user wants an Animation
    % Prompt the user to save the animation as a .mp4 file
    [filename,~] = uiputfile('*.mp4','Save Hydrogel Animation');
    transparency = 0.6; % Default animation transparency
    shade = 0.3;    % Default animation shading 
    framerate = 1;   % Default animation framerate
    % Default animation view point
    xpos = 15;
    ypos = 15;
    zpos = 3;
    settings = Animation_Properties; % Call the Animation_Properties app
    pause(10)      % Pause to all user to modulate animation parameters
    if isvalid(settings)
        if settings.OKButton.Value == 1 % If the 'OK' button is pressed
            % Save the changed animation parameters
            transparency = settings.TransparencySpinner.Value;
            shade = settings.ShadingSpinner.Value;
            framerate = settings.FramerateSpinner.Value;
            xpos = settings.XSlider_2.Value;
            ypos = settings.YSlider.Value;
            zpos = settings.ZSlider.Value;
        end 
    end
    
    delete(settings)    % Delete the saved settings
    idx = 0;    % Counter
    fig1 = figure;      % Delcare new figure
    for k = 1:numel(tSol)    % Iterate through all time points in model
        clf    % Clear fig1
        idx = idx + 1;  % Increment the counter
        [x,y,z] = sphere;   % Make a sphere (geometry will be modified)
        drug = zeros(size(x));  % Preallocate space  
        for j = 1:R0    % Iterate through all spatial slices (r-dir)
            rho = r*j/R0*10^3; % Convert the radius to mm
            % Modify sphere geometry based on radius of slice
            x1 = x*rho;
            y1 = y*rho;
            z1 = z*rho;
            drug(:) = DrugSol(k,j); % Extract [Drug] values
            % Make a 3D plot of sphere
            surf(x1,y1,z1,drug,'EdgeAlpha',transparency);
            hold on
        end 
        shading interp
        hold off
        alpha(shade)
        c = colorbar;
        caxis([0 DrugSol(1)])  
        c.Label.String = 'Molar Free Drug Conc. (mM)';
        t_k = tSol(k);
        title(['Hydrogel at t = ',num2str(t_k),' hrs'])
        xlabel('X-direction (mm)')
        ylabel('Y-direction (mm)')
        zlabel('Z-direction (mm)')
        xlim([-2*r*10^3 2*r*10^3])
        ylim([-2*r*10^3 2*r*10^3])
        zlim([-2*r*10^3 2*r*10^3])
        view([xpos ypos zpos]) 
        fig2 = fig1;    % Save an empty figure as the previous figure in 
                    % case fig1 is closed prematurely
                          
        % Force MATLAB to draw the image at this point
        drawnow 
        
        if isvalid(fig1) ~= 0  % If fig1 is not prematurely closed
            % Flipbook of images to save as video struc
            flipBook(idx) = getframe(fig1, [10, 10, 520, 400]); 
        else 
            fig1 = fig2;    % Use fig2 to re-establish fig1
        end
    end 
    
    savefig('fig7')    % Automatically save figure in working folder 
    
    if filename == 0   % If the filename is empty use default filename
        filename = '1D_Sphere_Hydrogel_Animation.mp4';
    end
    
    % Save animation as an .mp4 file 
    myWriter = VideoWriter(filename,'MPEG-4');
    myWriter.FrameRate = framerate;    % Specify animation framerate
    
    % Open the VideoWriter object, write the movie, and close the file
    open(myWriter);
    writeVideo(myWriter, flipBook);
    close(myWriter);
end 
end
