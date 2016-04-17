% Slapbass synthesis
% Programed by Wen Xiao

N = 44100;  % Implementation time number
fs = 44100; % Sample rate
c = 347.23; % Speed of sound (m/sec)

% Fret Position

fp = 0.6;

% Fret Height

fret_h = -0.6;

% Length of string (m)

L = 0.864;

% Delay Line Length

delay = round( L * fs / c);
delay1 = round (fp * delay);
delay2 = delay - delay1;

% Reading position

read_pos = 0.3;
read_pos = round(read_pos * delay);
if read_pos < delay1
    rm = read_pos+1;      % Minus direction reading pointer
    rp = delay1-read_pos; % Plus direction reading pointer
    read_bound = delay1;
else
    rm = (read_pos-delay1)+1;      % Minus direction reading pointer
    rp = delay2-(read_pos-delay1); % Plus direction reading pointer
    read_bound = delay2;
end

% Plucking/Striking position

pl = 0.5;
pl = round(delay * pl);

% Initialization
% Two sections are established. The junction position is where the fret
% locates. Each section has two delay lines.

%dp = 0.5*ones(1, pl);       % Stricking Initialization
dp = linspace(0,0.5,pl);   % Plucking(Displacement) Initialization
dp = [dp linspace(0.5,0,delay-pl)];
dp1 = dp(1:delay1);
dp2 = dp((delay1+1):(delay1+delay2));
dp1 = fliplr(dp1);
dp2 = fliplr(dp2);

%dm = 0.5*ones(1, pl);       % Stricking Initialization
dm = linspace(0,0.5,pl);   % Plucking(Displacement) Initialization
dm = [dm linspace(0.5,0,delay-pl)];
dm1 = dm(1:delay1);
dm2 = dm((delay1+1):(delay1+delay2));

pointer1 = 1;  % Refreshing Pointer at section 1
pointer2 = 1;  % Refreshing Pointer at section 2

obs_pos1 = delay1; % Junction sample of minus delay line in section 1
obs_pos2 = delay2; % Junction sample of plus delay line in section 2

st = 0; % Status button to see if the string is touching the fret(1 is touching and 0 is not).

p = zeros(1,N); % Output Buffer

z1 = 0; % Filter initial condition
z2 = 0; % Filter initial condition

displayp = zeros(1, delay); % Buffer for animation
displaym = zeros(1, delay); % Buffer for animation

for i = 1:N,
    
    % Read the samples at the given position and pass it through a lowpass
    % filter
    p(i) = dm1(rm) + dp1(rp);
    [p(i) , z2] = filter([1 1], 3, p(i), z2);
    
    % Store the differences of positive at the junction point
    tp = dp1(pointer1) - dp2(obs_pos2);
    tm = dm1(obs_pos1) - dm2(pointer2);
    
    % Store the current value of reading points 
    dm1c = dm1(pointer1);
    dp1c = dp1(pointer1);
    dm2c = dm2(pointer2);
    dp2c = dp2(pointer2);
    
    % Judge if the string is touching the fret and do the corresponding
    % calculation
    if dp1(pointer1) + dm1(obs_pos1) > fret_h
        
        % Judge if it is the status the string is going to leave the fret
        % and make the compensation if so
        if st == 1
            p_diff = dp1(pointer1) - dp1(pointer1+1);
            dp1 = dp1 - tp - p_diff;
            dm1 = dm1 - tm + p_diff;
                
            dm1c = dm1(pointer1); % Refresh the stored current value
            dp1c = dp1(pointer1); % Refresh the stored current value
            st = 0;               % Refresh the string status
        end
        % Do the calculation when the string leaves the fret
        dp1(pointer1) = -dm1c;
        dm1(pointer1) = dm2c;
        dp2(pointer2) = dp1c;
        [dm2(pointer2), z1] = filter(-0.49, [1 -0.5], dp2c, z1);
    
    else
        % Do the calculation when the string touches the fret
        dp1(pointer1) = -dm1c;
        dm1(pointer1) = fret_h - dp1c;
        dp2(pointer2) = fret_h - dm2c;
        [dm2(pointer2), z1] = filter(-0.49, [1 -0.5], dp2c, z1);
        
        st = 1; % Refresh the string status
    
    end
    
    pointer1 = pointer1 + 1;       % increment "pointer" 
    if pointer1 > delay1           % check delay-line limits 
        pointer1 = 1;
    end
    
    pointer2 = pointer2 + 1;       % increment "pointer" 
    if pointer2 > delay2           % check delay-line limits 
        pointer2 = 1;
    end
    
    rm = rm + 1;                   % increment "pointer" 
    if rm > read_bound             % check delay-line limits 
        rm = 1;
    end
    
    rp = rp + 1;                   % increment "pointer" 
    if rp > read_bound             % check delay-line limits 
        rp = 1;
    end
    
    obs_pos1 = obs_pos1 +1;        % increment "pointer" 
    if obs_pos1 > delay1           % check delay-line limits 
        obs_pos1 = 1;
    end
    
    obs_pos2 = obs_pos2 +1;        % increment "pointer" 
    if obs_pos2 > delay2           % check delay-line limits 
        obs_pos2 = 1;
    end
    
    % Here is a small animation to show the evolution of string vibration.
    % The blue line is denoted to positive direction wave and the red line is
    % demoted to negative direction wave. Yellow line is the sum of
    % positive one and negative one. We can pause the animation by changing
    % the "if" condition and pause()command.
    
    %{
    if rem(i, 5) == 0
    %if i == 24.6*5
    %if i == 24.6*5
    num1 = pointer1;
    for n= 1:delay1
        displayp(n) = dp1(num1);
        displaym(n) = dm1(num1);
        num1 = num1 + 1;
        if num1 > delay1;
            num1 = 1;
        end
    end
    num2 = pointer2;
    for n = 1:delay2
        displayp(delay1 + n) = dp2(num2);
        displaym(delay1 + n) = dm2(num2);
        num2 = num2 + 1;
        if num2 > delay2;
            num2 = 1;
        end
    end
    displayp(1:delay1) = fliplr(displayp(1:delay1));
    displayp(delay1+1:delay1+delay2) = fliplr(displayp(delay1+1:delay1+delay2));
    display_all = displayp + displaym;
    plot(displayp,'.','MarkerSize',20,'MarkerFaceColor','b')
    %plot(displayp,'Linewidth',3,'MarkerFaceColor','b')
    hold on 
    plot(displaym,'.','MarkerSize',20,'MarkerFaceColor','r')
    %plot(displaym,'Linewidth',3,'MarkerFaceColor','r')
    hold on 
    plot(display_all,'.','MarkerSize',20,'MarkerFaceColor','y')
    %plot(display_all,'Linewidth',3,'MarkerFaceColor','y')
    hold on
    plot([delay1 delay1],[-1,fret_h],'k','Linewidth',2)
    hold off
    
    pause(1);
    %pause();
    end
    %}
end

%sound(p,44100)
plot(p(1:5000));