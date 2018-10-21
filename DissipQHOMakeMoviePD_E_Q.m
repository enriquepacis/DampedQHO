function varargout = DissipQHOMakeMoviePD_E_Q( varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

t = varargin{1};
PDEt = varargin{2};
PDQt = varargin{3};
Qval = varargin{4};
Qpoints = varargin{5};
VQ = varargin{6};
Eexp = varargin{7};
Qexp = varargin{8};

NReqArgIn = 8;

nt = length(t);
Nv = length(Qval);

maxPDQ = max(max(PDQt));
maxPDE = max(max(PDEt));

ntPlay = 60; % number of points in animation
waittime = 0.005;
MovieName = 'myMovie';
FPS = 8;
ShowTimePlots = 1;
ScalingFactorPDQ = 3;
CaptureFrames = [];

%% Override defaults

args = varargin(NReqArgIn+1:end);
while length(args) >= 2
    prop = args{1};
    val = args{2};
    args = args(3:end);
    switch prop
        case 'MovieName'
            MovieName = val;
            
        case 'Frames'
            ntPlay = val;
            
        case 'FPS'
            FPS = val;
            
        case 'ScalingFactorPDQ'
            ScalingFactorPDQ = val;
            
        case 'CaptureFrames'
            CaptureFrames = val;
            
        otherwise
            error(['DiisipQHOMakMoviePD_E_Q.M: ', prop, ...
                ' is an invalid property specifier.']);
    end
end

%% Secondary Calculations
Qmax = max(Qexp);
if Qmax == 0
    Qdiv = 1;
else
    Qdiv = Qmax;
end


if nt > ntPlay
    dnsmpl = floor(nt/ntPlay);
else
    dnsmpl = 1;
    ntPlay = nt;
end

dns_vect = dnsmpl*(0:ntPlay-1) + 1;

td = t( dns_vect ); % downsampled time vector

ntd = length(td);
disp(['Movie will have ', num2str(ntd), ' frames.']);


MovieFig = figure;

resize_prompt_string = 'Resize figure and press any key to continue.';
disp(resize_prompt_string);
vocalize(resize_prompt_string);
pause;

% scrSz = get(0, 'ScreenSize');
% set(MovieFig, 'Position', [0 0 scrSz(3:4)]);

movie_fig_pos = get(MovieFig, 'Position');
movie_rect = [0 0 movie_fig_pos(3:4)];

if ShowTimePlots
    % NewFig = figure;
    MarginL = 0.1; MarginR = 0.075; MarginC = 0.075;
    MarginT = 0.05; MarginB = 0.1;
    
    XL = MarginL;  XR = 0.4;
    dXL = XR - XL - MarginC;
    dXR = 1 - XR - MarginR;
    dXfull = 1 - MarginR - MarginL;
    
    YTop = 0.55; YBot = MarginB;
    
    dYTop = 1-YTop-MarginT;
    
    PDE_axes = axes('Position', [XL YTop dXL dYTop]);
    PDEQ_axes = axes('Position', [XR YTop dXR dYTop]);
    
    % PDEvt_axes = axes('Position', [XR MarginB dXfull 0.15]);
    PDEvt_axes = axes('Position', [XR MarginB dXR 0.15]);
    %PDQvt_axes = axes('Position', [XR MarginB + 0.175 dXfull 0.15]);
    PDQvt_axes = axes('Position', [XR MarginB + 0.175 dXR 0.15]);
    
    DblDotAx = axes('Position', [XL MarginB dXL 0.325]);
    
    AllAxes = [PDE_axes PDEQ_axes PDEvt_axes PDQvt_axes DblDotAx];
else
    PDE_axes = subplot(1,7,1:2);
    PDEQ_axes = subplot(1,7,4:7);
    
    
    AllAxes = [PDE_axes PDEQ_axes];
end

EnergyHi = ceil(2*sqrt(Nv))+ 0.5;

pct20 = round( ntd * 0.2 );
% tidx = str2double(get(handles.editTimeIdx, 'String'));
for sample_idx = 1:ntd
    
    didx = dns_vect(sample_idx); % data index
    
    % Plot the position distribution in handles.axes2
    axes(PDEQ_axes);
    cla(PDEQ_axes);
    plot(Qpoints, VQ, 'LineWidth', 2); grid on;
    hold on;
    
    ELadderXData = ones(Nv,1)*[Qval(1) Qval(end)];
    ELadderYData = (0.5+[0:Nv-1]')*[1 1];
    ELadder = plot([Qval(1) Qval(end)], ELadderYData, ...
        'LineWidth', 2, 'Color', [1 0 1]);
    
    EexpMark = plot([Qval(1) Qval(end)], Eexp(didx)*[1 1], ...
        'LineWidth', 3, 'Color', [0 0.8 0]);
    
    PDQBar = myFloatingBar( Qval, ...
        (ScalingFactorPDQ/maxPDQ)*PDQt(:,didx), ...
        Eexp(didx) );
    
    hold off;
    set(PDEQ_axes, 'FontSize', 24, 'FontName', 'Times', ...
        'GridLineStyle', '-');
    ylim([0, EnergyHi])
    
    % xlim([Qval(1) Qval(end)])
    % Position 0.5*(Nv+1) is the center, assuming Nv is odd
    Ncenter = 0.5*(Nv+1);
    dN = ceil(2*sqrt(Nv));
    xlim([Qval(Ncenter - dN) , ...
        Qval(Ncenter + dN )]);
    
    xlabel('$Q$ ($\AA$)', 'Interpreter', 'latex');
    ylabel('$E/\hbar \omega$', 'Interpreter', 'latex');    
    
    
    % Plot the energy distribution in PDE_axes
    axes(PDE_axes);
    cla(PDE_axes);
    barh(0.5+[0:Nv-1], PDEt(:,didx)); grid on;
    hold on;
    EexpMark2 = plot([0 1.05*maxPDE], Eexp(didx)*[1 1], ...
        'LineWidth', 3, 'Color', [0 0.8 0]);
    hold off;
    xlim([0 1.05*maxPDE])
    ylim([0, EnergyHi])
    set(PDE_axes, 'FontSize', 24, 'FontName', 'Times', 'GridLineStyle', '-');
    xlabel('Probability');
    ylabel('$E/\hbar \omega$', 'Interpreter', 'latex');
    
    
    if ShowTimePlots
        
        axes(PDEvt_axes);
        EexpData = plot(t, Eexp, 'LineWidth', 2); grid on;
        hold on;
        plot(t(didx)*[1 1], [0, 1.05*max(Eexp)], ...
            'LineWidth', 2, 'LineStyle', '--', 'Color', [1 0 0]);
        hold off;
        xlim([t(1) t(dns_vect(end))]);
        ylim([0, 1.05*max(Eexp)]);
        set(PDEvt_axes, 'FontSize', 24, 'FontName', 'Times', ...
            'GridLineStyle', '-');
        xlabel('$t/\tau$', 'Interpreter', 'latex');
        ylabel('$\left< E \right>/\hbar \omega$', 'Interpreter', 'latex');
        
        axes(PDQvt_axes);
        pcolor(t, Qval, real(PDQt)); shading interp; lighting phong;
        set(PDQvt_axes, 'FontSize', 24, 'FontName', 'Times', ...
            'GridLineStyle', '-', 'XTick', []);

        ylabel('$Q$ ($\AA$)', 'Interpreter', 'latex');
        hold on;
        plot(t(didx)*[1 1], [Qval(Ncenter - dN) , ...
        Qval(Ncenter + dN )], ...
            'LineWidth', 2, 'LineStyle', '--', 'Color', [1 1 1]);
        hold off;
        xlim([t(1) t(dns_vect(end))]);
        ylim([Qval(Ncenter - dN) , ...
            Qval(Ncenter + dN )]);
        
        axes(DblDotAx);
        
        drawDoubleDot_v2( DblDotAx, 0.5*[1 1], Qexp(didx)/Qmax, ...
            'ChargeColor', 'w');
        
    end
    
    drawnow;
    
    DissipQHO_EnergyMovie(sample_idx) = getframe(MovieFig, movie_rect);
    
    if ~mod(sample_idx, pct20) && sample_idx ~= ntd
        pct_complete = round(100*(sample_idx/ntd));
        vocalize(['Movie is ', num2str(pct_complete), ' percent complete.']);
    end
    
end

close(MovieFig);

videoFile = VideoWriter(MovieName, 'MPEG-4');
open(videoFile);
writeVideo(videoFile, DissipQHO_EnergyMovie);
close(videoFile);

% Old movie2avi syntax
% movie2avi(DissipQHO_EnergyMovie, MovieName, 'FPS', FPS);

varargout{1} = 0;

end

