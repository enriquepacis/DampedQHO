function varargout = DissipativeQHOTool(varargin)
% DISSIPATIVEQHOTOOL MATLAB code for DissipativeQHOTool.fig
%      DISSIPATIVEQHOTOOL, by itself, creates a new DISSIPATIVEQHOTOOL or raises the existing
%      singleton*.
%
%      H = DISSIPATIVEQHOTOOL returns the handle to a new DISSIPATIVEQHOTOOL or the handle to
%      the existing singleton*.
%
%      DISSIPATIVEQHOTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISSIPATIVEQHOTOOL.M with the given input arguments.
%
%      DISSIPATIVEQHOTOOL('Property','Value',...) creates a new DISSIPATIVEQHOTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DissipativeQHOTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DissipativeQHOTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DissipativeQHOTool

% Last Modified by GUIDE v2.5 01-Aug-2014 14:29:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DissipativeQHOTool_OpeningFcn, ...
                   'gui_OutputFcn',  @DissipativeQHOTool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DissipativeQHOTool is made visible.
function DissipativeQHOTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DissipativeQHOTool (see VARARGIN)

% Choose default command line output for DissipativeQHOTool
handles.output = hObject;

MainFig = gcf;
hbar = 6.582119E-16;     % [eV*s] Reduced Plank constant
c = 2.99792E17;           % [nm/s]
kB = 8.617332478E-5;   % [eV/k] Boltzmann constant
setappdata(MainFig, 'hbar', hbar);
setappdata(MainFig, 'c', c);
setappdata(MainFig, 'kB', kB);

FontName = 'Times';
FontSizeLabel = 18;
setappdata(MainFig, 'FontSizeLabel', FontSizeLabel);
setappdata(MainFig, 'FontSizePlotLabel', 24)
setappdata(MainFig, 'FontName', FontName)

N = str2double(get(handles.editDimension, 'String'));

projectorCombination = normalizevector(rand([N,1])*2 - 1);
setappdata(MainFig, 'ProjComb', projectorCombination);

convertAxes2Label(handles.axesT1Label, '$\tau/T_1$', ...
    'FontSize', FontSizeLabel, 'FontName', FontName, ...
    'HorizontalAlignment', 'right');

convertAxes2Label(handles.axesT2Label, '$\tau/T_2$', ...
    'FontSize', FontSizeLabel, 'FontName', FontName, ...
    'HorizontalAlignment', 'right');

convertAxes2Label(handles.axesFreqLabel, '$f$ (cm$^{-1})$', ...
    'FontSize', FontSizeLabel, 'FontName', FontName, ...
    'HorizontalAlignment', 'right');

convertAxes2Label(handles.axesMassLabel, '$m$ (amu)', ...
    'FontSize', FontSizeLabel, 'FontName', FontName, ...
    'HorizontalAlignment', 'right');

Path.home = pwd;
Path.data = pwd;
setappdata(MainFig, 'Path', Path);


UpdateQHO(handles);

popupSelectDissipation_Callback( handles.popupSelectDissipation, [], handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DissipativeQHOTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DissipativeQHOTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ========================================================================
% ========================vvv UTILITIES vvv===============================
% ========================================================================

% --------------------------------------------
function SaveNameBase  = MakeSaveName(handles)
MainFig = gcf;

MassString = [ num2str_nodot(textboxval(handles.editMass)),  'amu'];
FreqString = [ num2str_nodot(textboxval(handles.editFreq)),  'cminv'];
NString = [ num2str_nodot(textboxval(handles.editDimension)) ];
TimeString = [ 'nT', num2str_nodot( textboxval( handles.editDuration) ) ];

figure(MainFig);
DissipationCase = get(handles.popupSelectDissipation, 'Value');
switch DissipationCase
    case 1 % NONE
        EnvCplg = 'EnvNone';
    case 2 % Relaxation
        EnvCplg = ['EnvR_Ta', ...
            num2str_nodot( 1/textboxval(handles.editAlpha1 ) )];
    case 3 % Dephasing
        EnvCplg = ['EnvDph_Tb', ...
            num2str_nodot( 1/textboxval(handles.editAlpha2 ) )];
    case 4 % Dephasing with dissipation 
        EnvCplg = ['EnvRDph_Ta', ...
            num2str_nodot( 1/textboxval(handles.editAlpha1 ) ), ...
            'Tb', ...
            num2str_nodot( 1/textboxval(handles.editAlpha2 ) )];
    case 5 % Thermalization
        EnvCplg = ['EnvThrmT', num2str_nodot(textboxval(handles.editTemperature)) , ...
            'K_Ta', ...
            num2str_nodot( 1/textboxval(handles.editAlpha1 ) ), ...
            'Tb', ...
            num2str_nodot( 1/textboxval(handles.editAlpha2 ) )];
end

SaveNameBase = ['DissQHO_', NString, '_', TimeString, '_', ...
    FreqString, '_', MassString, '_', EnvCplg] ;



% --------------------------------------------
function PrintFigure(handles, MainFig, TargetFig, savename)

if strcmp(get(handles.PrintFigsMI, 'Checked'), 'on')
    
    LastSaveDir = getappdata(MainFig, 'LastSaveDir')
    
    if ~isempty(LastSaveDir)
        SavePathStart = LastSaveDir;
    else
        Path = getappdata(MainFig, 'Path');
        SavePathStart = Path.data;
    end
    
    LastSavePath = printeps(TargetFig, SavePathStart, savename)
    
    if LastSavePath ~= 0
        setappdata(MainFig, 'LastSaveDir', LastSavePath);
    else
        disp('User pressed CANCEL.');
    end
    
end

%-vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv----------------------------------------
function fighandle = GetFigHandle(handles)
% Returns the handle to the last figure in gFigList for use with the new
% graphic. The last figure either is a new figure, or it is a used figure,
% depending on the 'Checked' status of handles.ReuseFigsMI.

MainFig = gcf;
FigList = getappdata(MainFig, 'FigList');

ReuseFig = get(handles.ReuseFigsMI, 'Checked');
switch ReuseFig;
    case 'on'
        if length(FigList) % cFigList is not empty
            if ishandle(max(FigList)) % use the last figure in cFigList
                figure(max(FigList));
                clf(max(FigList));
                fighandle = max(FigList);
            else                       % make/add a new figure to cFigList
                fighandle = figure;
                FigList = [FigList fighandle];
                
            end
        else
            fighandle = figure;
            FigList = [FigList fighandle];
            
        end
    case 'off'
        fighandle = figure;
        FigList = [FigList fighandle];
end

setappdata(MainFig, 'FigList', FigList);

%-^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^----------------------------------------
%- FIGURE LIST MANAGEMENT FUNCTION ----------------------------------------
%--------------------------------------------------------------------------


% ------------------------------------------------------------------------
function UpdateQHO(handles)
% UPDATEQHO( handles ) updates the data representing the QHO. The dimension
% is obtained, and the Hamiltonian (H), the position operator (Q), the
% momentum operator (P), all are obtained.
%
MainFig = gcf;
hbar = getappdata(MainFig, 'hbar'); % [eV*s] Reduced Plank constant
c = getappdata(MainFig, 'c'); % Speed of light [nm/s]
c_cm = c*(1E-7); % Speed of light [cm/s]

%% Get constants from MainFig
N = str2double(get(handles.editDimension, 'String'));

%% Default QHO Parameters
f = str2double(get(handles.editFreq, 'String')) * c_cm; % Oscillator frequency. Given as 298 cm^(-1), converted to Hz
w = 2*pi*f; % Oscillator frequency. [rad/s]
hbw = hbar*w;
mc2 = 931.49E6*str2double(get(handles.editMass, 'String'));  % QHO mass*c^2. Given as 5.58 amu, converted to eV
m = mc2/c^2; % QHO mass. [eV*s^2/nm^2]

%% Establish operators
I_QHO = eye(N);   % Identity operator
d = sqrt(1:N-1);
a = zeros(N) + diag(d, +1);              % [] Annihilation operator
ad = a';                                  % [] Creation operator
H = hbar*w*(ad*a + 0.5*I_QHO);
Q = sqrt( hbar / ( 2 * m * w) ) * (ad + a);    % [nm] Position
P = 1i * sqrt( (m * hbar* w )/2 ) ...
    * (ad - a);
mw2 = m*w^2; %  [eV/nm^2];

ReportString = ['$\hbar \omega = ', num2str(1E3*hbw, '%0.4g'), '$ meV; $m \omega^2 = ', ...
    num2str(mw2, '%0.4g'), '$ eV/nm$^2$;'];
convertAxes2Label(handles.axesReport, ReportString, 'FontSize', 16);

[psi_stat, E_stat] = eig(H);
[Uqe, Eqe] = eig(Q);
UQE = Uqe'; % A transformation matrix for a eigenbasis vector to the position basis
            % psi_Q = UQE * psi_E

setappdata(MainFig, 'm', m);
setappdata(MainFig, 'f', f);
setappdata(MainFig, 'H', H);
setappdata(MainFig, 'a', a);
setappdata(MainFig, 'Q', Q);
setappdata(MainFig, 'P', P);
setappdata(MainFig, 'UQE', UQE);
setappdata(MainFig, 'Qval', diag(Eqe));
setappdata(MainFig, 'Eval', diag(E_stat));

%% Establish nine options for the initial condition
psi0 = zeros(N, 9); % Eight initial wavefunctions

[psi_stat, E_stat] = eig(H); % Find the eigenstates


% Eigenstate (ground state)
psi0(:,1) = psi_stat(:,1);

% Eigenstate (first excited state)
psi0(:,2) = psi_stat(:,2);

% Eigenstate (second excited state)
psi0(:,3) = psi_stat(:,3);

% COHERENT STATES: IC #7-9
% Coherent state, displaced by alpha = -2
alpha = -2;
D = expm( alpha * ad - conj(alpha)*a);
psi0(:,7) = D*psi_stat(:,1); % Displace the ground state

% Coherent state, displaced by alpha = 2
alpha = 2;
D = expm( alpha * ad - conj(alpha)*a);
psi0(:,8) = D*psi_stat(:,1); % Displace the ground state

% Coherent state, displaced by alpha = 1
alpha = 2;
D = expm( alpha * ad - conj(alpha)*a);
psi0(:,9) = D*psi_stat(:,1); % Displace the ground state

% ARBITRARY STATES: IC #4-6
% Arbitrary state #1
c1 = 0.75; c2 = sqrt(1-c1^2);
psi0(:,4) = normalizevector(c1* psi_stat(:,3) - ...
    c2 * psi_stat(:,1));

% Arbitrary state #2: linear combination of eigenstate and coherent state
c1 = 0.75; c2 = sqrt(1-c1^2);
psi0(:,5) = normalizevector(c1* psi_stat(:,4) - ...
    c2 * psi0(:,7));

% Arbitrary state #3: linear combination of eigenstate and coherent state
c1 = 0.75; c2 = sqrt(1-c1^2);
psi0(:,6) = normalizevector(c1* psi_stat(:,5) - ...
    c2 * psi_stat(:,8));

setappdata(MainFig, 'psi0', psi0);

CalculateEvolution( handles )


% ------------------------------------------------------------------------
function CalculateEvolution(handles)

MainFig = gcf;

hbar = getappdata(MainFig, 'hbar'); % [eV]
N = str2double( get(handles.editDimension, 'String') ); % SHO dimension
SimulationLength = str2double( get(handles.editDuration, 'String') ); % SHO periods
T = 1/getappdata(MainFig, 'f'); % QHO period [s]
H = getappdata(MainFig, 'H'); % [eV]
Q = getappdata(MainFig, 'Q'); % [nm]
UQE = getappdata(MainFig, 'UQE'); % Transformation from eigenbasis to position basis
alpha1 = str2double( get( handles.editAlpha1, 'String' ) );
alpha2 = str2double( get( handles.editAlpha2, 'String' ) );
f = getappdata(MainFig, 'f'); % [eV]
w = 2*pi*f;

[V, E_stat] = eig(H);      % V are stationary states
es = diag(E_stat);         % Energy of stationary states
Vinv = inv(V);

%% Get the initial state
psi0 = getappdata(MainFig, 'psi0');
psi0selection = get(handles.popupSelectInitialState, 'Value');
rho0 = outerprod(psi0(:,psi0selection), psi0(:,psi0selection));

%% Time evolution for non-dissipative case

CalcStart = tic;
textStatusMsg(handles.textStatus, ...
    'Calculating time evolution. Please wait ...', [1 0 0]);

% waitfig = dialog('WindowStyle', 'modal', 'Name', 'Please wait ...');
figure(MainFig);
DissipationCase = get(handles.popupSelectDissipation, 'Value');
switch DissipationCase
    case 1 % NONE
        TimeSamples = str2double( get(handles.editSamples, 'String') );
        t = linspace(0, SimulationLength*T, SimulationLength*TimeSamples);
        nt = length(t);
        
        set(handles.textTimeIdxMax, 'String', num2str(nt) );
        
        UpdateEditboxDiscreteFromSlider( handles.sliderTimeIdx, ...
            handles.editTimeIdx, handles.textTimeIdxMin, ...
            handles.textTimeIdxMax);
        
        Qdistrib = zeros(N, nt);
        Edistrib = zeros(N, nt);
        Eexp = zeros(1, nt);
        Qexp = zeros(1, nt);
        SVN = zeros(1, nt);
        for tind = 1:nt
            rhot = V*diag(exp(-1i*es*t(tind)/hbar))*Vinv*rho0*...
                V*diag(exp(1i*es*t(tind)/hbar))*Vinv;
            
            Edistrib(:, tind) = real( diag(rhot) );
            Qdistrib(:, tind) = real( diag( UQE*rhot*(UQE') ) );
            SVN(tind) = entropyVonNeumann( rhot );
            Eexp(tind) = real(trace(H*rhot));
            Qexp(tind) = real(trace(Q*rhot));
            
        end % END [ for tind = 1:length(t) ]
        
        % disp('Coherent time evolution is calculated.')
        textStatusMsg(handles.textStatus, ...
            'Calculated coherent time evolution', [0 0 0]);
        
    case 2 % Dissipation only
               
        a = getappdata(MainFig, 'a');
        L = sqrt(alpha1/T)*a;
        
        DissipType = 'relaxation only';
        
    case 3 % Dephasing only
        
        DephasingOnlyMode = get(handles.popupSelectDephasing, 'Value');
        
        switch DephasingOnlyMode
            case 1
                L = sqrt(alpha2/T)*(1/(hbar*w))*H;
                
                DissipType = 'dephasing only - Hamiltonian';
            case 2
                ProjComb = getappdata(MainFig, 'ProjComb');
                M = diag(ProjComb);
                szPC = size(ProjComb)
                szM = size(M)
                L = sqrt(alpha2/T)*M;
                
                DissipType = 'dephasing only - Random combination of projectors.';
        
        end
        
    case 4 % Dephasing and Dissipation
        N = str2double( get(handles.editDimension, 'String') );
        a = getappdata(MainFig, 'a');
        
        L = zeros(N, N, 2);
        L(:,:,1) = sqrt(alpha1/T)*a; % Relaxation (dissipation) operator
        L(:,:,2) = sqrt(alpha2/T)*(1/(hbar*w))*H; % Dephasing operator
        
        DissipType = 'dephasing with dissipation';
    case 5 % Thermalization
        N = str2double( get(handles.editDimension, 'String') );
        a = getappdata(MainFig, 'a');
        %         figure;
        %         spy(a)
        ad = a';
        Temperature = str2double(get(handles.editTemperature, 'String')); 
        kB = getappdata(MainFig, 'kB')
        kBT = kB*Temperature;
        L = zeros(N, N, 2);
        
        L(:,:,1) = sqrt(alpha1/T)*a;
        L(:,:,2) = sqrt(alpha1/T)*exp(- 0.5* hbar*w/kBT )*ad;
        
        DissipType = 'thermalization';

    otherwise        
        disp('Error: the dissipative evolution is not yet written.')
end

if get(handles.popupSelectDissipation, 'Value') > 1
    
    rho_v_0 = rho0(:);
    tmax = T*SimulationLength;
    
    statusMessage = ['Calculating density matrix per Lindblad equation (', ...
        DissipType, ') ...'];
    textStatusMsg(handles.textStatus, statusMessage, [1 0 0]);
    ODETic = tic;
    [ t, rho_time ] = ode45(@LindbladN, [0 tmax], rho_v_0, [], H, L, hbar );
    ODETime = toc(ODETic);
    statusMessage = ['Density matrix calculation complete. (', ...
        num2str(ODETime), ' s)'];
    textStatusMsg(handles.textStatus, statusMessage, [0 0 0]);
    
    ParameterTic = tic;
    statusMessage = ['Calculating parameters of interest from ', ...
        'density matrix ...'];
    textStatusMsg(handles.textStatus, statusMessage, 'k');
    
    nt = length(t);
    
    Qdistrib = zeros(N, nt);
    Edistrib = zeros(N, nt);
    SVN = zeros(1, nt);
    Eexp = zeros(1, nt);
    Qexp = zeros(1, nt);
    
    for tind = 1:nt
        rhot = reshape( rho_time(tind,:), [N N] );
        
        Edistrib(:, tind) = real( diag(rhot) );
        Qdistrib(:, tind) = real( diag( UQE*rhot*(UQE') ) );
        SVN(tind) = entropyVonNeumann( rhot ); % real( trace(rhot) );
        Eexp(tind) = real(trace(H*rhot));
        Qexp(tind) = real(trace(Q*rhot));

        
    end % END [ for tidx = 1:nt ]
    
    ParameterTime = toc(ParameterTic);
    statusMessage = ['Done calculating parameters of interest (', ...
        num2str(ParameterTime), ' s)'];
    textStatusMsg(handles.textStatus, statusMessage, 'k');
    
    disp('Dissipative time evolution is calculated.')
    
    set(handles.textTimeIdxMax, 'String', num2str(nt) );
    
    UpdateEditboxDiscreteFromSlider( handles.sliderTimeIdx, ...
        handles.editTimeIdx, handles.textTimeIdxMin, ...
        handles.textTimeIdxMax);
end % END [ if get(handles.popupSelectDissipation, 'Value') > 1 ]

set(handles.textTimeIdxMax, 'String', num2str(nt));

setappdata(MainFig, 'Qdistrib', Qdistrib);
setappdata(MainFig, 'Edistrib', Edistrib);
setappdata(MainFig, 'SVN', SVN);
setappdata(MainFig, 't', t);
setappdata(MainFig, 'Eexp', Eexp);
setappdata(MainFig, 'Qexp', Qexp);

CalcTime = toc(CalcStart);
FinishedReport = ['Calculated time evolution in ', num2str(CalcTime), ...
    ' s.'];
textStatusMsg(handles.textStatus, FinishedReport, 'k');

% close(waitfig);
figure(MainFig);

set(handles.sliderTimeIdx, 'Value', 0);
UpdateEditboxDiscreteFromSlider( handles.sliderTimeIdx, handles.editTimeIdx, ...
    handles.textTimeIdxMin, handles.textTimeIdxMax )
UpdatePlots(handles)

% ------------------------------------------------------------------------
function UpdatePlots(handles)

UpdateStaticPlot(handles);
UpdateDynamicPlot(handles);

% ------------------------------------------------------------------------
function UpdateStaticPlot(handles)

MainFig = gcf;
hbar = getappdata(MainFig, 'hbar');
FontName = getappdata(MainFig, 'FontName');
FontSizePlotLabel = getappdata(MainFig, 'FontSizePlotLabel');


% Update handles.axes1 (static axes)
switch get(handles.selectAxes1Plot, 'Value')
    case 1 % Potential with QHO Energy Ladder
        plotPotentialAndQHOLadder( handles.axes1 );
    case 2 % Time-varying position
        plotPositionVsTime( handles, handles.axes1 );
        
    case 3 % Trace of the RDM
        % plotTraceDensityMatrix( handles, handles.axes1 );
        plotEntropy( handles, handles.axes1 );
    case 4 % <E>
        plotEnergyTrace( handles, handles.axes1 );
        
        
end


% ------------------------------------------------------------------------
function UpdateDynamicPlot(handles)

MainFig = gcf;
hbar = getappdata(MainFig, 'hbar');
w = 2*pi*getappdata(MainFig, 'f');
FontName = getappdata(MainFig, 'FontName');
FontSizePlotLabel = getappdata(MainFig, 'FontSizePlotLabel');
N = str2double(get(handles.editDimension, 'String'));

Qdistrib = getappdata(MainFig, 'Qdistrib');
Edistrib = getappdata(MainFig, 'Edistrib');
Qval = getappdata(MainFig, 'Qval');  % [nm]
Eval = getappdata(MainFig, 'Eval');  % [nm]
Qpoints = linspace(Qval(1), Qval(end), 200);  % [nm]
Qlim = [Qval(1), Qval(end)];
Qexp = getappdata(MainFig, 'Qexp');  % [nm]
Qmax = max(Qval);

tidx = str2double(get(handles.editTimeIdx, 'String'));

% Plot the position distribution in handles.axes2
axes(handles.axes2);
bar(Qval*10, Qdistrib(:,tidx)); grid on;
set(handles.axes2, 'FontName', FontName, 'FontSize', FontSizePlotLabel, ...
    'GridLineStyle', '-');
maxProb = max(max(Qdistrib));
ylim([0 1.05*maxProb]);
xlbl = xlabel('$Q$ ($\AA$)'); set(xlbl, 'Interpreter', 'latex');
ylbl = ylabel('Probability');

% Plot the energy distribution in handles.axes3
axes(handles.axes3);
barh(Eval/(hbar*w), Edistrib(:,tidx)); grid on;
Eexpval = mean(Edistrib(:,tidx));
set(handles.axes3, 'FontName', FontName, 'FontSize', FontSizePlotLabel, ...
    'GridLineStyle', '-');
maxProb = max(max(Edistrib));
xlim([0 1.05*maxProb]);
ylim([0 round(N/2)]);
% hold on;
% EexpTrace = plot([0 1.05], (Eexpval/(hbar*w))*[1 1], 'LineWidth', 2, ...
%     'Color', [1 0 0 ]);
% hold off;
ylbl = ylabel('$E/\hbar \omega$');
set(ylbl, 'Interpreter', 'latex');
xlbl = xlabel('Probability');
% drawnow;

UpdateCursor( handles )
 

axes(handles.axes4);
drawDoubleDot_v2( handles.axes4, 0.5*[1 1], Qexp(tidx)/Qmax, ...
    'ChargeColor', 'w');

% ------------------------------------------------------------------------
function UpdateCursor( handles )

MainFig = gcf;

TimeCursor = getappdata(MainFig, 'TimeCursor');
if isempty(TimeCursor)
   JobType = 'drawCursor';  
else
    if ~ishandle(TimeCursor)
        JobType = 'drawCursor';
    else
        JobType = 'moveCursor';
    end
end

% JobType

if get( handles.selectAxes1Plot, 'Value' ) ~= 1
    T = 1/getappdata(MainFig, 'f'); % SHO period (S);
    Qval = getappdata(MainFig, 'Qval');  % [nm]
    SimulationLength = str2double(get(handles.editDuration, 'String'));
    t = getappdata(MainFig, 't');
    Eexp = getappdata(MainFig, 'Eexp');  % [eV]
    hbar = getappdata(MainFig, 'hbar');  % [eV*s]
    w =2*pi/T;

    tidx = str2double(get(handles.editTimeIdx, 'String'));
    SVN = getappdata(MainFig, 'SVN');
    
end



axes(handles.axes1);

switch JobType
    case 'drawCursor'
        
        hold on;

        switch get( handles.selectAxes1Plot, 'Value' )
            case 1 % Potential - do nothing
                
            case 2 % Position - draw the vertical white line
                
                
                TimeCursor = plot( (t(tidx)/T)*[1 1],  ... 
                    10*[min(Qval) max(Qval)], 'Color', [1 1 1], ...
                    'LineWidth',  2, 'LineStyle', '--');
                
            case 3 % Entropy - draw a red circle
                
                % Add a cursor
                
                St = SVN(tidx);
                tT = t(tidx)/T;
                TimeCursor = plot(tT, St);
                set(TimeCursor, 'Color', [1 0 0], 'LineWidth', 3, 'LineStyle', 'none', ...
                    'Marker', 'o');
                
            case 4 % Energy - draw a red circle
                Et = Eexp(tidx);
                tT = t(tidx)/T;
                TimeCursor = plot(tT, Et/(hbar*w));
                set(TimeCursor, 'Color', [1 0 0], 'LineWidth', 3, 'LineStyle', 'none', ...
                    'Marker', 'o');
        end
        
        hold off;

        setappdata(MainFig, 'TimeCursor', TimeCursor);
    case 'moveCursor'
        
        switch get( handles.selectAxes1Plot, 'Value' )
            case 1 % Potential - do nothing
                
            case 2 % Position - move the vertical white line
                St = SVN(tidx);
                tT = t(tidx)/T;
                set(TimeCursor, 'XData', (t(tidx)/T)*[1 1], ...
                    'YData', 10*[min(Qval) max(Qval)] );
                
            case 3 % Entropy - move a red circle
                % Add a cursor
                St = SVN(tidx);
                tT = t(tidx)/T;
                
                set(TimeCursor, 'XData', (t(tidx)/T), 'YData', St );
                
            case 4 % Energy - move a red circle
                Et = Eexp(tidx)/(hbar*w);
                tT = t(tidx)/T;
                set(TimeCursor, 'XData', tT, 'YData', Et );
        end
        
end
   
% ------------------------------------------------------------------------
function plotPositionVsTime( handles, targetAxes )

MainFig = gcf;

FontName = getappdata(MainFig, 'FontName');
FontSizePlotLabel = getappdata(MainFig, 'FontSizePlotLabel');
hbar = getappdata(MainFig, 'hbar'); % [eV*s] Reduced Plank constant
UQE = getappdata(MainFig, 'UQE'); % Transformation matrix: energy basis to position
f = getappdata(MainFig, 'f');
T = 1/getappdata(MainFig, 'f'); % SHO period (S);
Qval = getappdata(MainFig, 'Qval');  % [nm]
SimulationLength = str2double(get(handles.editDuration, 'String'));

switch get(handles.popupSelectDissipation, 'Value')
    case 1 % No dissipation
        % DISSIPATION IS DISABLED
        
        TimeSamples = str2double( get(handles.editSamples, 'String') );
        t = linspace(0, SimulationLength*T, SimulationLength*TimeSamples);

    otherwise
        t = getappdata(MainFig, 't');
end


axes(targetAxes);

Qdistrib = getappdata(MainFig, 'Qdistrib');

% sztT = size(t/T)
% szQval = size(Qval*10)
% szQdist = size(real(Qdistrib))
pcolor(t/T, Qval*10, real(Qdistrib));
shading interp
lighting phong
set(targetAxes, 'FontName', FontName, 'FontSize', FontSizePlotLabel);
title('Probability');

Qlim = [Qval(1) Qval(end)];
ylim(Qlim*10) % Angstroms 
xlim([0 SimulationLength]) % Unitless
ylbl = ylabel('$Q$ ($\AA$)');
set(ylbl, 'Interpreter', 'latex');

xlbl = xlabel('$t/\tau$');
set(xlbl, 'Interpreter', 'latex');


% ------------------------------------------------------------------------
function plotPositionInitial( handles, targetAxes )

MainFig = gcf;

FontName = getappdata(MainFig, 'FontName');
FontSizePlotLabel = getappdata(MainFig, 'FontSizePlotLabel');
hbar = getappdata(MainFig, 'hbar'); % [eV*s] Reduced Plank constant
UQE = getappdata(MainFig, 'UQE'); % Transformation matrix: energy basis to position

axes(targetAxes);

Qval = getappdata(MainFig, 'Qval');  % [nm]

psi0 = getappdata(MainFig, 'psi0');

psi0selection = get(handles.popupSelectInitialState, 'Value');
rho0 = outerprod(psi0(:,psi0selection), psi0(:,psi0selection));

Qdistro = diag( UQE*rho0*(UQE'));

bar(Qdistro);

% ------------------------------------------------------------------------
function plotSHOPotential( targetAxes )

MainFig = gcf;

FontName = getappdata(MainFig, 'FontName');
FontSizePlotLabel = getappdata(MainFig, 'FontSizePlotLabel');
hbar = getappdata(MainFig, 'hbar'); % [eV*s] Reduced Plank constant

H = getappdata(MainFig, 'H'); % [eV]
[psi_stationary, E_QHO] = eig(H);

Q = getappdata(MainFig, 'Q'); % [nm]
Qval = getappdata(MainFig, 'Qval');  % [nm]
Qpoints = linspace(Qval(1), Qval(end), 200);  % [nm]
Qlim = [Qval(1), Qval(end)];

m = getappdata(MainFig, 'm');  % QHO mass. [eV*s^2/nm^2]
w = 2*pi*getappdata(MainFig, 'f'); % QHO angular frequency [rad/s]

VQ = 0.5* m * w^2 * Qpoints.^2;   % [eV]

axes(targetAxes);

Vplot = plot(Qpoints*10, VQ/(hbar*w), 'LineWidth', 2); % Unitless vs Angstroms
set(targetAxes, 'FontName', FontName, 'FontSize', FontSizePlotLabel, ...
    'GridLineStyle', '-');
xlim(Qlim*10) % Angstroms 
ylim([0 8]) % Unitless
ylbl = ylabel('$E/\hbar \omega$');
set(ylbl, 'Interpreter', 'latex');

xlbl = xlabel('$Q$ ($\AA$)');
set(xlbl, 'Interpreter', 'latex');

% ------------------------------------------------------------------------
function plotEntropy( handles, targetAxes )

MainFig = gcf;

FontName = getappdata(MainFig, 'FontName');
FontSizePlotLabel = getappdata(MainFig, 'FontSizePlotLabel');
T = 1/getappdata(MainFig, 'f'); % SHO period (S);
Qval = getappdata(MainFig, 'Qval');  % [nm]
SimulationLength = str2double(get(handles.editDuration, 'String'));
t = getappdata(MainFig, 't');

SVN = getappdata(MainFig, 'SVN');
axes(targetAxes);
plot(t/T, SVN, 'LineWidth', 2); grid on;
set(targetAxes, 'FontName', FontName, 'FontSize', FontSizePlotLabel, ...
    'GridLineStyle', '-');
% ylim( [ 0 1.05] )
xlbl = xlabel('$t/\tau$');
set(xlbl, 'Interpreter', 'latex');
ylbl = ylabel('$S$ (bits)');
set(ylbl, 'Interpreter', 'latex');

ylim([0 max(1.05*[1, max(SVN)])]);
xlim([0 t(end)/T]);

% % Add a cursor
% tidx = str2double(get(handles.editTimeIdx, 'String'));
% St = SVN(tidx);
% tT = t(tidx)/T;
% hold on;
% SVNCursor = plot(tT, St);
% set(SVNCursor, 'Color', [1 0 0], 'LineWidth', 3, 'LineStyle', 'none', ...
%     'Marker', 'o');
% hold off;
% ------------------------------------------------------------------------
function plotEnergyTrace( handles, targetAxes )

MainFig = gcf;

FontName = getappdata(MainFig, 'FontName');
FontSizePlotLabel = getappdata(MainFig, 'FontSizePlotLabel');
T = 1/getappdata(MainFig, 'f'); % SHO period (S);
Eexp = getappdata(MainFig, 'Eexp');  % [eV]
hbar = getappdata(MainFig, 'hbar');  % [eV*s]
w =2*pi/T;

t = getappdata(MainFig, 't');

%szEexp = size(Eexp)
%szt = size(t)
axes(targetAxes);
Etrace = plot(t/T, Eexp/(hbar*w), 'LineWidth', 2);
grid on;
set(targetAxes, 'FontName', FontName, 'FontSize', FontSizePlotLabel, ...
    'GridLineStyle', '-');
ylim([ 0 1.125*max(Eexp/(hbar*w))] )
xlbl = xlabel('$t/\tau$');
set(xlbl, 'Interpreter', 'latex');
ylbl = ylabel('$\left< E \right>/\hbar \omega$');
set(ylbl, 'Interpreter', 'latex');

% ------------------------------------------------------------------------
function plotPotentialAndQHOLadder( targetAxes )

MainFig = gcf;

FontName = getappdata(MainFig, 'FontName');
FontSizePlotLabel = getappdata(MainFig, 'FontSizePlotLabel');
hbar = getappdata(MainFig, 'hbar'); % [eV*s] Reduced Plank constant

H = getappdata(MainFig, 'H'); % [eV]
[psi_stationary, E_QHO] = eig(H);

Q = getappdata(MainFig, 'Q'); % [nm]
Qval = getappdata(MainFig, 'Qval');  % [nm]
Qpoints = linspace(Qval(1), Qval(end), 200);  % [nm]
Qlim = [Qval(1), Qval(end)];

m = getappdata(MainFig, 'm');  % QHO mass. [eV*s^2/nm^2]
w = 2*pi*getappdata(MainFig, 'f'); % QHO angular frequency [rad/s]

VQ = 0.5* m * w^2 * Qpoints.^2;   % [eV]

axes(targetAxes);

Vplot = plot(Qpoints*10, VQ/(hbar*w), 'LineWidth', 2); % Unitless vs Angstroms
hold on;
Evals = diag(E_QHO)*[1 1]/(hbar*w);
Eplot = plot(Qlim*10, Evals); % Unitless vs Angstroms
hold off;
set(Eplot, 'Color', 'm', 'LineWidth', 2); grid on;
set(targetAxes, 'FontName', FontName, 'FontSize', FontSizePlotLabel, ...
    'GridLineStyle', '-', 'YGrid', 'off');
xlim(Qlim*10) % Angstroms 
ylim([0 8]) % Unitless
ylbl = ylabel('$E/\hbar \omega$');
set(ylbl, 'Interpreter', 'latex');

xlbl = xlabel('$Q$ ($\AA$)');
set(xlbl, 'Interpreter', 'latex');


% ========================================================================
% ========================^^^ UTILITIES ^^^===============================
% ========================================================================


function editDimension_Callback(hObject, eventdata, handles)
% hObject    handle to editDimension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDimension as text
%        str2double(get(hObject,'String')) returns contents of editDimension as a double

MainFig = gcf;
N = str2double(get(hObject,'String'));
projCombination = normalizevector(rand([N,1])*2 - 1)
setappdata(MainFig, 'ProjComb', projCombination);

UpdateQHO( handles );

% --- Executes during object creation, after setting all properties.
function editDimension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDimension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in selectAxes1Plot.
function selectAxes1Plot_Callback(hObject, eventdata, handles)
% hObject    handle to selectAxes1Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selectAxes1Plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectAxes1Plot

cla(handles.axes1);
UpdatePlots(handles)

% --- Executes during object creation, after setting all properties.
function selectAxes1Plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectAxes1Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editDuration_Callback(hObject, eventdata, handles)
% hObject    handle to editDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDuration as text
%        str2double(get(hObject,'String')) returns contents of editDuration as a double
Periods = str2double(get(hObject,'String'));
Samples = str2double(get(handles.editSamples,'String'));

switch get(handles.popupSelectDissipation, 'Value')
    case 1 % No dissipation
        nt = round(Periods*Samples);
        set(handles.textTimeIdxMax, 'String', num2str( nt ) );
        UpdateEditboxDiscreteFromSlider( handles.sliderTimeIdx, ...
            handles.editTimeIdx, handles.textTimeIdxMin, ...
            handles.textTimeIdxMax );
    otherwise
        disp('Error: callback for the dissipative case for the editDuration button must be written.')
        
end

UpdateQHO( handles );

% --- Executes during object creation, after setting all properties.
function editDuration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSamples_Callback(hObject, eventdata, handles)
% hObject    handle to editSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSamples as text
%        str2double(get(hObject,'String')) returns contents of editSamples as a double
Periods = str2double(get(handles.editDuration,'String'));
Samples = str2double(get(hObject,'String'));

switch get(handles.popupSelectDissipation, 'Value')
    case 1 % No dissipation
        nt = round(Periods*Samples);
        set(handles.textTimeIdxMax, 'String', num2str( nt ) );
        UpdateEditboxDiscreteFromSlider( handles.sliderTimeIdx, ...
            handles.editTimeIdx, handles.textTimeIdxMin, ...
            handles.textTimeIdxMax );
    otherwise
        disp('Error: callback for the dissipative case for the editDuration button must be written.')
        
end

UpdateQHO( handles );

% --- Executes during object creation, after setting all properties.
function editSamples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxEnableDissipation.
function checkboxEnableDissipation_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxEnableDissipation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxEnableDissipation
CalculateEvolution( handles )

% --- Executes on selection change in popupSelectInitialState.
function popupSelectInitialState_Callback(hObject, eventdata, handles)
% hObject    handle to popupSelectInitialState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupSelectInitialState contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupSelectInitialState

CalculateEvolution( handles );
% plotPositionInitial( handles, handles.axes2 );

% --- Executes during object creation, after setting all properties.
function popupSelectInitialState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupSelectInitialState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderTimeIdx_Callback(hObject, eventdata, handles)
% hObject    handle to sliderTimeIdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

UpdateEditboxDiscreteFromSlider(hObject, handles.editTimeIdx, ...
    handles.textTimeIdxMin, handles.textTimeIdxMax)

UpdateDynamicPlot(handles)

% --- Executes during object creation, after setting all properties.
function sliderTimeIdx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderTimeIdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editTimeIdx_Callback(hObject, eventdata, handles)
% hObject    handle to editTimeIdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTimeIdx as text
%        str2double(get(hObject,'String')) returns contents of editTimeIdx as a double
MaxIdx = str2double( get(handles.textTimeIdxMax, 'String') );

IdxEntry = str2double( get(hObject,'String') );

if IdxEntry > MaxIdx
    NewIdx = MaxIdx;
else
    NewIdx = round(IdxEntry);
end
set(hObject, 'String', num2str(NewIdx));

UpdateSlider( handles.sliderTimeIdx, handles.editTimeIdx, handles.textTimeIdxMin, ...
    handles.textTimeIdxMax );

UpdateDynamicPlot(handles)


% --- Executes during object creation, after setting all properties.
function editTimeIdx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTimeIdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in selectAxes2Plot.
function selectAxes2Plot_Callback(hObject, eventdata, handles)
% hObject    handle to selectAxes2Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selectAxes2Plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectAxes2Plot

UpdateDynamicPlot(handles)

% --- Executes during object creation, after setting all properties.
function selectAxes2Plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectAxes2Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editAlpha1_Callback(hObject, eventdata, handles)
% hObject    handle to editAlpha1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editAlpha1 as text
%        str2double(get(hObject,'String')) returns contents of editAlpha1 as a double
CalculateEvolution(handles)

% --- Executes during object creation, after setting all properties.
function editAlpha1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAlpha1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTemperature_Callback(hObject, eventdata, handles)
% hObject    handle to editTemperature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTemperature as text
%        str2double(get(hObject,'String')) returns contents of editTemperature as a double
CalculateEvolution(handles)

% --- Executes during object creation, after setting all properties.
function editTemperature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTemperature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupSelectDissipation.
function popupSelectDissipation_Callback(hObject, eventdata, handles)
% hObject    handle to popupSelectDissipation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupSelectDissipation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupSelectDissipation

FontName = 'Times';
FontSizeLabel = 18;
BlankString = ' ';

switch get(hObject, 'Value')
    case 1 % none
        set(handles.editSamples, 'Enable', 'on', 'Visible', 'on');
        set(handles.textTimeSampling, 'Visible', 'on');
        
        set(handles.editAlpha1, 'Visible', 'off');
        set(handles.editAlpha2, 'Visible', 'off');
        set(handles.editTemperature, 'Visible', 'off');
        set(handles.popupSelectDephasing, 'Visible', 'off');
        
        set(handles.textTempLabel, 'Visible', 'off');
        
        convertAxes2Label(handles.axesT1Label, BlankString, ...
            'FontSize', FontSizeLabel, 'FontName', FontName);
        
        convertAxes2Label(handles.axesT2Label, BlankString, ...
            'FontSize', FontSizeLabel, 'FontName', FontName);
        
        
    case 2 % relaxation
        set(handles.editSamples, 'Enable', 'on', 'Visible', 'on');
        set(handles.textTimeSampling, 'Visible', 'off');
        
        set(handles.editAlpha1, 'Visible', 'on');
        set(handles.editAlpha2, 'Visible', 'off');
        set(handles.editTemperature, 'Visible', 'off');
        set(handles.popupSelectDephasing, 'Visible', 'off');
        
        
        set(handles.textTempLabel, 'Visible', 'off');

        convertAxes2Label(handles.axesT1Label, '$\tau/T_1$', ...
            'FontSize', FontSizeLabel, 'FontName', FontName);
        
        convertAxes2Label(handles.axesT2Label, BlankString, ...
            'FontSize', FontSizeLabel, 'FontName', FontName);
                
       
    case 3 % dephasing only
        set(handles.editSamples, 'Enable', 'on', 'Visible', 'on');
        set(handles.textTimeSampling, 'Visible', 'off');

        set(handles.editSamples, 'Visible', 'off');
        set(handles.editAlpha1, 'Visible', 'off');
        set(handles.editAlpha2, 'Visible', 'on', 'Enable', 'on');
        set(handles.editTemperature, 'Visible', 'off');
        set(handles.popupSelectDephasing, 'Enable', 'on', ...
            'Visible', 'off');
        
                
        set(handles.textTempLabel, 'Visible', 'off');

        convertAxes2Label(handles.axesT1Label, BlankString, ...
            'FontSize', FontSizeLabel, 'FontName', FontName);
        
        convertAxes2Label(handles.axesT2Label, '$\tau/T_2$', ...
            'FontSize', FontSizeLabel, 'FontName', FontName);
        
    case 4 % dephasing with relaxation
        set(handles.editSamples, 'Enable', 'on', 'Visible', 'on');
        set(handles.textTimeSampling, 'Visible', 'off');

        set(handles.editSamples, 'Visible', 'off');
        set(handles.editAlpha1, 'Visible', 'on');
        set(handles.editAlpha2, 'Visible', 'on');
        set(handles.editTemperature, 'Visible', 'off');
        set(handles.popupSelectDephasing, 'Enable', 'on', ...
            'Visible', 'off');

        set(handles.textTempLabel, 'Visible', 'off');
        
        convertAxes2Label(handles.axesT1Label, '$\tau/T_1$', ...
            'FontSize', FontSizeLabel, 'FontName', FontName);
        
        convertAxes2Label(handles.axesT2Label, '$\tau/T_2$', ...
            'FontSize', FontSizeLabel, 'FontName', FontName);
        
    case 5 % thermalization
        set(handles.editSamples, 'Enable', 'on', 'Visible', 'on');
        set(handles.textTimeSampling, 'Visible', 'off');
        
        set(handles.editSamples, 'Visible', 'off');
        set(handles.editAlpha1, 'Visible', 'on');
        set(handles.editAlpha2, 'Visible', 'off');
        set(handles.editTemperature, 'Visible', 'on', 'Enable', 'on');
        set(handles.popupSelectDephasing, 'Enable', 'off', ...
            'Visible', 'off');
        
        set(handles.textTempLabel, 'Visible', 'on');
        
        convertAxes2Label(handles.axesT1Label, '$\tau/T_1$', ...
            'FontSize', FontSizeLabel, 'FontName', FontName);
        
        convertAxes2Label(handles.axesT2Label, BlankString, ...
            'FontSize', FontSizeLabel, 'FontName', FontName);

        
    otherwise
        set(handles.editSamples, 'Visible', 'off');
        
end

CalculateEvolution(handles)

% --- Executes during object creation, after setting all properties.
function popupSelectDissipation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupSelectDissipation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushAnimate.
function pushAnimate_Callback(hObject, eventdata, handles)
% hObject    handle to pushAnimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MainFig = gcf;
hbar = getappdata(MainFig, 'hbar');

set(hObject, 'Enable', 'off');
set(handles.pbStop, 'Enable', 'on');

t = getappdata(MainFig, 't');
Edistrib = getappdata(MainFig, 'Edistrib');
Qdistrib = getappdata(MainFig, 'Qdistrib');

T = 1/getappdata(MainFig, 'f');
nt = length(t);

ntPlay = 60; % number of points in animation
waittime = 0.005;

if nt > ntPlay
    dnsmpl = floor(nt/ntPlay);
else
    dnsmpl = 1;
    ntPlay = nt;
end

dns_vect = dnsmpl*(0:ntPlay-1) + 1;

td = t( dns_vect ); % downsampled time vector

ntd = length(td);
td_idx = 0;
StopLoop = 0;
LoopFinished = 0;
while (td_idx < ntd + 1) && ~StopLoop && ~LoopFinished
    
    StopLoop = strcmp( get(handles.pbStop, 'Enable'), 'off');
    
    td_idx = td_idx + 1;
    
    set(handles.editTimeIdx, 'String', num2str( dns_vect(td_idx) ));
    
    UpdateSlider(handles.sliderTimeIdx, handles.editTimeIdx, ...
        handles.textTimeIdxMin, handles.textTimeIdxMax);
    
    UpdateDynamicPlot( handles );
    
    pause(waittime);
    
    if td_idx == ntd
        td_idx = 0;
        LoopFinished = 1;
    end
    
end

set(handles.pbStop, 'Enable', 'off')
set(hObject, 'Enable', 'on')

% ln_dns = length(dns_vect)

% for tidx = 1:ntPlay
%     set(handles.editTimeIdx, 'String', num2str( dns_vect(tidx) ));
%     
%     UpdateSlider(handles.sliderTimeIdx, handles.editTimeIdx, ...
%         handles.textTimeIdxMin, handles.textTimeIdxMax);
%     
%     UpdateDynamicPlot( handles );
%     
%     pause(waittime);
% end % END [ for tidx = 1:ntPlay ]






% 
% 
% 
% if nt > ntPlay
%     idx_dns = dnsmpl*(0:numPoints-1) + 1;
%     
%     % ntr = length(idx_dns) % reduced number of time samples
%     td = t(dnsmpl*(0:numPoints-1) + 1); % downsampled time vector
% 
%     
%     [td, Edistribd] = DownsampleData(t, Edistrib, ntPlay);
%     [td, Qdistribd] = DownsampleData(t, Qdistrib, ntPlay);
% else
%     td = t;
%     Edistribd = Edistrib;
%     Qdistribd = Qdistrib;
% end
% 
% ntPlay = length(td);
% 
% for tidxp = 1:ntPlay
%     
% end % END [ for tidxp = 1:ntPlay ]



function editAlpha2_Callback(hObject, eventdata, handles)
% hObject    handle to editAlpha2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editAlpha2 as text
%        str2double(get(hObject,'String')) returns contents of editAlpha2 as a double
CalculateEvolution(handles)

% --- Executes during object creation, after setting all properties.
function editAlpha2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAlpha2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupSelectDephasing.
function popupSelectDephasing_Callback(hObject, eventdata, handles)
% hObject    handle to popupSelectDephasing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupSelectDephasing contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupSelectDephasing
CalculateEvolution(handles)


% --- Executes during object creation, after setting all properties.
function popupSelectDephasing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupSelectDephasing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFreq_Callback(hObject, eventdata, handles)
% hObject    handle to editFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFreq as text
%        str2double(get(hObject,'String')) returns contents of editFreq as a double
UpdateQHO( handles );


% --- Executes during object creation, after setting all properties.
function editFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMass_Callback(hObject, eventdata, handles)
% hObject    handle to editMass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMass as text
%        str2double(get(hObject,'String')) returns contents of editMass as a double
UpdateQHO( handles );


% --- Executes during object creation, after setting all properties.
function editMass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbStop.
function pbStop_Callback(hObject, eventdata, handles)
% hObject    handle to pbStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'Enable', 'off')
set(handles.pushAnimate, 'Enable', 'on');


% --------------------------------------------------------------------
function OptionsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to OptionsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function VisualizationMenu_Callback(hObject, eventdata, handles)
% hObject    handle to VisualizationMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function VisualizeEnergySubmenu_Callback(hObject, eventdata, handles)
% hObject    handle to VisualizeEnergySubmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function EnableAdvancedMenuMI_Callback(hObject, eventdata, handles)
% hObject    handle to EnableAdvancedMenuMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CheckedStatus = get(hObject, 'Checked');
AdvancedMenu = [handles.VisualizationMenu];

switch CheckedStatus
    case 'on'
        set(hObject, 'Checked', 'off');
        
        set(AdvancedMenu, 'Visible', 'off');
        
    case 'off'
        set(hObject, 'Checked', 'on');
        
        set(AdvancedMenu, 'Visible', 'on');
end

% --------------------------------------------------------------------
function VisualizeQHOEnergyMI_Callback(hObject, eventdata, handles)
% hObject    handle to VisualizeQHOEnergyMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MainFig = gcf;

hbar = getappdata(MainFig, 'hbar'); % [eV]
% N = str2double( get(handles.editDimension, 'String') ); % SHO dimension
% SimulationLength = str2double( get(handles.editDuration, 'String') ); % SHO periods
T = 1/getappdata(MainFig, 'f'); % QHO period [s]
H = getappdata(MainFig, 'H'); % [eV]
Q = getappdata(MainFig, 'Q'); % [nm]
UQE = getappdata(MainFig, 'UQE'); % Transformation from eigenbasis to position basis
alpha1 = str2double( get( handles.editAlpha1, 'String' ) ); % = \tau/T_1
alpha2 = str2double( get( handles.editAlpha2, 'String' ) );
f = getappdata(MainFig, 'f'); % [eV]
w = 2*pi*f;

t = getappdata(MainFig, 't');

Eexp = real(getappdata(MainFig, 'Eexp'));
E = Eexp - 0.5*hbar*w;

m = E/E(1);

NewFig = figure;
NewAx = axes;
plot(t/T, m, 'LineWidth', 2, 'Color', 'b'); grid on;
set(NewAx, 'GridLineStyle', '-', 'FontName', 'Times', 'FontSize', 24);
CouplingType = get(handles.popupSelectDissipation, 'Value');
if CouplingType == 2
    hold on;
    plot((1/alpha1)*[1 1], [0 1], [0 t(end)/T], exp(-1)*[1 1], ...
        'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
    t1label = text(1/alpha1, -0.075, '$\frac{T_1}{\tau}$', ...
        'Interpreter', 'latex', 'FontSize', 20, 'Color', 'r')
    expneg1label = text(t(end)/T, exp(-1), '$e^{-1}$', ...
        'Interpreter', 'latex', 'FontSize', 20, 'Color', 'r')
    hold off;
end
xlabel('$t/\tau$', 'Interpreter', 'latex');
ylabel('$f$', 'Interpreter', 'latex')

% GaussianCheck(t/T, m)


% --------------------------------------------------------------------
function VisualizeMultiFunctionPlotMI_Callback(hObject, eventdata, handles)
% hObject    handle to VisualizeMultiFunctionPlotMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MainFig = gcf;

NewFig = GetFigHandle( handles );
NewAx = axes;

figure(MainFig);

hbar = getappdata(MainFig, 'hbar');
FontName = getappdata(MainFig, 'FontName');
FontSizePlotLabel = getappdata(MainFig, 'FontSizePlotLabel');

% Update handles.axes1 (static axes)
switch get(handles.selectAxes1Plot, 'Value')
    case 1 % Potential with QHO Energy Ladder
        plotPotentialAndQHOLadder( NewAx );
        PlotType = 'V_and_E';
        
    case 2 % Time-varying position
        plotPositionVsTime( handles, NewAx );
        PlotType = 'ProbQ_v_t';
    case 3 % Trace of the RDM
        % plotTraceDensityMatrix( handles, handles.axes1 );
        plotEntropy( handles, NewAx );
        PlotType = 'SVN_v_t';
    case 4 % <E>
        plotEnergyTrace( handles, NewAx );
        PlotType = 'E_v_t';
        
end

figure(MainFig);
savename = [MakeSaveName(handles), '_', PlotType, '.eps'];

figure(NewFig);
PrintFigure(handles, MainFig, NewFig, savename);

% --------------------------------------------------------------------
function VisualizationOptions_Callback(hObject, eventdata, handles)
% hObject    handle to VisualizationOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function PrintFigsMI_Callback(hObject, eventdata, handles)
% hObject    handle to PrintFigsMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
toggleMenuItemCheck(hObject)

% --------------------------------------------------------------------
function ReuseFigsMI_Callback(hObject, eventdata, handles)
% hObject    handle to ReuseFigsMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
toggleMenuItemCheck(hObject)


% --------------------------------------------------------------------
function CloseFigsMI_Callback(hObject, eventdata, handles)
% hObject    handle to CloseFigsMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FigList = getappdata(gcf, 'FigList');

nFigs = length(FigList)
if nFigs > 0
    for FigIdx = nFigs:-1:1
        if ishandle(FigList(FigIdx))
            close(FigList(FigIdx))
            disp(['Closed handle ', num2str(FigList(FigIdx))]);
        end

        if length(FigList) > 1
            FigList = FigList(1:FigIdx-1);
        else
            FigList = [];
        end

    end % END: for FigIdx = 1:nFigs
end % END: if nFigs > 0

setappdata(gcf, 'FigList', FigList);


% --------------------------------------------------------------------
function MovieSubMenu_Callback(hObject, eventdata, handles)
% hObject    handle to MovieSubMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function MovieProbDistEnergyPositionMI_Callback(hObject, eventdata, handles)
% hObject    handle to MovieProbDistEnergyPositionMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MainFig = gcf;

Path = getappdata(MainFig, 'Path');
hbar = getappdata(MainFig, 'hbar');

t = getappdata(MainFig, 't');
Edistrib = getappdata(MainFig, 'Edistrib');
Qdistrib = getappdata(MainFig, 'Qdistrib');

Eexp = getappdata(MainFig, 'Eexp');
Qexp = getappdata(MainFig, 'Qexp');

maxPDQ = max(max(Qdistrib));

T = 1/getappdata(MainFig, 'f');
nt = length(t);

H = getappdata(MainFig, 'H');
Q = getappdata(MainFig, 'Q'); % [nm]
Qval = getappdata(MainFig, 'Qval');  % [nm]

Qpoints = linspace(Qval(1), Qval(end), 200);  % [nm]
Qlim = [Qval(1), Qval(end)];

m = getappdata(MainFig, 'm');  % QHO mass. [eV*s^2/nm^2]
w = 2*pi*getappdata(MainFig, 'f'); % QHO angular frequency [rad/s]

VQ = 0.5* m * w^2 * Qpoints.^2;   % [eV]

UQE = getappdata(MainFig, 'UQE'); % Transformation from eigenbasis to position basis
alpha1 = str2double( get( handles.editAlpha1, 'String' ) );
alpha2 = str2double( get( handles.editAlpha2, 'String' ) );
f = getappdata(MainFig, 'f'); % [eV]
w = 2*pi*f;

[V, E_stat] = eig(H);      % V are stationary states
es = diag(E_stat);         % Energy of stationary states
Vinv = inv(V);

savename = [MakeSaveName(handles), '_Movie_PD_EQ.avi'];
save_tgt = fullfile(Path.data, savename);

CaptureFrames = str2num(get(handles.editCaptureFrames, 'String'));

MovieFrames = textboxval( handles.editNumMovFrames );
FrameRate = textboxval( handles.editFrameRate );

DissipQHOMakeMoviePD_E_Q( t/T, Edistrib, Qdistrib, 10*Qval, ...
    10*Qpoints, VQ/(hbar*w), Eexp/(hbar*w), Qexp, 'MovieName', save_tgt, ...
    'Frames', MovieFrames, 'FPS', FrameRate, 'CaptureFrames', CaptureFrames);


function editNumMovFrames_Callback(hObject, eventdata, handles)
% hObject    handle to editNumMovFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumMovFrames as text
%        str2double(get(hObject,'String')) returns contents of editNumMovFrames as a double


% --- Executes during object creation, after setting all properties.
function editNumMovFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumMovFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFrameRate_Callback(hObject, eventdata, handles)
% hObject    handle to editFrameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFrameRate as text
%        str2double(get(hObject,'String')) returns contents of editFrameRate as a double


% --- Executes during object creation, after setting all properties.
function editFrameRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFrameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCaptureFrames_Callback(hObject, eventdata, handles)
% hObject    handle to editCaptureFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCaptureFrames as text
%        str2double(get(hObject,'String')) returns contents of editCaptureFrames as a double


% --- Executes during object creation, after setting all properties.
function editCaptureFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCaptureFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
