function drawDoubleDot_v2( varargin ) % (haxc,P,Q)
%   O--O
%
%   circles represent vibrational antisymmetric breathing mode
%   P=[P1, P2] 
%   Qn=Q/Q0  is in the range [-2:2]
%   Let radius be proportional to represented quantity
%   r(3) for the charge
%   R(2) for the vibrational mode
%
%P=[0.8, 0.05, 0.15];
%P=[0.05, 0.05, 0.8];
%Q= 2;
%
% v2 of this function includes a modification by Erik Blair which allows
% the charge to be drawn 

haxc = varargin{1};
P = varargin{2};
Q = varargin{3};

DotColor = 'w';
ChargeColor = 'r';

% FUNCTIONAL ARGUMENTS
args = varargin(4:end);
while length(args) >= 2
    prop = args{1};
    val = args{2};
    args = args(3:end);
    switch prop
        case 'DotColor'
            DotColor = val;
        case 'ChargeColor'
            ChargeColor = val;
    end % END [ switch prop ]
end % END [ while length(args) >= 2 ]

cla(haxc)
%% maxima and minima
rmin=0;
rmax=1;
Rmin=1.1;
Rmax=2.0;
Rmid=(Rmin+Rmax)/2;
deltaR=(Rmax-Rmin)/2;

%% find r and R
r=rmin + (rmax-rmin)*P;
R(1)=Rmid - (Q)*deltaR;
R(2)=Rmid + (Q)*deltaR;

% plot
XYmax=5;
xL=-1;
xR=+1;
axis(haxc,[-XYmax  XYmax -XYmax  XYmax]);
axis(haxc,'square');

Dx=2.5;
rc=0.5;
DrawCircle(haxc,-Dx,0, R(1),'k',1,DotColor);
DrawCircle(haxc,-Dx,0, r(1),ChargeColor,1,ChargeColor);


DrawCircle(haxc,+Dx,0, R(2),'k',1,DotColor);
DrawCircle(haxc,+Dx,0, r(2),ChargeColor,1,ChargeColor);

axis(haxc,'off')


