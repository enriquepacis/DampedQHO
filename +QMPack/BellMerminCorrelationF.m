function BellC=BellMerminCorrelationF(rho)
%        Author: C.S. Lent
%       Returns value of Bell quantity PAB + PAC + PBC
%       for two entangled two-state systems described by
%       density matrix rho
%          note: the 2-system basis is written writtent 
%          the a1,a0 basis set
%             |1 1>, |1 0>, |0 1>, |0 0>  
%
%      A value of BellC < 1 is a Bell violation, signalling
%      non-classical correlations between the two systems.
%      There is no underlying classical description that could
%      produce a Bell violation.

%% set angle between three basis states
theta=60;   % gives maximal Bell violation for pure Bell state:
%         Psi=[1; 0; 0; 1]/sqrt(2); 
%         Rho=Psi*Psi';
%
%% Explanation:
%  
%  Classically, each property A,B,C would have definite values
%  of either 1 or 0 for each particle
%  so a  particle "state" would be descibed as e.g., 
%       (A,B,C)=(1,1,0)
%  
%   We don't know the state of the two particles. 
%   We do know that the two particles are identical, that is, they
%       have the same state.
%   That same state could be (A,B,C) with 8 probabilities
% 
%   P(0,0,0), P(0,0,1), P(0,1,0)
%   P(0,1,1), P(1,0,0), P(1,0,1)        
%   P(1,1,0), P(1,1,1)
%   
%   These would have 
% 1=P(0,0,0)+P(0,0,1)+P(0,1,0)+P(0,1,1)+P(1,0,0)+P(1,0,1)+P(1,1,0)+P(1,1,1)
% 
%   PAA = Psame(A,B) is the probability that one particle
%        has the same A value as the other particle's B value.
%       
%         Since the  particles are identical, PAA=1=PBB=PCC
%
%    PAB = Psame(A,B)=P(0,0,0)+P(0,0,1)+P(1,1,0)+P(1,1,1)
%    PAC = Psame(A,C)=P(0,0,0)+P(0,1,0)+P(1,0,1)+P(1,1,1)
%    PBC = Psame(B,C)=P(0,0,0)+P(1,0,0)+P(0,1,1)+P(1,1,1)
%
%  PAB+PAC+PBC = 
%   P(0,0,0)+P(0,0,1)+P(0,1,0)+P(0,1,1)+P(1,0,0)+P(1,0,1)+P(1,1,0)+P(1,1,1)
%      + 2*P(0,0,0)+2*P(1,1,1)
%   = 1 + 2*P(0,0,0)+2*P(1,1,1) > 1
%      
%   So if there is a fact-of-the-matter regarding the state
%        of each particle PAB+PAC+PBC > 1
%
%   For random uniformly-distributed classical pairs P=1.5
%
%   We calculate the same quantity quantum mechanically
%    for the system in state described by the (reduced) density
%    matrix rho
%   

%% define 3 separate basis sets, rotated by theta degrees
a0=[0; 1];
a1=[1; 0];

R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)];

b0=R*a0;
b1=R*a1;

c0=R*b0;
c1=R*b1;

%% define projection operators on each basis vector
%   (the "0" vector and the "1" vector)
%   These operators correspond to the observable A, B, and C
%   a0 is an eigenstate of A with eigenvalue 0
%   a1 is an eigenstate of A with eigenvalue 1
%   b0 is an eigenstate of B with eigenvalue 0
%      etc.

PA0=a0*a0' ;
PB0=b0*b0';
PC0=c0*c0';

PA1=a1*a1' ;
PB1=b1*b1';
PC1=c1*c1';
%% define operators that project both systems on
%  the same (0 or 1) basis vector
% these operators are observables 
%   e.g. PAB= measuring the same value (0 or 1) when
%    measuring A on system 1 and B on system 2

PsameAA=kron(PA0,PA0)+kron(PA1,PA1);
PsameBB=kron(PB0,PB0)+kron(PB1,PB1);
PsameCC=kron(PC0,PC0)+kron(PC1,PC1);

PsameAB=kron(PA0,PB0)+kron(PA1,PB1);
PsameAC=kron(PA0,PC0)+kron(PA1,PC1);
PsameBC=kron(PB0,PC0)+kron(PB1,PC1);

%% calculate probabilities of getting same (1 or 0)
%  
PAA=trace(rho*PsameAA);
PBB=trace(rho*PsameBB);
PCC=trace(rho*PsameCC);

PAB=trace(rho*PsameAB);
PAC=trace(rho*PsameAC);
PBC=trace(rho*PsameBC);

%% return Bell-Mermin Correlation function
BellC=PAB+PAC+PBC;

