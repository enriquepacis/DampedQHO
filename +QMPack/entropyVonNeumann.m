function varargout = entropyVonNeumann( varargin )
% ENTROPYVONNEUMANN calculates the Von Neumann entropy of a density matrix
%      RHO.
%
% SVN = entropyVonNeumann ( RHO ) calculates the Von Neuman entropy
%      SVN for RHO.
%
% This is based on Equation (2.347) on page 89 of REF [1].
%
%   REFERENCES
%   [1] Mahler G. and V.A. Weberus, QUANTUM NETWORKS: DYNAMICS OF OPEN
%       NANOSTRUCTURES (Second Edition), Springer, Berlin, 1998. 
%
%
% By E. P. Blair
%    University of Notre Dame
%    152033R FEB 2012
%

rho = varargin{1};

verbose = 0;

ev_raw = eig(rho);
ev_real = real(ev_raw); % remove imaginary artifacts of calculation
max_imag = max(abs(imag(ev_raw))); % find the biggest imaginary artifact removed
min_ev = min(ev_real);  
all_ev = abs(ev_real);
ev = all_ev(all_ev~=0); % remove those eigenvalues which are precisely zero
% these eigenvalues caues problems: log2(0) = Inf,
% and 0*Inf = NaN ... causing entropy to be NaN.

if verbose
    if min_ev < -1E-9
        disp(['entropyVonNeumann.m: Warning - negative probabilities as low as ', ...
            num2str(min_ev),' were removed.']);
    end
    
    if max_imag > 1E-12
        disp(['entropyVonNeumann.m: Warning - imaginary probabilities as high as +/-', ...
            num2str(max_imag), ' were removed.']);
        rho
        ev_raw_T = ev_raw'
    end
end % END [ if verbose ... ]

szrho = size(rho);
dimrho = log2(szrho(1));

SVN = real(-sum(ev.*log2(ev)));

if isnan(SVN)
    disp('entropyVonNeumann.m: Warning - calcuated entropy is NaN.')
    if dimrho < 5
        rho
    end
    ev
    max(imag(ev))
    min(ev)
    max(ev)
    
    l2ev = log2(ev)
    
    [tf zero_ev_idx] = ismember(0, ev);
    zero_ev_idx
    nzeroev= length(zero_ev_idx);
    disp(['There are ' num2str(nzeroev), ' zero eigenvalues.'])
    for k = 1: nzeroev
        kqx = l2ev(zero_ev_idx(k))
        skqx = ev(zero_ev_idx(k))*kqx
    end
    
elseif isinf(SVN)
        disp('entropyVonNeumann.m: Warning - calcuated entropy is INF.')
    
end
    

varargout{1} = SVN;

