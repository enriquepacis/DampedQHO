function varargout = SUsparse(varargin)
% LAMBDA = SUsparse(N)
%   SUsparse returns N-by-N generator matrices LAMBDA for the Special Unitary
%      Group SU(N). There will be S = N^2-1 such generators, so LAMBDA is
%      a N-by-N-by-S array. The returs are Matlab sparse matrices
%
% LAMBDA = SU(N, V)
%      The matrix V is a N-by-N column-concatenated matrix consisting of
%      N N-by-1 orthonormal vectors spanning a N-dimentsional Hilbert
%      space. V is an optional argument, and if unspecified, is an identity
%      matrix.
%
% [LAMBDA, P] = SU(N, V)
%
% [LAMBDA, P, U] = SU(N, V)
%
% [LAMBDA, P, U, V] = SU(N, V)
%
% [LAMBDA, P, U, V, W, ew] = SU(N, V)
%
% REFERENCES
%  [1] Mahler G. and V.A. Weberus, QUANTUM NETWORKS: DYNAMICS OF OPEN
%      NANOSTRUCTURES (Second Edition), Springer, Berlin, 1998. 
%
% Created by
%    Erik Blair
%    University of Notre Dame
%    012025ROCT2011
%
time_SUstart = tic;
Nwait = 32;
Nsave = 0;
mustCalc = 1;
VerboseReporting = 0;

switch nargin
    case 0
        error('SUsparse.M: Too few input arguments.')
    case 1
        N = varargin{1};
        V = eye(N);
    case 2
        N = varargin{1};
        V = varargin{2};
    otherwise
        error('SUsparse.M: Too many input arguments.');
end % END:switch nargin

storefilename = ['SUsparsedata_N', num2str(N), '.mat'];

S = N^2 - 1;

if N >= Nsave
    % check to see if previously-calculated data file is already stored
    stored_file_data = dir('SUsparsedata_N*.mat');
    
    if length(stored_file_data) == 0 % No data has been saved
        if VerboseReporting
            disp('No SUsparsedata_N*.mat files have been saved.');
        end
        mustCalc = 1;
    else                             % Some data has been saved
        NameList = {};
        for k = 1:length(stored_file_data)
            NameList{k} = stored_file_data(k).name;
        end
        
        if ismember(storefilename, NameList)
            if VerboseReporting
                disp(['SUsparsedata_N', num2str(N), '.mat HAS been PREVIOUSLY SAVED.'])
            end
            mustCalc = 0;
        else
            if VerboseReporting
                disp(['SUsparsedata_N', num2str(N), '.mat has NOT been saved.'])
            end
            mustCalc = 1;
        end
    end % END: if length(stored_file_data) == 0  ... else ...
 
end % END: if N >= Nsave


if mustCalc
    % Create N^2 projection P_{j,k} Operators
    P = {};
    
    if N >= Nwait
        wbp = waitbar(0, 'SUsparse.M: Please wait. Calculating P ...');
    end
    
    for a = 1:N
        for b = 1:N
            P{a,b} = sparse(outerprod(V(:,a), V(:,b)));
        end % for b = 1:N
        if N >= Nwait
            waitbar(a/N, wbp)
        end
    end % for a = 1:N
    
    if N >= Nwait
        close(wbp)
    end
    
    
    % I wrote out some cases and determined that there will be the following
    % numbers of u, v, and w matrices:
    %     u: sum(1:N-1)
    %     v: sum(1:N-1)
    %     w: N-1
    R = sum(1:N-1); % the number of u or v matrices
    u = zeros(N,N,R);
    v = u;
    w = zeros(N,N,N-1);
    
    % Create the lambda matrices
    %    Populate the u and v matrices first
    uv_count = 1;
    if N >= Nwait
        wbuv = waitbar(0, 'SUsparse.M: Please wait. Calculating u and v ...');
    end
    for k = 2:N
        for j = 1:k-1
            u(:,:,uv_count) = P{j,k}+P{k,j};
            v(:,:,uv_count) = 1i*(P{j,k}-P{k,j});
            uv_count = uv_count+1;
        end % END: for j = 1:k-1
        
        if N >= Nwait
            waitbar(k/N, wbuv)
        end
    end % END: for k = 1:Nk
    
    if N >= Nwait
        close(wbuv)
    end
    
    if N >= Nwait
        wbw = waitbar(0, 'SU.M: Please wait. Calculating w ...');
    end
    
    %    Now populate the w matrices first
    for l = 1:N-1
        tempX = zeros(N,N,l+1);
        for m = 1:l
            tempX(:,:,m) = P{m,m};
        end
        tempX(:,:,l+1) = -l*P{l+1, l+1};
        w(:,:,l) = -sqrt(2/(l*(l+1)))*sum(tempX,3);
        
        
        if N >= Nwait
            waitbar(l/(N-1), wbw)
        end
    end % for l = 1:N-1
    
    if N >= Nwait
        close(wbw)
    end
    
    
    X = zeros(N,N,S);
    X(:,:,1:R) = u;
    X(:,:,R+1:2*R) = v;
    X(:,:,2*R+1:end) = w;
    
    LAMBDA = {};
    for k = 1:S
        LAMBDA{k} = sparse(X(:,:,k));
    end % END: for k = 1:S
    
    % Cellarrays to store the u, v, and w operators. Refer to Eq 2.79 in
    %   Reference [1] (page 46).
    U = {};     V = {};     W = {};
    for k = 1:R
        U{k} = u(:,:,k);
        V{k} = v(:,:,k);
    end
    for k = 1:N-1
        W{k} = w(:,:,k);
    end
    
    % A matrix to store the eigenvalues of w operators. Refer to Eq. 2.109 in
    %   Reference [1] (page 51).  Here, l is the subscript and denotes a row of
    %   matrix ew; nu is the superscript and dentoes a column of the matrix ew.
    %
    ew = zeros(N-1, N);
    for l = 1:N-1
        
        for nu = 1:N
            
            if nu <= l
                ew(l, nu) = -sqrt(2/(l*(l+1)));
            elseif nu == l+1
                ew(l, nu) = l*sqrt(2/(l*(l+1)));
            else
                ew(l, nu) = 0;
            end % END: if nu <= l ... elseif nu == l+1 ...
            
        end % END: for nu = 1:N
        
    end % END: for l = 1:N-1
    
    save(storefilename, 'LAMBDA', 'P', 'U', 'V', 'W', 'ew', '-v7.3');
    if VerboseReporting
        disp(['SUsparse.M: Saved data to ''', storefilename, '''']);
    end

    
else
    load(storefilename);
    if VerboseReporting
        disp(['SUsparse.M: Data loaded from ''', storefilename, '''']);
    end
end

switch nargout
    case 0
        varargout{1} = LAMBDA;
    case 1 % LAMBDA = SU(N)
        varargout{1} = LAMBDA;
    case 2 % [LAMBDA, P] = SU(N)
        varargout{1} = LAMBDA;
        varargout{2} = P;
    case 3 % [LAMBDA, P, U] = SU(N)
        varargout{1} = LAMBDA;
        varargout{2} = P;
        varargout{3} = U;
    case 4 % [LAMBDA, P, U] = SU(N)
        varargout{1} = LAMBDA;
        varargout{2} = P;
        varargout{3} = U;
        varargout{4} = V;
    case 5 % [LAMBDA, P, U, V, W] = SU(N)
        varargout{1} = LAMBDA;
        varargout{2} = P;
        varargout{3} = U;
        varargout{4} = V;
        varargout{5} = W;
        
    case 6 % [LAMBDA, P, U, V, W] = SU(N)
        varargout{1} = LAMBDA;
        varargout{2} = P;
        varargout{3} = U;
        varargout{4} = V;
        varargout{5} = W;
        varargout{6} = ew;
end % switch nargout
SUruntime = toc(time_SUstart);

if VerboseReporting
    disp(['SUsparse(', num2str(N), ') took ', sec2hrstr(SUruntime)]);
end