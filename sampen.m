% ----------------------------------------------------------------------- %
%                           H    Y    D    R    A                         %
% ----------------------------------------------------------------------- %
% Function 'sampen' computes the Sample Entropy of a given signal.        %
%                                                                         %
%   Input parameters:                                                     %
%       - signal:       Signal vector with dims. [1xN]                    %
%       - m:            Embedding dimension (m < N).                      %
%       - r:            Tolerance (percentage applied to the SD).         %
%       - dist_type:    (Optional) Distance type, specified by a string.  %
%                       Default value: 'chebychev'                        %
%                                                                         %
%   Output variables:                                                     %
%       - value:        SampEn value.                                     %
% ----------------------------------------------------------------------- %
%   Script information:                                                   %
%       - Version:      1.0.                                              %
%       - Author:       V. Martínez-Cagigal                               %
%       - Date:         21/09/2018                                        %
% ----------------------------------------------------------------------- %
%   References:                                                           %
%       [1]     Richman, J. S., & Moorman, J. R. (2000). Physiological    %
%               time-series analysis using approximate entropy and sample %
%               entropy. American Journal of Physiology-Heart and         %
%               Circulatory Physiology, 278(6), H2039-H2049.              %
% ----------------------------------------------------------------------- %
function value = sampen(signal, m, r, dist_type)

    % Error detection and defaults
    if nargin < 3, error('Not enough parameters.'); end
    if nargin < 4
        dist_type = 'chebychev';
    end
    if ~isvector(signal)
        error('The signal parameter must be a vector.');
    end
    if m > length(signal)
        error('Embedding dimension must be smaller than the signal length (m<N).');
    end

    % Useful parameters
    signal = signal(:)';
    N = length(signal);
    sigma = std(signal);
    N_m = N - m;

    % Calculul distantelor pentru m (fara pdist)
    d_m = zeros(N_m*(N_m-1)/2, 1);
    idx = 0;
    for ii = 1:N_m-1
        for jj = ii+1:N_m
            idx = idx + 1;
            d_m(idx) = max(abs(signal(ii:ii+m-1) - signal(jj:jj+m-1)));
        end
    end

    if isempty(d_m)
        value = Inf;
    else
        % Calculul distantelor pentru m+1
        d_m1 = zeros(N_m*(N_m-1)/2, 1);
        idx = 0;
        for ii = 1:N_m-1
            for jj = ii+1:N_m
                idx = idx + 1;
                d_m1(idx) = max(abs(signal(ii:ii+m) - signal(jj:jj+m)));
            end
        end

        % Compute A si B
        B = sum(d_m  <= r*sigma);
        A = sum(d_m1 <= r*sigma);

        % Sample Entropy value - Richman & Moorman 2000
        if B == 0
            value = NaN;
        else
            value = -log((A/B)*((N-m+1)/(N-m-1)));
        end
    end

    % Upper bound daca value este Inf
    if isinf(value)
        value = -log(2/((N-m-1)*(N-m)));
    end
end