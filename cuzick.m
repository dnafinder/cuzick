function STATS = cuzick(x,varargin)
% CUZICK: Perform the Cuzick's test on trend.
% This function provides a Wilcoxon-type test for trend across a group of
% three or more independent random samples.
%
% Assumptions:
% - Data must be at least ordinal
% - Groups must be selected in a meaningful order i.e. ordered
% If you do not choose to enter your own group scores then scores are
% allocated uniformly (1 ... k) in order of selection of the k groups.
% The null hypothesis of no trend across the groups T will have mean E(T),
% variance Var(T) and the null hypothesis is tested using the normalised
% test statistic z.
% A logistic distribution is assumed for errors. Please note that this test
% is often more powerful than applying the Wilcoxon rank-sum /
% Mann-Whitney test pairwise between more than two groups of data.
% Cuzick J. A Wilcoxon-Type Test for Trend. Statistics in Medicine
% 1985;4:87-89.
%
% Syntax:
%     STATS = cuzick(x)
%     STATS = cuzick(x, score)
%     STATS = cuzick(x, score, 'Display', DISPLAY)
%
%     Inputs:
%           X - N-by-2 data matrix:
%               X(:,1) = observations
%               X(:,2) = integer group labels
%               Group labels must be consecutive integers 1,2,...,k without gaps.
%
%           SCORE - optional vector of numeric scores (real values) associated
%                   with each group, length k. If omitted or empty, default
%                   scores are 1,2,...,k, where k is the number of groups.
%                   These scores encode the ordering (and spacing) of groups
%                   for the trend test. Scores need not be integers, but
%                   must not all be equal.
%
%           'Display' - logical flag (true/false), default true.
%                      If true, prints the group summary and the Cuzick
%                      statistics. If false, no output is printed and only
%                      the STATS structure is returned.
%
%     Outputs:
%           STATS - structure with fields:
%               STATS.GroupTable  : table with Group, Score, Samples, Ranks_sum
%               STATS.L           : sum of score(i)*n_i across groups
%               STATS.T           : sum of score(i)*R_i across groups
%               STATS.E           : expected value of T under H0
%               STATS.Var         : variance of T under H0
%               STATS.z           : normalised test statistic
%               STATS.pvalue      : one-tailed p-value (right tail)
%               STATS.tail        : 'right'
%               STATS.TiesFactor  : ties adjustment factor (2*t from tiedrank)
%
%   Example:
% Mice were inoculated with cell lines, CMT 64 to 181, which had been
% selected for their increasing metastatic potential. The number of lung
% metastases found in each mouse after inoculation are quoted below:
%
%                                 Sample
%                   ---------------------------------
%                      64   167  170  175  181
%                   ---------------------------------
%                      0    0    2    0    2
%                      0    0    3    3    4
%                      1    5    6    5    6
%                      1    7    9    6    6
%                      2    8    10   10   6
%                      2    11   11   19   7
%                      4    13   11   56   18
%                      9    23   12   100  39    
%                           25   21   132  60
%                           97
%                   ---------------------------------
%
%       Data matrix must be:
%    d=[0 0 1 1 2 2 4 9 0 0 5 7 8 11 13 23 25 97 2 3 6 9 10 11 11 12 21 ...
%       0 3 5 6 10 19 56 100 132 2 4 6 6 6 7 18 39 60];
%    g=[ones(1,8) 2.*ones(1,10) 3.*ones(1,9) 4.*ones(1,9) 5.*ones(1,9)];
%    x=[d' g'];
%
%           Calling on Matlab the function: STATS = cuzick(x)
% (in this case, the groups are automatically scored from 1 to 5)
%
%           Answer is:
%
% CUZICK'S TEST FOR NON PARAMETRIC TREND ANALYSIS
% --------------------------------------------------------------------------------
%     Group    Score    Samples    Ranks_sum
%     _____    _____    _______    _________
% 
%     1        1         8            79    
%     2        2        10           256    
%     3        3         9           229    
%     4        4         9         246.5    
%     5        5         9         224.5    
% 
% Ties factor: 366
% --------------------------------------------------------------------------------
%  
% CUZICK'S STATISTICS
% --------------------------------------------------------------------------------
%      L       T        E       Var       z       one_tailed_p_values
%     ___    ______    ____    _____    ______    ___________________
% 
%     136    3386.5    3128    14973    2.1125    0.01732   
%                         
% With these data we are interested in a trend in one direction only,
% therefore, we can use a one sided test for trend. We have shown a
% statistically significant trend for increasing number of metastases across
% these malignant cell lines in this order.
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo.75@gmail.com
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008). Cuzick's test: A Wilcoxon-Type Test for Trend.
% https://github.com/dnafinder/cuzick

% Input Error handling
p = inputParser;
addRequired(p,'x',@(y) validateattributes(y,{'numeric'}, ...
    {'real','finite','nonnan','nonempty','ncols',2}));
addOptional(p,'score',[],@(s) isempty(s) || ...
    (isnumeric(s) && isvector(s) && all(isreal(s(:))) && ...
     all(isfinite(s(:))) && ~all(isnan(s(:)))));
addParameter(p,'Display',true, ...
    @(d) (islogical(d) || (isnumeric(d) && isscalar(d) && (d==0 || d==1))));
parse(p,x,varargin{:});
x       = p.Results.x;
score   = p.Results.score;
Display = logical(p.Results.Display);
clear p

% Check that group labels are integers and consecutive 1..k
assert(all(x(:,2) == fix(x(:,2))), ...
    'cuzick:InvalidGroupLabels', ...
    'All elements of column 2 must be whole numbers (integer group labels).');

groups = x(:,2);
ug     = unique(groups);
k      = numel(ug); % number of groups

assert(min(ug)==1 && max(ug)==k && k==max(groups), ...
    'cuzick:ConsecutiveGroups', ...
    'Group labels in column 2 must be consecutive integers from 1 to k without gaps.');

% Check/define score (numeric scores for each group)
if isempty(score)
    score = 1:k;
else
    score = score(:).'; % ensure row vector
    assert(numel(score)==k, ...
        'cuzick:InvalidScoreLength', ...
        'Length of score must match the number of groups (k).');
    % Ensure scores are not all equal (degenerate trend)
    assert(numel(unique(score)) > 1, ...
        'cuzick:ConstantScore', ...
        'score must not be constant: use at least two distinct values.');
end

tr = repmat('-',1,80); % divisor

% Group sizes
ni = crosstab(groups); % elements for each group
N  = sum(ni);          % total elements

% Ranks and ties from pooled data
[r,t] = tiedrank(x(:,1)); % r = ranks, t = ties adjustment factor

% Preallocate vectors
R = zeros(1,k); % sum of ranks per group
L = zeros(1,k); % sum of score(i)*n_i per group
T = zeros(1,k); % sum of score(i)*R_i per group

for I = 1:k
    R(I) = sum(r(groups==I)); % sum of ranks of each group
    L(I) = score(I) * ni(I);
    T(I) = score(I) * R(I);
end

GroupTable = table((1:k)', score', ni, R', ...
    'VariableNames',{'Group','Score','Samples','Ranks_sum'});

if Display
    disp('CUZICK''S TEST FOR NON PARAMETRIC TREND ANALYSIS')
    disp(tr)
    disp(GroupTable)
    if t > 0
        fprintf('Ties factor: %d\n', 2*t);
    end
    disp(tr); 
    disp(' ')
end

% For the null hypothesis of no trend across the groups;
% T will have mean E(T), variance Var(T) and the null hypothesis is tested
% using the normalised test statistic z.
Lval = sum(L);
Tval = sum(T);
Et   = Lval * (N+1) / 2; % mean of T under H0
Vart = ((N*sum(score.*L) - Lval^2) * (N+1) / 12) - t/6; % variance of T under H0
z    = abs(Tval - Et) / sqrt(Vart); % z statistic

% One-tailed p-value from standard normal approximation (right tail)
p = 1 - normcdf(z);

if Display
    disp('CUZICK''S STATISTICS')
    disp(tr)
    disp(table(Lval,Tval,Et,Vart,z,p, ...
        'VariableNames',{'L','T','E','Var','z','one_tailed_p_values'}))
end

% Build output structure
STATS = struct();
STATS.GroupTable = GroupTable;
STATS.L          = Lval;
STATS.T          = Tval;
STATS.E          = Et;
STATS.Var        = Vart;
STATS.z          = z;
STATS.pvalue     = p;
STATS.tail       = 'right';
STATS.TiesFactor = 2*t;
end
