function Z = linkage_corr(Y, method)
%LINKAGE Create hierarchical cluster tree.
%   Z = LINKAGE(Y) creates a hierarchical cluster tree, using the single
%   linkage algorithm.  The input Y is a correlation matrix.  Y may also be
%   a more general similarity
%   matrix conforming to the output format of PDIST.
%
%   Z = LINKAGE(Y, METHOD) creates a hierarchical cluster tree using
%   the specified algorithm. The available methods are:
%
%      'single'    --- farthest similarity
%      'complete'  --- nearest similarity
%      'average'   --- unweighted average similarity 
%      'weighted'  --- weighted average similarity 
%      
%   Cluster information will be returned in the matrix Z with size m-1
%   by 3, where m is the number of observations in the original data.
%   Column 1 and 2 of Z contain cluster indices linked in pairs
%   to form a binary tree. The leaf nodes are numbered from 1 to
%   m. They are the singleton clusters from which all higher clusters
%   are built. Each newly-formed cluster, corresponding to Z(i,:), is
%   assigned the index m+i, where m is the total number of initial
%   leaves. Z(i,1:2) contains the indices of the two component
%   clusters which form cluster m+i. There are m-1 higher clusters
%   which correspond to the interior nodes of the output clustering
%   tree. Z(i,3) contains the corresponding linkage distances between
%   the two clusters which are merged in Z(i,:), e.g. if there are
%   total of 30 initial nodes, and at step 12, cluster 5 and cluster 7
%   are combined and their distance at this time is 1.5, then row 12
%   of Z will be (5,7,1.5). The newly formed cluster will have an
%   index 12+30=42. If cluster 42 shows up in a latter row, that means
%   this newly formed cluster is being combined again into some bigger
%   cluster.
%
%   The hausdorff method can produce a cluster tree that is
%   not monotonic.  This occurs when the similarity of the union of two
%   clusters to a third cluster is more than the similarity of either
%   individual cluster to that third cluster.  In such a case, sections of
%   the dendrogram change direction.  This is an indication that another
%   method should be used.
%  
%
%   See also PDIST, INCONSISTENT, COPHENET, DENDROGRAM, CLUSTER,
%   CLUSTERDATA, KMEANS, SILHOUETTE.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.18.4.11 $

% Check size of the input vector
[k, n] = size(Y);
m = ceil(sqrt(2*n)); % m = (1+sqrt(1+8*n))/2, but works for large n

% Selects appropiate method
if nargin == 1 % set default switch to be 's'
    method = 'si';
    s=1;
else
    okmethods = {'single','nearest',...
                 'complete','farthest',...
                 'average','upgma',...
                 'weighted','wpgma'};
    methodkeys = {'si','si','co','co','av','av','we','we'};
    s = strmatch(lower(method), okmethods);
    if isempty(s)
        error('stats:linkage:BadMethod','Unknown method name: %s.',method);
    elseif length(s)>1
        error('stats:linkage:BadMethod','Ambiguous method name: %s.',method);
    else
        method = methodkeys{s};
    end
end

% The recursive distance updates for these methods only make sense when Y
% contains Euclidean distances (which will be squared).

  
% old linkage function 
Z = linkageold(Y,method);


% check for monotonicity warning
if any(diff(Z(:,3))>0)
     warning('stats:linkage:NonMonotonicTree',...
             'Non-monotonic cluster tree.');
end

%%%%%%%%%%%%%%%%%%%%%%%%% OLD LINKAGE FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Z = linkageold(Y, method)
%LINKAGEOLD Create hierarchical cluster tree using only m code.

[k, n] = size(Y);
m = ceil(sqrt(2*n)); % (1+sqrt(1+8*n))/2, but works for large n
if isa(Y,'single')
   Z = zeros(m-1,3,'single'); % allocate the output matrix.
else
   Z = zeros(m-1,3); % allocate the output matrix.
end

% during updating clusters, cluster index is constantly changing, R is
% a index vector mapping the original index to the current (row, column)
% index in Y.  N denotes how many points are contained in each cluster.
N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; % since m is changing, we need to save m in n.
R = 1:n;

% Square the distances so updates are easier.  The cluster heights will be
% square-rooted back to the original scale after everything is done.
if ~isempty(strmatch(method,['ce';'me';'wa']))
   Y = Y .* Y;
end
for s = 1:(n-1)
   if strcmp(method,'av')
      p = (m-1):-1:2;
      I = zeros(m*(m-1)/2,1);
      I(cumsum([1 p])) = 1;
      I = cumsum(I);
      J = ones(m*(m-1)/2,1);
      J(cumsum(p)+1) = 2-p;
      J(1)=2;
      J = cumsum(J);
      W = N(R(I)).*N(R(J));
      [v, k] = max(Y./W);
   else
      [v, k] = max(Y);
   end

   i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
   j = k - (i-1)*(m-i/2)+i;

   Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A

   % Update Y. In order to vectorize the computation, we need to compute
   % all the indices corresponding to cluster i and j in Y, denoted by I
   % and J.
   I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables
   U = [I1 I2 I3];
   I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
   J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];

   switch method
   case 'si' % single linkage
      Y(I) = max(Y(I),Y(J));
   case 'co' % complete linkage
      Y(I) = min(Y(I),Y(J));
   case 'av' % average linkage
      Y(I) = Y(I) + Y(J);
   case 'we' % weighted average linkage
      Y(I) = (Y(I) + Y(J))/2;
   otherwise 
      error('method is not supported')
   end
   J = [J i*(m-(i+1)/2)-m+j];
   Y(J) = []; % no need for the cluster information about j.

   % update m, N, R
   m = m-1;
   N(n+s) = N(R(i)) + N(R(j));
   R(i) = n+s;
   R(j:(n-1))=R((j+1):n);
end

Z(:,[1 2])=sort(Z(:,[1 2]),2);