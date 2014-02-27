%minimum variance short selling portfolio
function [PortRisk,PortWts,mineigenvalue]=minvarport(ExpCovariance,ConSet)
azioni=length(ExpCovariance);
%METHOD MATLAB
if (nargin<3 | isempty(ConSet) )
   % A default constraint matrix ConSet will be created if this is not entered.
   ConSet = portcons('default', azioni);
end
[Aineq, Bineq, Aeq, Beq, LB, UB] = ineqparse(ConSet);
% Set the options: ('LargeScale' mode os turned off because it can cause 
%                   some warnings to be thrown.)
localmaxiter = size(Bineq,1) + size(LB,1);      % <- change
localmaxiter = max(localmaxiter,azioni);
localmaxiter = 10 * max(localmaxiter,20);

options = optimset('display', 'off', 'largescale', 'off', 'MaxIter', localmaxiter);
W0 = ones(azioni, 1)/azioni;
F = zeros(azioni, 1);

PortWts= quadprog(ExpCovariance, F, Aineq, Bineq, Aeq, Beq, LB, UB, W0, options);
PortWts=PortWts';
PortRisk=sqrt(PortWts(1,:)*ExpCovariance*PortWts(1,:)');
mineigenvalue=min(min(eig(ExpCovariance)));
%----------------------------------------------------------------
% Auxiliary function(s)
%----------------------------------------------------------------

function [A,b,Aeq,beq,LB,UB] = ineqparse(Ain, bin)
%INEQPARSE Find inequalities, equalities, and bounds implied by Ain*x <= bin.
%  Identifies equalities specified as Arow*x <= bval, -Arow*x <= -bval.
%  Parses out duplicate entries and all zero equations.
%  Finds bound constraints among inequalities.
%
%  [A,b,Aeq,beq,LB,UB] = ineqparse(Ain, bin)
%  [A,b,Aeq,beq,LB,UB] = ineqparse([Ain, bin])
%
%  The function does not catch linear combinations of inequalities which
%  together imply an equality constraint.
%
%  See also PORTCONS.
%
%----------------------------------------------------------------------
%
% % Test the logic to parse out inequalites with single-entry equations
% % 1 : redundant equaltiy
% % 2 : equality
% % 3 : redundant upper bound
% % 4 : lower bound
% % 5 : upper bound
% Astart = [1 2 -1 1 3 -1 -2 -4 5 3]'
% Ain = full(sparse(1:length(Astart),abs(Astart),Astart))
% [A,b,Aeq,beq,LB,UB] = ineqparse(Ain, zeros(length(Astart),1))
%
% % Catch rows which are multiples
% m = 1:length(Astart)
% Ain = diag(m)*Ain
% [A,b,Aeq,beq,LB,UB] = ineqparse(Ain, zeros(length(Astart),1))
%
% % Degenerate case with equality, lower, upper bounds
% C = portcons('default',3,'AssetLims',0,[0.5 0.6 0.7],3)
% [A,b,Aeq,beq,LB,UB] = ineqparse(C)
%
% % Case with a general inequality constraint
% C = portcons('default',3,'AssetLims',0,[0.5 0.6 0.7],3, ...
%              'Custom',[0.1 0.2 0.3],0.40)
% [A,b,Aeq,beq,LB,UB] = ineqparse(C)
%

% J. Akao 8/22/99

% find usage ineqparse(ConSet)
if nargin==1
  bin = Ain(:,end);
  Ain = Ain(:,1:end-1);
end

[NumEquations, NumVars] = size(Ain);
if any( size(bin)~=[NumEquations, 1] )
  error('finance:portopt:mismatchAandB','dimensions of A and b are inconsistent')
end

% Pull out degenerate rows
I = all(Ain==0,2);
if(any(I))
   warning('Degenerate rows found in constraint matrix. Eliminating these contraints');
   Ain(I,:) = [];
   bin(I) = [];
end

% Constraint rows
ConRows = [Ain, bin];

% Form numerator and denominator dot products.
%
% row I and row J are the same direction when:
%   rowI*rowJ' == sqrt(rowI*rowI')*sqrt(rowJ*rowJ') 
%        numIJ == denIJ
%
% row I and row J are the opposite direction when:
%   rowI*rowJ' == - sqrt(rowI*rowI')*sqrt(rowJ*rowJ') 
%        numIJ == - denIJ

% square (rowI*rowJ') but keep the sign
numIJsqrt = ConRows*ConRows';
numIJ = sign(numIJsqrt).*(numIJsqrt.*numIJsqrt);

% form (rowI*rowI') times (rowJ*rowJ')
rowKdot = dot(ConRows, ConRows, 2);
[rowIdot, rowJdot] = meshgrid(rowKdot, rowKdot);
denIJ = rowIdot .* rowJdot;

% record which equations are negations or duplicates
% denIJ is always positive
% take the upper triangular part only
%
% isdupIJ [NumEqs x NumEqs] row I is a positive multiple of row J
% isnegIJ [NumEqs x NumEqs] row I is a negative multiple of row J
reltol = 1000*eps;
isdupIJ = ( denIJ*(1-reltol) <  numIJ ) & (  numIJ < denIJ*(1+reltol) );
isnegIJ = ( denIJ*(1-reltol) < -numIJ ) & ( -numIJ < denIJ*(1+reltol) );

isdupIJ = triu(isdupIJ, 1);
isnegIJ = triu(isnegIJ, 1);

% search through the equations and clean out equalities and duplicates.
% store the equalities separately.
%
% ConEqs  [NumEqs   x NumVars+1] : [Aeq, beq]
% ConRows [NumInEqs x NumVars+1] : [A, b]
ConEqs = zeros(0, NumVars+1);

i=1;
while (i < size(ConRows,1) )
  % find negations and duplicates of this row
  RowIsNeg = isnegIJ(i,:);
  RowIsDup = isdupIJ(i,:);
  
  % negations and duplicates should be removed from the inequality list
  IndRemove = RowIsNeg | RowIsDup;
  
  if any(RowIsNeg)
    % add the row to the equality list
    ConEqs = [ConEqs; ConRows(i,:)];
    
    % remove the row from the inequality list along with negs and dups
    IndRemove(i) = 1;
  else
    % equation i has been left in
    i = i + 1;
  end
  
  % remove equations from the inequality list
  ConRows(IndRemove,:) = [];
  isnegIJ = isnegIJ(~IndRemove, ~IndRemove);
  isdupIJ = isdupIJ(~IndRemove, ~IndRemove);
end

% Break up into left and right hand sides
Aeq = ConEqs(:,1:NumVars);
beq = ConEqs(:,NumVars+1);
A = ConRows(:,1:NumVars);
b = ConRows(:,NumVars+1);

% search through the inequalities and find bounds
% SingleValue * x(Ind) <= b(Ind)
%
% IndSingle   [NumInEqs x 1] true if only 1 non-zero value in row of A
% SingleValue [NumInEqs x 1] only valid for IndSingle == 1
%
% VarNum       [NumInEqs x NumVars] column of each entry of A
% SingleVarNum [NumInEqs x 1] column of first non-zero entry in A
%
IndSingle   = sum(A~=0 , 2) == 1;
SingleValue = sum(A    , 2);

IndLower = IndSingle & ( SingleValue < 0 );
IndUpper = IndSingle & ( SingleValue > 0 );

VarNum = (1:NumVars);
VarNum = VarNum(ones(size(A,1),1),:);
VarNum(A==0) = Inf;
SingleVarNum = min(VarNum,[],2);

if any(IndLower)
  LB = -Inf*ones(NumVars,1);

  % find the variable and the bound value
  VarNum = SingleVarNum(IndLower);
  BVal = b(IndLower)./SingleValue(IndLower);

  % apply the most restrictive bound to each variable
  UniqVar = unique(VarNum);
  if length(UniqVar)==length(VarNum)
    % no variables have multiple parallel bounds
    LB(SingleVarNum(IndLower)) = BVal;
  else
    for Var=UniqVar(:)'
      LB(Var) = max( BVal( VarNum==Var ) );
    end
  end
    
else
  LB = [];
end

if any(IndUpper)
  UB = Inf*ones(NumVars,1);
  
  % find the variable and the bound value
  VarNum = SingleVarNum(IndUpper);
  BVal = b(IndUpper)./SingleValue(IndUpper);

  % apply the most restrictive bound to each variable
  UniqVar = unique(VarNum);
  if length(UniqVar)==length(VarNum)
    % no variables have multiple parallel bounds
    UB(SingleVarNum(IndUpper)) = BVal;
  else
    for Var=UniqVar(:)'
      UB(Var) = min( BVal( VarNum==Var ) );
    end
  end
    
else
  UB = [];
end

% remove lower or upper bound inequalities
A(IndSingle,:) = [];
b(IndSingle) = [];