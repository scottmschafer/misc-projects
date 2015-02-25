function stats = anova(data,betmodel,betnames,wimodel,winames)

%usage: stats = anova(data,betmodel,betnames,wimodel,winames)
%
% data is a cell array - first level is between, second level is within.
% If the numel of the first level = 1, a standard anova is performed.  The
% betmodel is a matrix of the codes for the between subjects comparisons (including
% intercept) and wimodel stores the within subjects codes for each subject
% as a cell array
%
% betmodel must always be specified, even if it is only the intercept



fprintf('Source\tb\tSS\tdf\tMS\tF\tPRE\tp\n')
BetErr = 0;
%%first, fit the between model
if size(betmodel,2) > 1
    N = numel(data);
    betdata = zeros(N,1);
    %between model is more than just the intercept
    for n = 1:N
        betdata(n) = mean(data{n});
    end
    [SSEA,SSEC,SSR,BA,df,F,PRE,pA,pC,p] = fitmodel(betdata,betmodel);
    %Begin ANOVA table output
    n = numel(betnames)+1;
    fprintf('%s\t\t%2.2f\t%2.0f\t%2.2f\t%2.2f\t%2.2f\t%2.3f\n','Model',SSR{n+1},df{n+1},SSR{n+1}/df{n+1},F{n+1},PRE{n+1},p{n+1})
    for n = 1:numel(betnames)
        fprintf('%s\t%2.2f\t%2.2f\t%2.0f\t%2.2f\t%2.2f\t%2.2f\t%2.3f\n',betnames{n},BA(n+1),SSR{n+1},df{n+1},SSR{n+1}/df{n+1},F{n+1},PRE{n+1},p{n+1})
    end
    fprintf('Error\t\t%2.2f\t%2.0f\t%2.2f\n',SSEA,N-pA,SSEA/(N-pA))
    fprintf('Total\t\t%2.2f\t%2.0f\t%2.2f\n',SSEC{end},N-1,SSEC{end}/(N-1))
    fprintf('_______________________________________________________\n')
    keyboard
    BetErr = SSEC{end};
    stats.bet.SSEA = SSEA;
    stats.bet.SSEC = SSEC;
    stats.bet.b = BA;
    stats.bet.df = df;
    stats.bet.F = F;
    stats.bet.PRE = PRE;
    stats.bet.pA = pA;
    stats.bet.pC = pC;
    stats.bet.p = p;
end

WithinErr = 0;
if nargin > 3
    %within subjects comparisons!
    %first, compute the within subjects contrasts
    widata = zeros(numel(data),size(wimodel{1},2));
    for n = 1:numel(data)
        for j = 1:size(wimodel{n},2)
            widata(n,j) = data{n}'*wimodel{n}(:,j);
        end
    end
    fprintf('Within\n')
    for j = 1:size(wimodel{n},2)
        N = numel(widata(:,j));
        [SSEA,SSEC,SSR,BA,df,F,PRE,pA,pC,p] = fitmodel(widata(:,j),betmodel);
        fprintf('%s\t%2.2f\t%2.2f\t%2.0f\t%2.2f\t%2.2f\t%2.2f\t%2.3f\n',winames{j},BA(1),SSR{1},df{1},SSR{1}/df{1},F{1},PRE{1},p{1})
        WithinErr = WithinErr + SSR{1};
        for n = 1:numel(betnames)
            fprintf('%s*%s \t%2.2f\t%2.2f\t%2.0f\t%2.2f\t%2.2f\t%2.2f\t%2.3f\n',winames{j},betnames{n},BA(n+1),SSR{n+1},df{n+1},SSR{n+1}/df{n+1},F{n+1},PRE{n+1},p{n+1})
            WithinErr = WithinErr + SSR{n+1};
        end
        fprintf('%s Err\t\t%2.2f\t%2.0f\t%2.2f\n',winames{j},SSEA,N-pA,SSEA/(N-pA))
        WithinErr = WithinErr + SSEA;
        stats.wit(j).name = winames{j};
        stats.wit(j).between_intnames = betnames;
        stats.wit(j).SSEA = SSEA;
        stats.wit(j).SSEC = SSEC;
        stats.wit(j).b = BA;
        stats.wit(j).df = df;
        stats.wit(j).F = F;
        stats.wit(j).PRE = PRE;
        stats.wit(j).pA = pA;
        stats.wit(j).pC = pC;
        stats.wit(j).p = p;
    end
    fprintf('Total\t\t%2.2f\n',WithinErr)
    fprintf('_______________________________________________________\n')
end
TotErr = WithinErr + BetErr;
fprintf('Total\t\t%2.2f\n\n\n',TotErr)

end

function [SSEA,SSEC,SSR,BA,df,F,PRE,pA,pC,p] = fitmodel(dat,model)
    n = numel(dat);
    [BA, ~, ErrA] = regress(dat,model);
    SSEA = sum(ErrA.^2);
    pA = size(model,2);
    numdat = numel(dat);
    for k = 1:pA
        if k == 1
            [BC{k},~,ErrC{k}] = regress(dat,model(:,2:end));
        elseif k < pA
            [BC{k},~,ErrC{k}] = regress(dat,[model(:,1:k-1) model(:,k+1:end)]);
        else
            [BC{k},~,ErrC{k}] = regress(dat,model(:,1:k-1));
        end
        SSEC{k} = sum(ErrC{k}.^2);
        pC{k} = pA - 1;
    end
    [BC{k+1},~,ErrC{k+1}] = regress(dat,model(:,1));    %I'm not sure if this part is correct
    SSEC{k+1} = sum(ErrC{k+1}.^2);                      %
    pC{k+1} = 1;                                        %
    for k = 1:numel(BC)
        SSR{k} = SSEC{k}-SSEA;
        PRE{k} = SSR{k}/SSEC{k};
        df{k} = pA-pC{k};
        F{k} = (PRE{k}/df{k})/((1-PRE{k})/(numdat-pA));  
        p{k} = fpval(F{k},df{k},n-pA);
    end
end


function p = fpval(x,df1,df2)
%FPVAL F distribution p-value function.
%   P = FPVAL(X,V1,V2) returns the upper tail of the F cumulative distribution
%   function with V1 and V2 degrees of freedom at the values in X.  If X is
%   the observed value of an F test statistic, then P is its p-value.
%
%   The size of P is the common size of the input arguments.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also FCDF, FINV.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2010/11/08 02:37:46 $

if nargin < 3,
    error(message('stats:fpval:TooFewInputs'));
end

xunder = 1./max(0,x);
xunder(isnan(x)) = NaN;
p = fcdf(xunder,df2,df1);
end