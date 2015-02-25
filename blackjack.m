function [Qvalue,combact,check,action,playercardnames,dealercardnames] = blackjack(decknum,alpha,epsilon,N,reshuffle,H17)
%usage: [Qvalue,action,check,playercardnames,dealercardnames] = blackjack(decknum,alpha,epsilon,N,reshuffle,H17)
%updates as TD(0)
%uses epsilon-greedy policy, but initializes with random policy
%default reshuffle = .8, H17 = 1
%allows for splitting (up to four times per hand) and doubling options
%didn't do TD(1+) because transitions vary based on deck penetration
%
%output:
% one table of player vs. dealer cards and recommended action
% H = hit, S = stand, D = double down, P = Split
% 
% first table distinguishes between first and subsequent actions, second
% table combines them.  If the model choices are different for first and
% subsequent choices (i.e. Hit for one and Stand for the other) this is
% represented by a blank in the graph.  Alternatively, if the model
% recommends doubling down, this is shown by either a DS or DH, indiciating
% you should double down if possible, otherwise hit or stand.
%
% The final surface plot shows a comparison between the model predictions
% and the 'true' actions a player should take.  This 'true' table is
% derived from H17 rules, with a deck penetration of .8 and 6 decks in
% play.  There are two horizontal green lines on the graph, marking out
% sections.  The top section has all pair hands, the middle section has all
% ace hands, and the lower section has all general hands.  The value of the
% dealer's hand increases from left to right (just like in the combined
% table).
%
%
% example usage:
% [Qvalue,action,check,playercardnames,dealercardnames] = blackjack(6,.05,.1,50000,.8,1);

if nargin < 5
    H17 = 1;
end
if nargin < 4
    reshuffle = .8;
end

[action, policy, Qvalue, playercardnames,truecorrect] = defineCardStates; %set basic policy, define statevals

cardtot = 52*decknum; %figure out the total number of cards available
cardindex = cardtot/13; %find separation points between card vals
basicdeck = zeros(cardtot,1); %initialize deck

% put cards into deck
% 1 = Ace, 10 = 10,J,Q,K, 2:9 = 2:9
for n = 1:10
    basicdeck((n-1)*cardindex + (1:cardindex),:) = n;
end
for n = 11:13
    basicdeck((n-1)*cardindex + (1:cardindex),:) = 10;
end


decklimit = round(cardtot*reshuffle); %card count prior to reshuffling (may be used for future card counting app)
[cardcount,shuffledeck] = shuffle(cardtot,basicdeck); %new card positions in the shuffled deck
learningiter = 0; %counter for number of experienced trials


while learningiter < N
    if  cardcount > decklimit
        [cardcount,shuffledeck] = shuffle(cardtot,basicdeck);
    end
    actioncount = 0;
    splitcount = 0;  %this will be a running split count
    [playercards,dealercards,cardcount] = deal(shuffledeck,cardcount); %deal cards to player and dealer
    cards = {playercards dealercards};
    resolved = false;
    statesplayed = {};
    actionschosen = [];
    currentcardset = 1;
    cardsetplayed = [];
    subcardsets = {[]};
    splitorder = 0;  %this will keep track of new splits
    mastercardsets = [1 0 0 0 0 0];
    reward = [0 0 0 0 0 0];
    while ~resolved %specifically, while the first hand dealt is unresolved
        [resolved,cards,currPi,currAct,currState] = getstate(cards,policy,action);
        [selAct,cards,cardcount,resolved,actioncount,splitcount] = play(cards,currPi,currAct,shuffledeck,cardcount,resolved,actioncount,splitcount);
        statesplayed{end+1} = currState; %#ok
        if ~isempty(selAct)
            actionschosen(end+1) = selAct; %#ok
        end
        cardsetplayed(end+1) = currentcardset; %#ok
        if splitcount > splitorder %then some finagling with the split order needs to happen
            subcardsets{currentcardset}(end+1) = splitcount+2; %this keeps track of where the splits came from, so we can assign credit later
            splitorder = splitorder+1;
        end
    end
    if splitcount > 0 %do splits
        for currentcardset = [1 3:5] %set 6 cannot have subsets, so is not included (set 2 is the dealer)
            if mastercardsets(currentcardset) > 0 %if the current card set has any subsets
                for n = 1:numel(subcardsets{currentcardset}) %go through all the card sets that were generated from card set n
                    subcardset = subcardsets{currentcardset}(n); %get the current subset
                    subcardsets{subcardset} = []; %set the subsets for the subset of set n to []
                    resolved = false;
                    tempcards = cards;
                    tempcards{subcardset} = cards{1};  %adjust the cards var to make play work easier
                    tempcards{1} = cards{subcardset};
                    while ~resolved  %go until the subset is resolved
                        [resolved,tempcards,currPi,currAct,currState] = getstate(tempcards,policy,action);
                        [selAct,tempcards,cardcount,resolved,actioncount,splitcount] = play(tempcards,currPi,currAct,shuffledeck,cardcount,resolved,actioncount,splitcount);
                        statesplayed{end+1} = currState; %#ok
                        actionschosen(end+1) = selAct; %#ok
                        cardsetplayed(end+1) = subcardset; %#ok
                        if splitcount > splitorder %then some finagling with the split order needs to happen
                            subcardsets{subcardset}(end+1) = splitcount+2; %this keeps track of where the splits came from, so we can assign credit later
                            splitorder = splitorder+1;
                            mastercardsets(subcardset) = 1;
                        end
                    end
                    cards = tempcards;
                    cards{1} = tempcards{subcardset};
                    cards{subcardset} = tempcards{1};%return the cards var to its initial state, plus changes to the subset
                end
            end
        end
    end %split loop
    
    if ~isempty(actionschosen) %if dealer had blackjack, completely irrelevant to play choices
        
        %get player scores
        playerscore = [0 0 0 0 0 0 0];
        for n = [1 3:6]
            trials = find(cardsetplayed == n);
            if trials
                playerscore(n) = score(cards{n});
            end
        end
        
        if sum(playerscore) > 0 %if at least one hand did not bust
            [cards,cardcount] = dealerhit(cards,shuffledeck,cardcount,H17); %now that play has resolved, let the dealer hit (only if all hands did not bust!)
        end
        
        dealerscore = score(cards{2});
        
        
        for n = [1 3:6]
            trials = find(cardsetplayed == n);
            if trials
                lastact = actionschosen(trials(end));
                if playerscore(n) > dealerscore
                    if strcmp(lastact,'D')
                        reward(n) = 2;
                    else
                        reward(n) = 1;
                    end
                elseif playerscore(n) == 0 || playerscore(n) < dealerscore
                    if strcmp(lastact,'D')
                        reward(n) = -2;
                    else
                        reward(n) = -1;
                    end
                else %playerscore(n) == dealerscore
                    reward(n) = 0;
                end
            end %if trials
        end
        
        %now that all basic rewards are calculated, get net rewards for split
        %decisions
        splitreward = [0 0 0 0 0 0];
        for n = [6:-1:3 1]
            if mastercardsets(n) > 0 && splitorder > 0
                splitreward(n) = reward(n) + sum(splitreward(subcardsets{n}));
            else
                splitreward(n) = reward(n);
            end
        end
        [policy,Qvalue] = update(reward,splitreward,actionschosen,statesplayed,cardsetplayed,alpha,Qvalue,policy,epsilon);
        
        learningiter = learningiter+1;
    end % end if statement that checked for dealer blackjack, learning iterations do not progress in this case either
end

dealercardnames = {'Ace';'2';'3';'4';'5';'6';'7';'8';'9';'10'};

[check combact] = generateOutput(playercardnames,Qvalue,action,truecorrect,N);
end %main function

%======================================================================

function [cardcount shuffledeck] = shuffle(cardtot,basicdeck)
randval = rand(cardtot,1);
[~,deckpos] = sort(randval);
shuffledeck = basicdeck(deckpos);
cardcount = 0;
end %shuffle function

%======================================================================

function [playercards,dealercards,cardcount] = deal(shuffledeck,cardcount)
playercards = shuffledeck([1 3]+[cardcount cardcount],1);
dealercards = shuffledeck([2 4]+[cardcount cardcount],1);
cardcount = cardcount + 4;
end %deal function

%======================================================================

function [cards,cardcount]=dealerhit(cards,shuffledeck,cardcount,H17)

donehitting = false;
[~,soft,truecount] = score(cards{2});
if truecount >= 18
    donehitting = true;
elseif truecount == 17 && soft && ~H17
    donehitting = true;
end

while ~donehitting
    cards{2}(end+1) = shuffledeck(cardcount+1);
    cardcount = cardcount+1;
    [~,soft,truecount] = score(cards{2});
    if truecount >= 18
        donehitting = true;
    elseif truecount == 17 && soft && ~H17
        donehitting = true;
    end
end

end %dealerhit function

%======================================================================

function [cardscore,soft,truecount] = score(cardnums)

soft = false;
if sum(cardnums) >= 12
    cardscore = sum(cardnums);
elseif ~isempty(find(cardnums == 1,1))
    acepos = find(cardnums == 1,1);
    cardnums(acepos) = 11;
    cardscore = sum(cardnums);
    soft = true;
else
    cardscore = sum(cardnums);
end

truecount = cardscore;

if cardscore > 21  %for busting!
    cardscore = 0;
end

end %score function

%======================================================================

function [myact,cards,cardcount,resolved,actioncount,splitcount] = ...
    play(cards,currPi,currAct,shuffledeck,cardcount,resolved,actioncount,splitcount)

if ~resolved    %if dealer blackjack, no action possible
    myact = chooseact(currPi,currAct);
    if splitcount >= 4 && strcmp(myact,'P') %not going to allow more than 4 splits.  The nested ifs above are ugly enough.
        currPi = currPi(1:3)./sum(currPi(1:3));
        currAct = 'HSD';
        myact = chooseact(currPi, currAct);
    end
    switch myact
        case 'H'
            cards{1}(end+1) = shuffledeck(cardcount+1);
            cardcount = cardcount+1;
            playerscore = score(cards{1});
            if playerscore == 0
                resolved = true;
            end
            actioncount = actioncount+1;
        case 'S'
            resolved = true;
            actioncount = actioncount+1;
        case 'D'
            cards{1}(end+1) = shuffledeck(cardcount+1);
            cardcount = cardcount+1;
            resolved = true;
            actioncount = actioncount+1;
        case 'P'
            switch splitcount
                case 0
                    cards = {[cards{1}(1); shuffledeck(cardcount+1)] cards{2} [cards{1}(2); shuffledeck(cardcount+2)]};
                case 1
                    cards = {[cards{1}(1); shuffledeck(cardcount+1)] cards{2} cards{3} [cards{1}(2); shuffledeck(cardcount+2)]};
                case 2
                    cards = {[cards{1}(1); shuffledeck(cardcount+1)] cards{2} cards{3} cards{4} [cards{1}(2); shuffledeck(cardcount+2)]};
                case 3
                    cards = {[cards{1}(1); shuffledeck(cardcount+1)] cards{2} cards{3} cards{4} cards{5} [cards{1}(2); shuffledeck(cardcount+2)]};
            end
            cardcount = cardcount+2;
            splitcount = splitcount+1;
            if cards{1}(1) == 1 %Ace!  Can only split once, and take a single card
                resolved = true;
            end
    end %switch
else
    myact = []; %no action available - dealer blackjack
end %resolved if statement
end %play function

%========================================================================

function myact = chooseact(currPi,currAct)
pthresh = currPi(1);
for n = 2:numel(currPi)
    pthresh(n) = pthresh(n-1)+currPi(n);
end
actp = rand;
n = 1;
while n <= numel(currPi)
    if actp <= pthresh(n)
        myact = currAct(n);
        n = numel(currPi) + 1;
    end
    n = n+1;
end

end %chooseact function

%========================================================================

function [resolved,cards,currPi,currAct,currState] = getstate(cards,policy,action)
resolved = false;

%only allow 4 splits

%if dealer has blackjack ...
if (cards{2}(2) == 1 && cards{2}(1) == 10) || (cards{2}(1) == 1 && cards{2}(2) == 10) %dealer has blackjack, push for player blackjack, otherwise loss.
    if sum(cards{1} == 11) && (cards{1}(1) == 1 || cards{1}(1) == 10)
        %                 reward = 0; technically this is true, but since no action
        %                 was taken, there is no learning.  Same for the else
        %                 clause
        resolved = true;
    else
        %                 reward = -1;
        resolved = true;
    end
end

%figure out current state
if numel(cards{1})==2 %first round
    if cards{1}(1)==cards{1}(2) %split state
        pstate = cards{1}(1);
    elseif cards{1}(1) == 1 || cards{1}(2) == 1 %ace state
        pstate = max(cards{1}) + 9;
    else
        pstate = sum(cards{1}) + 15; %general state
    end
else %not first round, no split or double available
    if ~isempty(find(cards{1}==1,1)) && sum(cards{1})<=11 %multi-card ace state
        pstate = sum(cards{1})+35;
    else
        pstate = sum(cards{1})+41;
    end
end
dstate = cards{2}(2);

%get current policy, actions, and state
currState = [pstate,dstate];
currPi = policy{pstate,dstate};
currAct = action{pstate,dstate};
end %getstate function

%====================================================================

function [policy,Qvalue] = update(reward,splitreward,actionschosen,statesplayed,cardsetplayed,alpha,Qvalue,policy,epsilon)
for n = 1:numel(actionschosen)
    p = statesplayed{n}(1);
    d = statesplayed{n}(2);
    switch actionschosen(n)
        case 'H'
            R = reward(cardsetplayed(n));
            a = 1;
        case 'S'
            R = reward(cardsetplayed(n));
            a = 2;
        case 'D'
            R = reward(cardsetplayed(n));
            a = 3;
        case 'P'
            R = splitreward(cardsetplayed(n));
            a = 4;
    end
    
    Qvalue{p,d}(a) = Qvalue{p,d}(a) + alpha*(R-Qvalue{p,d}(a));
    totalopts = numel(Qvalue{p,d});
    eprob = epsilon/totalopts;
    [~,maxQ] = max(Qvalue{p,d});
    policy{p,d} = eprob.*ones(1,totalopts);
    policy{p,d}(maxQ) = policy{p,d}(maxQ) + 1 - epsilon;
    
    if p > 10 && a < 3 && p < 37 %then updates can apply to the not first round case as well
        policy{p,d}(maxQ) = .8 - epsilon + eprob;%give a small preference for 'D' selection so it can be learned better
        policy{p,d}(3) = policy{p,d}(3) + .2;
        p = p+26;
        Qvalue{p,d}(a) = Qvalue{p,d}(a) + alpha*(R-Qvalue{p,d}(a));
        totalopts = numel(Qvalue{p,d});
        eprob = epsilon/totalopts;
        [~,maxQ] = max(Qvalue{p,d});
        policy{p,d} = eprob.*ones(1,totalopts);
        policy{p,d}(maxQ) = policy{p,d}(maxQ) + 1 - epsilon; 
    elseif p > 36 %apply that information to first round case Q's
        p = p - 26;
        Qvalue{p,d}(a) = Qvalue{p,d}(a) + alpha*(R-Qvalue{p,d}(a));
        totalopts = numel(Qvalue{p,d});
        eprob = epsilon/totalopts;
        [~,maxQ] = max(Qvalue{p,d});
        policy{p,d} = eprob.*ones(1,totalopts);
        policy{p,d}(maxQ) = policy{p,d}(maxQ) + 1 - epsilon;
    end
end

end %update function

%========================================================================
function [check combact] = generateOutput(playercardnames,Qvalue,action,truecorrect,N)

act = cell(62,10);

%long table output suppressed.  With proper generalization, all suggestions
%are consistent.

% %output tabular results
% fprintf('SPLIT TABLE BY FIRST AND SUBSEQUENT ROUNDS\n\n')
% fprintf('\tAce\t2\t3\t4\t5\t6\t7\t8\t9\t10\t\n')
for n = 1:62
%     fprintf('%s\t',playercardnames{n})
    for j = 1:10
        [~,maxQ] = max(Qvalue{n,j});
        act{n,j} = action{n,j}(maxQ);
        if Qvalue{n,j}(1) == Qvalue{n,j}(2)
            act{n,j} = 'N';
        end
%         fprintf('%s\t',act{n,j})
    end
%     fprintf('\n')
end
% fprintf('\n\n\n')

combact = act(1:10,:);

for n = 11:36 %combine across first round/subsequent rounds (accounting for double option)
    for j=1:10
        same = strcmp(act{n,j},act{n+26,j});
        if ~same
            if strcmp(act{n,j},'N')
                combact{n,j} = act{n+26,j};
            elseif strcmp(act{n+26,j},'N')
                combact{n,j} = act{n,j};
            elseif strcmp(act{n,j},'D')
                combact{n,j} = [act{n,j} act{n+26,j}]; %i.e. double if possible, otherwise...
            else 
                combact{n,j} = ' ';
            end
        else
            combact{n,j} = act{n,j};
        end
    end
end

check = zeros(36,10);
fprintf('SUGGESTED ACTIONS FOR ALL HANDS\n\n')
fprintf('\tAce\t2\t3\t4\t5\t6\t7\t8\t9\t10\t\n')
for n = 1:36
    fprintf('%s\t',playercardnames{n})
    for j = 1:10
        fprintf('%s\t',combact{n,j})
        if strcmp(combact{n,j},truecorrect{n,j})
            check(n,j) = 1;
        end
    end
    fprintf('\n')
end
fprintf('\n\n\n')

pcorrect = 100*mean(mean(check));

figure;
imagesc(check)
colormap gray
colorbar
hold on
plot([-.5:10.5],ones(size([-.5:10.5]))*10.5,'g','LineWidth',2)
plot([-.5:10.5],ones(size([-.5:10.5]))*19.5,'g','LineWidth',2)
title(['N = ' num2str(N) ', pcorrect = ' num2str(pcorrect) '%'],'FontSize',22)
fprintf('The mean percentage correct is %4.2f%%\n',pcorrect)

end
