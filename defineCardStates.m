function [action, policy,Qvals,pcards,T,statevals] = defineCardStates

%early (lower number) states are characterized by more potential actions
%values learned in early states generalize to later, but without the change
%in the unavailable action

%rows are player's hand, columns are dealer's showing card
splitspace = zeros(10,10); %availableactions = Hit, Stand, Double, Split: pairs
firstroundspace = zeros(26,10); %availableactions = Hit, Stand, Double: card sums 5 and up, ace-nonace pairs
finalspace = zeros(26,10); %availableactions = Hit, Stand: card sums 5 and up, ace-2xnonace triplets
statevals = [splitspace;firstroundspace;finalspace];  %the current iteration doesn't use these, but no reason to delete them

%automatically generalize across similar sums
%generalize across similar states as well

for n = 1:10
    pcards{n} = [num2str(n) '-' num2str(n)];
    for j = 1:10
        policy{n,j} = [.25 .25 .25 .25];
        action{n,j} = 'HSDP';
        Qvals{n,j} = [0 0 0 0];
    end
end
for n = 11:19
    pcards{n} = ['A-' num2str(n-9)];
    for j = 1:10
        policy{n,j} = [1/3 1/3 1/3];
        action{n,j} = 'HSD';
        Qvals{n,j} = [0 0 0];
    end
end
for n = 20:36
    pcards{n} = num2str(n-15);
    for j = 1:10
        policy{n,j} = [1/3 1/3 1/3];
        action{n,j} = 'HSD';
        Qvals{n,j} = [0 0 0];
    end
end
for n = 37:45
    pcards{n} = ['A-' num2str(n-35)];
    for j = 1:10
        policy{n,j} = [.5 .5];
        action{n,j} = 'HS';
        Qvals{n,j} = [0 0];
    end
end
for n = 46:62
    pcards{n} = num2str(n - 41);
    for j = 1:10
        policy{n,j} = [.5 .5];
        action{n,j} = 'HS';
        Qvals{n,j} = [0 0];
    end
end
pcards{1} = 'A-A';

%True correct choices
T(1,:) = {'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P'};
T(2,:) = {'H' 'P' 'P' 'P' 'P' 'P' 'P' 'H' 'H' 'H'};
T(3,:) = {'H' 'P' 'P' 'P' 'P' 'P' 'P' 'H' 'H' 'H'};
T(4,:) = {'H' 'H' 'H' 'H' 'P' 'P' 'H' 'H' 'H' 'H'};
T(5,:) = {'H' 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'H'};
T(6,:) = {'H' 'P' 'P' 'P' 'P' 'P' 'H' 'H' 'H' 'H'};
T(7,:) = {'H' 'P' 'P' 'P' 'P' 'P' 'P' 'H' 'H' 'H'};
T(8,:) = {'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P'};
T(9,:) = {'S' 'P' 'P' 'P' 'P' 'P' 'S' 'P' 'P' 'S'};
T(10,:) = {'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S'};
%pairs above
T(11,:) = {'H' 'H' 'H' 'H' 'DH' 'DH' 'H' 'H' 'H' 'H'};
T(12,:) = {'H' 'H' 'H' 'H' 'DH' 'DH' 'H' 'H' 'H' 'H'};
T(13,:) = {'H' 'H' 'H' 'DH' 'DH' 'DH' 'H' 'H' 'H' 'H'};
T(14,:) = {'H' 'H' 'H' 'DH' 'DH' 'DH' 'H' 'H' 'H' 'H'};
T(15,:) = {'H' 'H' 'DH' 'DH' 'DH' 'DH' 'H' 'H' 'H' 'H'};
T(16,:) = {'H' 'DS' 'DS' 'DS' 'DS' 'DS' 'S' 'S' 'H' 'H'};
T(17,:) = {'S' 'S' 'S' 'S' 'S' 'DS' 'S' 'S' 'S' 'S'};
T(18,:) = {'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S'};
T(19,:) = {'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S'};
%aces above
T(20,:) = {'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H'};
T(21,:) = {'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H'};
T(22,:) = {'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H'};
T(23,:) = {'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H' 'H'};
T(24,:) = {'H' 'H' 'DH' 'DH' 'DH' 'DH' 'H' 'H' 'H' 'H'};
T(25,:) = {'H' 'DH' 'DH' 'DH' 'DH' 'DH' 'DH' 'DH' 'DH' 'H'};
T(26,:) = {'DH' 'DH' 'DH' 'DH' 'DH' 'DH' 'DH' 'DH' 'DH' 'DH'};
T(27,:) = {'H' 'H' 'H' 'S' 'S' 'S' 'H' 'H' 'H' 'H'};
T(28,:) = {'H' 'S' 'S' 'S' 'S' 'S' 'H' 'H' 'H' 'H'};
T(29,:) = {'H' 'S' 'S' 'S' 'S' 'S' 'H' 'H' 'H' 'H'};
T(30,:) = {'H' 'S' 'S' 'S' 'S' 'S' 'H' 'H' 'H' 'H'};
T(31,:) = {'H' 'S' 'S' 'S' 'S' 'S' 'H' 'H' 'H' 'H'};
T(32,:) = {'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S'};
T(33,:) = {'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S'};
T(34,:) = {'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S'};
T(35,:) = {'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S'};
T(36,:) = {'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S' 'S'};








