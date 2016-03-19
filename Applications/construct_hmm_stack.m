function construct_hmm_stack(recordInfo)
%
% example
%     construct_hmm_stack('recordSummary.txt')
%

%% Load data information
[name, b1, b2] = textread(recordInfo, '%s %f %f');

for i = 1 : length(name)
    inputNAME{i} = ['../data/' name{i}];
    input = load(['data/' name{i}]);
    b3(i,1) = 1;
    b4(i,1) = size(input,1);
end

beginNend = [b3 b4 b1 b2];
dataL = max(b2);

save('recordSummary_info.mat', 'inputNAME', 'beginNend')


%% Run the profile-HMM algorithm

cd code

% Criteria of convergence
iterTol = 0.1;  %% Iteration stops if log likelihood does not increase larger than iterTol%
iterMax = 10;   %% Iteration stops if the number of iteration is larget than iterMax.

contiIter = true;
run_for_communication('newStack', dataL, date, 1, 'step', 0);
i = 0;
while contiIter == true
    i = i + 1;
    
    for inputI = 1 : length(name)
        call_run_for_each_input(inputI, dataL)
    end
    LL(i) = run_for_communication('newStack', dataL, date, 1, 'step', 0);

    if i ~= 1
        if abs(LL(i)-LL(i-1))/LL(i-1) < iterTol/100
            contiIter = false;
        end
    end
    if i > iterMax
        contiIter = false;
    end
end

cd ..

