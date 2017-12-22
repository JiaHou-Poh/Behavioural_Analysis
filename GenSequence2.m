function[order, check_lags]=GenSequence2(n,rep,min,max,num_seq)
%% ----------------- Script description -----------------------%%
% Generate sequence for item presentation based on set criterion for
% repetition intervals. Takes in the following inputs:
% 1) n: Number of items
% 2) rep : Number of repetition per item
% 3) min : Minimum interval between repetitions
% 4) max: Maximum interval between repetitions (mean will be (min+max)/2)
% 5) num_seq: Number of sequence required
%
% Coded base on Sumei's script.
% Take note that script requires equal occurence of all lag interval
% between the range of min & max. So N*Rep needs to be multiples of lag
% values.
%
% Completed 25/5/2017 JH
%% Starting script
orig_stim = 1 : n;
nstim = length(orig_stim);

lag_val = min : max;

if mod(n,length(lag_val))~=0
    msg=['Requires number of conditions to be multiples of '...
        'the number of lag values (length(min:max))'];
    error(msg)
end

orig_lag = repmat(lag_val, 1, nstim/length(lag_val));

ntrials = nstim * rep;

x = 1;
num_solution = 0;
order = zeros(num_seq,ntrials);

while num_solution<num_seq
    all_stim = orig_stim(randperm(nstim));
    all_lag = orig_lag(randperm(length(orig_lag)));
   
    seq_idx = num_solution + 1;
    order(seq_idx,:) = 0;
    
    j = 1;
    
    fprintf('Attempt #%d for sequence #%d \n', x, seq_idx);
    
    for i = 1 : nstim
        stim = all_stim(1);
        lag = all_lag(1);
        
        if i == 1
            order(seq_idx,j) = stim;
            order(seq_idx,j+lag) = stim;
            
            j = j + 1;
        else
            while order(seq_idx,j) ~= 0
                j = j + 1;
            end
            
            order(seq_idx,j) = stim;
            
            lag_found = 0;
            
            cycle = length(all_lag);
            
            while ~lag_found && cycle>0
                if (j+lag) > length(order(seq_idx,:)) || (order(seq_idx,j+lag) ~= 0)
                    all_lag = [all_lag, lag];
                    all_lag(1) = [];
                    lag = all_lag(1);
                    cycle = cycle - 1;
                else
                    lag_found = 1;
                end
            end
            
            if ~lag_found && cycle == 0
                break
            end
            
            order(seq_idx,j+lag) = stim;
        end
        
        all_stim(1) = [];
        all_lag(1) = [];
    end
    
    if ~any(order(seq_idx,1:ntrials)==0) && length(order(seq_idx,:)) == ntrials
        num_solution = num_solution + 1;
        fprintf('Generated sequence #%d \n', seq_idx);
    else
        x = x + 1;
    end
end

%% Check lags
check_lags=zeros(num_seq,nstim);

for i = 1 : num_seq
    for j = 1 : nstim
        idx = find(order(i,:)==j);
        check_lags(i,j) = idx(2) - idx(1);
    end
end

check_lags = sort(check_lags,2);

end
        
        
        
        
        
        