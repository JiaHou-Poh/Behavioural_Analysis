function[order, cum_rep, check_seq,num_rep] = GenSequence3(n,rep,maxr,num_seq)
%% ----------------- Script description -----------------------%%
% Generate sequence for item presentation based on set criterion for
% how many consecutive repetition is allowed. Takes in the following inputs:
% 1) n: Number of items
% 2) rep : Number of repetition per item. Takes in a single integer or an
%           array of length n.
% 3) maxr: Maximum number of repetitions allowed
% 4) num_seq: Number of sequence required
%
% If the number of repetitions across condition is important, choose 
% sequences based on 'cum_rep'. Cum_rep provides the number of repetition
% in each of the generate sequence for each of the conditions.
% 
% Largest possible number of repetition for any unconstrained sequence is 
% (n*rep)-n. If maxr > rep, sequence will be unconstrained.
%
% Completed 16/ 10/ 2017 JH
%% Starting script
orig_stim = 1 : n;
nstim = length(orig_stim);
num_rep = zeros(1,nstim);

if length(rep) == 1
    ntrials = nstim * rep;
    full_array = repmat(orig_stim,1,rep);
    num_rep(1,orig_stim) = rep;
else
    if length(rep)~=nstim
        error('Repetition number length mismatch')
    else
        full_array = [];
        ntrials = sum(rep);
        num_rep = rep;
        
        for i = 1 : n
            temp_seq = repmat(i,1,rep(i));
            full_array=[full_array,temp_seq];
        end
    end
end
        

x = 1;
num_solution = 0;

order = zeros(num_seq,ntrials);

while num_solution < num_seq
    all_cond = full_array(randperm(ntrials));
    
    seq_idx = num_solution + 1;
    order(seq_idx,:) = 0;
    
    j = 1;
    
    if mod(x,10) == 0
        fprintf('Attempt #%d for sequence #%d \n', x, seq_idx);
    end
    
    for i = 1 : ntrials
        cond = all_cond(1);
        chks = zeros(maxr,1);
        
        if i < maxr+1
            order(seq_idx,j) = cond;
            
            j = j + 1;
            
        else
            while order(seq_idx,j) ~= 0
                j = j + 1;
            end
            
            order(seq_idx,j) = cond;
            
            rep_found = 0;
            
            cycle = length(all_cond);
            
            %Check for repeats
            for k = 1 : maxr
                chks(k,1) = order(seq_idx,j-k);
            end
            
            if length(unique(chks)) == 1; % I.e. all are the same
                while ~rep_found && cycle > 0
                    if order(seq_idx,j) == chks(1)
                        all_cond = [all_cond, cond];
                        all_cond(1) = [];
                        cond = all_cond(1);
                        cycle = cycle - 1;
                        
                        order(seq_idx,j) = cond;
                    else
                        rep_found = 1;
                    end
                end
                
                if ~rep_found && cycle==0
                    order(seq_idx,j) = 0;
                    break
                end
            end
        end
        all_cond(1) = [];
    end
    
    % Check if the number of trials match for each condition
    numtrial = ones(1,nstim);
    for i = 1 : n
       numtrial(1, i) = length(find(order(seq_idx,:)==i));
    end
     
    chk_num = numtrial - num_rep;
    chk_num = sum(chk_num);
    
    if ~any(order(seq_idx,1:ntrials) == 0) && length(order(seq_idx,:)) == ntrials && chk_num == 0;
        num_solution = num_solution + 1;
        x = 1;
        fprintf('Generated sequence #%d \n', seq_idx);
    else
        x = x + 1;
    end
end

%% Checking
% Check repetition for each condition
check_seq=zeros(num_seq,n);
for i = 1 : num_seq
    for j = 1 : n
        check_seq(i,j) = length(find(order(i,:)==j));
    end
end

% Check the cumulative repetition for each condition
cum_rep = zeros(num_seq,n);
for i = 1 : num_seq
    for j = 2 : ntrials
        id = order(i,j);
        if id == order(i,j-1)
           cum_rep(i,id) = cum_rep(i,id) + 1;
        end
    end
end

num_rep = zeros(num_seq,n);
for i = 1 : num_seq
    for j = 1 : n
        num_rep(i,j) = length(find(order(i,:)==j));
    end
end

    
    

