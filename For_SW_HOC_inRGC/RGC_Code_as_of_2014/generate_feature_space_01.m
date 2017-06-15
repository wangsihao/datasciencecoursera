% generate_feature_space_01.m
% JG: Uses 0 and 1 to compute probability space
%
% Makes a proba. distribution over N binary random vars.  There will be 2^N
% states.  "output:" 
% (1) prob_d  (col vect, length 2^N) -- proba of each state
% (2) feature_matrix (2^N rows, corresp each state)
% (3) feature_matrix_firstorder (2^N rows, corresp each state)

% feature_matrix.  Each row corresp. to a state.  Along each col, have
% value of a different feature in that state.  Number of cols is number
% of features.  These are the features that match in the corresp. max-ent
% model.  
% HERE -- we have  N + (N^2-N)/2 = (N^2+N)/2 features:  
%         first and second order moments for N variables) 

% feature_matrix_firstorder.  Fewer cols -- only fit features up to first
% order.  HERE -- take N first order moments.


% calls:  function max_ent_analysis_fn(prob_d,feature_matrix,feature_matrix_firstorder)
% ... which is based on IterativeScalingScript by J. Beggs
% Modified to maximize likelihood / minimize KL distance via conjugate
% gradient method, 11/12/08



%N=4 ;  % number of dimensions" (i.e. neurons)

%JG: most often N=3

%===============================
%%%% First ... list the **possible states** size(state) = #words x #neurons
% ---> THIS GOES INTO MATRIX "state"

state = zeros(2^N,N);  %matrix:  each row will be 1's and -1's -- in N slots
for j=0:(2^(N))-1  
    holder = dec2bin(j);  %dec2bin returns binary 1 and 0 value of number ...
    adjust = (N) - length(holder);   

    for ind=1:adjust
        holder = ['0' holder];   %pad with zeros if necessary up to length N
    end
    input=[];
    for ind=1:N
        input(ind) = str2num(holder((N+1)-ind));
    end
    %JG: can this flip be avoided if the padding is done [holder '0']?    
    state(j+1, :) = fliplr(input);
end


%%%%  compute a matrix of "features"
% Each row corresponds to one state (same order as above)
% Each col is the value of the feature in that state.
% We choose means and second moments as these features:
% First N cols are value of ith entry in state.
% Next (N^2-N)/2 are value of ith entry times value of jth entry, j > i.  Loop
% over i, then over increasing values of j.
% Total of N + (N^2-N)/2  = (N^2+N)/2 features per state (num columns)
    
feature_matrix=zeros( 2^N , (N^2+N)/2 ) ;
for state_index=1:2^N %go through all the rows
        %first fill in the first N cols
    for feature_index=1:N ;
        feature = state(state_index,feature_index);
        feature_matrix(state_index,feature_index) = feature ;
    end
        %Next do the next (N^2 - N)/2 cols
    for i=1:(N-1)
        for j=(i+1):N
            %JG: initial feature index is N
            feature_index=feature_index+1;
            feature = state(state_index,i)*state(state_index,j);
            feature_matrix(state_index,feature_index) = feature ;
        end
    end
end

    
%%%%  compute a matrix of "features"
% Each row in one state (same order as above)
% Each col is the value of the feature in that state.
% First N cols are value of ith entry in state.  That's it -- no other
% features.  Total of N features per state
    
feature_matrix_firstorder=zeros( 2^N , N ) ;
for state_index=1:2^N
        %first fill in the first N cols
    for feature_index=1:N ;
        feature = state(state_index,feature_index);
        feature_matrix_firstorder(state_index,feature_index) = feature ;
    end
end

%JG: only if N==3 there is a moment matrix
if N==3

moment_matrix=zeros( 2^N , (N^2+N)/2 + 1 ) ;
for state_index=1:2^N
        %first fill in the first N cols
    for moment_index=1:N ;
        moment = state(state_index,moment_index);
        moment_matrix(state_index,moment_index) = moment ;
    end
        %Next do the next (N^2 - N)/2 cols
        %JG: all possible pairwise products
    for i=1:(N-1)
        for j=(i+1):N
        moment_index=moment_index+1;
        moment = state(state_index,i)*state(state_index,j);
        moment_matrix(state_index,moment_index) = moment ;
        end
    end

    %Next do the last col
    %JG: one tripplet product
        moment_index=moment_index+1;
        moment = state(state_index,1)*state(state_index,2)*state(state_index,3);
        moment_matrix(state_index,moment_index) = moment ;
end

else 
    moment_matrix=[];
end


