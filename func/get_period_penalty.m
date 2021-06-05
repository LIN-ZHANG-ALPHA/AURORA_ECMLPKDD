
function [periods,penalty_vector]  =  get_period_penalty(Pmax,Penalty_type)

% Pmax is the max period in NMP
% get penalty for periods in NPM

if nargin < 2
    Penalty_type = 'square';
end

periods = [];
for i=1:Pmax
    k=1:i;k_red=k(gcd(k,i)==1);k_red=length(k_red);
    periods=cat(1,periods,i*ones(k_red,1));
end

switch Penalty_type
    case 'square'
        penalty_vector = periods.^2;
        
    case 'linear'   
        penalty_vector = 2 * periods; 
        
    otherwise
        fprintf('Error, no such data is found! Try again!\n')
end


