
function  D = difference_operator_matrix(A,order)

% A is spline dictionary

[m,n] = size(A);

if nargin < 2
    order = 'second';
end


switch order
    case 'first'
        D        =  zeros(n-1,n);
        D(1,1:2) =  [-1,1];
        
        for i = 2:n-1
            D(i,:) = circshift(D(1,:),i-1);
        end
        
    case 'second'
        
        D        =  zeros(n-2,n);
        D(1,1:3) =  [1,-2,1];
        
        for i = 2:n-2
            D(i,:) = circshift(D(1,:),i-1);
        end
        
    case 'third'
        
        D        =  zeros(n-3,n);
        D(1,1:4) =  [-1,3,-3,1]; 
        
        for i = 2:n-3
            D(i,:) = circshift(D(1,:),i-1);
        end
        
    otherwise
        warning('Unexpected order difference. No matrix created.')
end

%% intercept in spline as zero in D
% D_intercept =  zeros(size(D,1),1);
% D = [D_intercept,D];

