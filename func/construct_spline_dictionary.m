

function B = construct_spline_dictionary(signal_length,num_space,degree)

% degree    = 4; % cubic, set as default
% num_space = 10; % for equally spaced knot vector

if nargin < 2
    degree    = 4; % cubic, set as default
    num_space = 10; % for equally spaced knot vector
end

x  = signal_length;
s1 = min(x);
s2 = max(x);
step        = (s2-s1)/num_space;
knot_vector = s1:step:s2;

knot_vector =  [s1, s1, s1,knot_vector, s2, s2, s2];

[B,~] = bspline_basismatrix(degree,knot_vector,x);


%% add Intercept 

Intercept_p = ones(length(signal_length),1);
Intercept_n = -1*ones(length(signal_length),1);

B = [Intercept_p,B];

% for i = 1:w
%     x = X(:,i);
%
%     s1 = min(x); s2 = max(x);
%     step = (s2-s1)/num_space;
%
%     knot_vector = s1:step:s2;
%
%     [B,~] = bspline_basismatrix(degree,knot_vector,x);
%
%     spline_d = [spline_d,B];
% end