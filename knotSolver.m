function [P,Q] = knotSolver(order,alpha)
% Written by Magnus V. Paludan, 2023. 
% The function takes as input a knot order (e.g., [1 3 5 6 4 2]) and
% the value for alpha (which is the ratio between the loop and overlap
% resistance.

% We first calculate the positions where the knot loops under (even
% numbers in the ordering) and positions where the knot loops over (odd
% numbers in the ordering).

% The even numbers are for calculating the additional compressional
% resistance, while the odd number expansion is neglected.
eId = find(mod(order,2)==0);
oId = find(mod(order,2)==1);

% Initialise "x" (which is the relative flow rate) and a list of
% equations to be updated iteratively.
syms x
syms eqs [1 length(order)]
syms eqs2 [1 length(order)]

% All the odd label equations are 0 since the added resistance is
% ignored.
eqs(oId) = 0;
eqs2(oId) = 0;

% We loop through the even numbers in the knot order
for i=1:length(eId)
    % Find distance between pairs (e.g., 3 and 4 in the knot order).
    L = eId(i)-find(order==order(eId(i))-1);
    % If the distance is 1, that means that the overlap only
    % interacts with itself.
    if L == 1
        eqs(eId(i)) = @(x) funf(L*x);
        eqs2(eId(i)) = @(x) 1./funf(L*x);
        % If the distance is more than 1, we have to add the resistance of other overlaps.
    else
        eqs(eId(i)) = @(x) funf(L*x + 2*alpha*x*sum(eqs(eId(i)-L:eId(i)-1)));
        eqs2(eId(i)) = @(x) 1./(funf(L*x + 2*alpha*x*sum(eqs(eId(i)-L:eId(i)-1))));
    end
    % The array of "eqs2" is to evaluate when the conductivity
    % (1/R) goes to 0, instead of evaluating when the resistance
    % (R) goes to infinity. The array "eqs" contains the
    % resistances.
end
eqs;
Pfun = matlabFunction(2*alpha*(length(eId)/alpha + sum(eqs))*x);

Qs = 0:0.0001:1; 

P = [];
Q = []; 

Pe = 0; 
count = 1; 
while Pe<25
    Q(end+1) = Qs(count); 
    P(end+1) = Pfun(Qs(count));
    count =  count + 1; 
    Pe = P(end); 
end
end
function f = funf(x)
    f = (1-x)^(-3); 
end