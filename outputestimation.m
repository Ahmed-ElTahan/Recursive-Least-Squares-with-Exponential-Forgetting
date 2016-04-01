% This function is made by Ahmed ElTahan
%{


    Any system can be written as 

                                            z^(-d) B                y
                                G = -------------------- = -----------
                                                   A                    u

    This function is intended to find the output of a system like a
    transfer function given the input, the system which you can have the
    orders (na, nb, d), the previous output and the point at which you want
    to estimate the system at. This is like a difference equation
    estimation for a system.

    Notice that in difference equation you have "u" inputs as a vector and
    "y" may not be estimated yet and hence you can initialize "y" with
    zeros and each point estimated should be added to "y" as the present
    output depends on the previous outputs and the inputs, then, the
    function should be like that 

        y = zeros(1, length(u))
        for m=1:length(u)
            y(m) = outputestimation( A, B, d, u, y, m );
        end

    This function is going to be benifical if the A, B matrices are going
    to be changed over time such as in Recursive Least Squares algorithm
    such that each loop we are going to find "theta" which contain A, B.
    Hence, we can calculated the estimated ouput based on the estimated
    parameters using this function as 

        y_estm = zeros(1, length(u))
        for m=1:length(u)
            y_estm(m) = outputestimation( A, B, d, u, y_estm, m );
        end

%}


function [ y_output ] = outputestimation( A, B, d, u, y, m )
% orders
na = length(A) - 1;
nb = length(B) -1;
d;
nu = na + nb + 1;

% coefficients
A = -A(2:end); % to delete "1" at the first and -ve for the difference eqn
B;

for (i=1:na)
        if ((m-i)<=0) % zero or negative indices shall give the ouput zero value
            Y(i) = 0;
        else
            Y(i) = y(m-i);
        end
end
for (i = na+1 : nu)
        if ((m-(i-na)-d)<=0)    % zero or negative indices shall give the input zero value
            U(i-na) = 0;
        else
            U(i-na) = u((m-d-i + (na+1)));;
        end
    end

y_output = A*Y'+B*U'; % resulting "y" at point "m"

end

