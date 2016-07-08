% This function is made by Ahmed ElTahan

%{
        This function is intended to estimate the parameters of a dynamic
        system of unknown time varying parameters using the 
        Recursive Least Squares with Exponential Forgetting Method (RLS).
        After an experiment, we get the inputs, the outputs of the system. 
        The experiment is operated with sample time Ts seconds. 
        The system here is transfer function in the form of:

                                            y         z^(-d) Bsys
                                Gp = ------ = ----------------------
                                            u               Asys
                                    
                                 Asys * y = z^(-d) Bsys * u + e

    where:
    -- y : output of the system.
    -- u : control action (input to the system).
    -- e : white guassian noise (noise with zero mean).
    -- Asys = 1 + a_1 z^-1 + a_2 z^-2 + ... + a_na z^(-na).
    -- Bsys = b_0 + b_1 z^-1 + b_2 z^-2 + ... + b_nb z^(-nb).
    -- d : delay in the system.

    Function inputs
    u : input to the system in column vector form
    y : input of the system in column vector form
    na : order of the denominator polynomail
    nb : order of the numerator polynomail
    d : number represents the delay between the input and the output
    Ts : sample time.
    lambda : forgetting factor  usually {0<lambda<1}
    
    Function Output
    Theta_final : final estimated parameters.
    Gz_estm : pulse (discrete) transfer function of the estimated parameters
    1 figure for the history of the parameters that are being estimated
    2 figure to validate the estimated parameters on the given output
    using the instantaneous estimated parameters.
    3 figure to plot the input versus time.

    An example is added to illustrate how to use the funcrtion
%}

function [ Theta_final Gz_estm] = RecursiveLeastSquaresWEF( u, y, na, nb, d, Ts, lambda )

nu = na + nb + 1; % number of unkowns
N = length(y); % length of the vectors
n = max(na+1, nb+d+1); % used to avoid zeros at the begining and added 1 to avoid zero index in matlab
ne = N - n + 1; % number of equations that will used in estimation

%% Initialization
% Covariance matrix Initialization
alpha = 10E6;
for i = 1:(n-1)
    P{1,i} = alpha*eye(nu, nu);
end

% Initialization of estimated parameters
theta = cell(1, N);
for i = 1:(n-1)
    theta{1, i} = zeros(nu,1);
end

% Initialization of cell for legend the estimated paramters
Name = cell(nu, 1);

% Initialization of y_estm
y_estm = zeros(N, 1);


%% Recursive Least Squares Algorithm
% Filling the Phi matrix
% % for j = 1 : N
 for j = n : N
    
        for i = 1 : na % this for loop used to fill parts in the same line that pertains to the output "denomenator"
            if ((j-i)<=0)
                Phi(j, i) = 0;
            else
                Phi(j, i) = -y((j-i));
            end
        end
    
        
        for i = na+1 : nu % this for loop used to fill parts in the same line that pertains to the input "numerator" starts from the last postion column of "na +1 
            if ((j-(i-na)-d)<=0) % add na as we left the output going to the input and we start index after "na"
                Phi(j, i) = 0;
            else
                Phi(j, i) = u((j-d-i + (na+1)));
            end
        end
        
        % Phi = [phi(i) phi(i+1) ... phi(N)], Phi is the matrix N*nu but
        % the phi is the vector nu*1
        phi = Phi(j, :)'; % to get the j th row and then convert into column and that what we need in RLS not Phi
        
        % RLS algorithm
%         if ((j-1)>0) % to skip zero and negative indecies
            P{1, j} = 1/lambda*(P{1, j-1} - P{1, j-1}*phi*inv(lambda+phi'*P{1, j-1}*phi)*phi'*P{1, j-1}); % Covariance Matrix
            K = P{1, j}*phi; % Gain
            theta{1, j} = theta{1, j-1} + K*(y(j)-phi'*theta{1, j-1}); % estimated parameters
            A = 1;
            for i = 1: na
                A(i+1) = theta{1, j}(i);
            end
            for i = na+1:nu
                B(i-na) = theta{1, j}(i);
            end
            y_estm(j) = outputestimation( A, B, d, u, y_estm, j ); % Instantaneous Estimated Parameter System Output
%         end
        
end
    
Theta = cell2mat(theta); % converting the Total_loop_estimate from cell array to matrix, the final column is the same as the final estimate
Theta_final = Theta(:, N); % final estimation of parameters

 %% Parameters and transfer function preparation (final Parameters)
 z = tf('z');
 A = z^(d+nb);
 B = 0;
for i = 1:na
    a(i) = Theta_final(i);
    A = A + a(i)*z^(d-na+nb+na-i);
end
for i = na+1:nu
    b(i-na) = Theta_final(i);
    B = B + b(i-na)*z^(nb - (i-(na+1)));
end

Gz_estm = B/A;
[num, den] = tfdata(Gz_estm);
num = cell2mat(num);
den = cell2mat(den);
Gz_estm = tf(num, den, Ts);

%% Validation
final_time = (N-1)*Ts;
t = 0:Ts:final_time;

%original system
h = figure;
hold all
plot(t, y, 'LineWidth', 2)

% estimated system
plot(t, y_estm,'MarkerSize',3, 'Marker','square','color', 'r', 'LineWidth', 2);
legend(['Ordinary System'], 'Instantaneous Estimated Parameter System Output');
xlabel('Time (sec.)');
ylabel('Output');
title('Output from both Real and Estimated Systems')
grid on 
print(h,'-dpng','-r500','Output')

% Parameter History Plot
for i = 1 : nu
    h = figure(2);
    hold all
    plot(t, Theta(i, :), 'LineWidth', 2);
    if (i>=1 & i<=na)
        Name{i, 1} = ['a', num2str(i)]; % Polynomal A parameters
    else
        Name{i, 1} = ['b', num2str(i-na-1)]; % Polynomal B parameters
    end
end

grid on
xlabel('Time (sec.)')
ylabel('Estimated Parameters')
title('Estimated Parameters History')
legend(Name)
print(h,'-dpng','-r500','Estimated Parameter History');

% Input plot over time
h = figure;
plot(t, u, 'LineWidth', 2);
xlabel('Time (sec.)');
ylabel('Input');
title('Input Versus Time')
grid on 
print(h,'-dpng','-r500','Input')
end
