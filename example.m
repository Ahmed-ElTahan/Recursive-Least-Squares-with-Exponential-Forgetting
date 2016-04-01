%
%{
        This is an example to show how to use the RecursiveLeastSquaresWEF (RLS)
        function. At first, Imagine that I have a transfer function that I
        know its parameters and its orders (na, nb and d). By exciting its output with a certain input
        which is here is a square wave, I can get the input and output of
        this dynamic system as similar as an experiment. To simulate the change of parameters
        during the experiment I choose to change the parameters of the transfer function
        manually and suddenly after some period of time and see how the estimation will change. 
        Now I use these input and output vectors with the RLS algorithm to estimate the
        parameters
%}
clc; clear all; close all;

% Gathering the input and the output of the original system

%% Input Initialization
Ts = 0.1; % sample time
t = 0:Ts:37;
u = square(t, 50);
u = u';

%% Transfer function preparation
s = tf('s');
Gs = 1/(s^2+s+1);
Gz = c2d(Gs, 0.2, 'zoh');
[num, den] = tfdata(Gz);
num = cell2mat(num);
den = cell2mat(den);

na = length(den) -1;
iodel = get(Gz, 'iodelay'); %input output delay
d = max(find(num ==0))+iodel; % total delay
nb = length(num) - max(find(num ==0)) - 1;
num = num(max(find(num ==0))+1:end); % to remove zeros at the first (comes from delay)
nu = na + nb + 1;

% Calculating A, B of the system
Asys = 1;
for i = 1: na
    Asys(i+1) = den(i+1);
end
        
for i = na+1:nu
    Bsys(i-na) = num(i-na);
end
         


%% Output Initialization
y = zeros(length(t), 1);
for i = 1 : length(t)
    
    % change paramter of a1 and b2
    if (i >= 100 & i < 250)
        Asys(2) = -1.6;
        Bsys(2) = 0.1;
    end
    % change paramter of a1 and b2 with different values
    if (i >= 250)
        Asys(2) = -1.5;
        Bsys(2) = 0.5;
    end
    y(i) = outputestimation( Asys, Bsys, d, u, y, i );
    
end
y = y';



%% Estimation using RecursiveLeastSquaresWEF function
[ theta, Gz_estm ] = RecursiveLeastSquaresWEF( u, y, 2, 1, 1, Ts, 0.95) % 2nd order estimation
