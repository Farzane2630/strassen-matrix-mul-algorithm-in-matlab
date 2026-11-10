% Demonstration of Strassen's Matrix Multiplication

clc;
clear;

n = 7;
cutoff = 2;

rng(1); % Generates the same random numbers every time we execute the program

A = rand(n);
B = rand(n);


% disp(A);
% disp(B);

tic; %Starts a stopwatch.
C_strassen = strassen_multiply(A, B, cutoff);
t_strassen = toc; %Stops the stopwatch and stores the elapsed execution time

fprintf('\n========== FINAL RESULT ==========\n');
fprintf('\nStrassen Result:\n');
disp(C_strassen);



% Running MATLAB's built-in multiplication
%tic;
%C_builtin = A * B;
%t_builtin = toc;
%fprintf('Built-in Result:\n');
%disp(C_builtin);
%fprintf('\nStrassen Time: %.6f seconds\n', t_strassen);
%fprintf('Built-in Time: %.6f seconds\n', t_builtin);


