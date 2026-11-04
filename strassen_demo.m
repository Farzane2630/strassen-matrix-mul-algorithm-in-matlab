% Demonstration of Strassen's Matrix Multiplication

clc;
clear;

n = 10;
cutoff = 1;

rng(1);

A = rand(n);
B = rand(n);


disp(A);

disp(B);

tic;
C_strassen = strassen_multiply(A, B, cutoff);
t_strassen = toc;

fprintf('\n========== FINAL RESULT ==========\n');
fprintf('\nStrassen Result:\n');
disp(C_strassen);

tic;
C_builtin = A * B;
t_builtin = toc;

fprintf('Built-in Result:\n');
disp(C_builtin);

fprintf('\nStrassen Time: %.6f seconds\n', t_strassen);
fprintf('Built-in Time: %.6f seconds\n', t_builtin);

fprintf('Maximum Error: %.2e\n', ...
    max(abs(C_strassen(:) - C_builtin(:))));

