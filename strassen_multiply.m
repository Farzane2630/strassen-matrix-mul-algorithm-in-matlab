
function C = strassen_multiply(A, B, cutoff)

    n = size(A, 1);
    m = 2^ceil(log2(n));

    Ap = zeros(m);
    Bp = zeros(m);

    Ap(1:n, 1:n) = A;
    Bp(1:n, 1:n) = B;


    C = strassen_rec(Ap, Bp, cutoff);

    C = C(1:n, 1:n);

end


function C = strassen_rec(A, B, cutoff)

    n = size(A, 1);


    if n <= cutoff


        disp(A);
        disp(B);

        C = A * B;

        disp(C);

        return;

    end

    h = n / 2;

    % Divide matrix A into four parts
    A11 = A(1:h, 1:h);
    A12 = A(1:h, h+1:n);
    A21 = A(h+1:n, 1:h);
    A22 = A(h+1:n, h+1:n);

    % Divide matrix B into four parts
    B11 = B(1:h, 1:h);
    B12 = B(1:h, h+1:n);
    B21 = B(h+1:n, 1:h);
    B22 = B(h+1:n, h+1:n);




    % Calculate seven Strassen products


    M1 = strassen_rec(A11 + A22, B11 + B22, cutoff);

    M2 = strassen_rec(A21 + A22, B11, cutoff);

    M3 = strassen_rec(A11, B12 - B22, cutoff);

    M4 = strassen_rec(A22, B21 - B11, cutoff);

    M5 = strassen_rec(A11 + A12, B22, cutoff);

    M6 = strassen_rec(A21 - A11, B11 + B12, cutoff);

    M7 = strassen_rec(A12 - A22, B21 + B22, cutoff);

    % Combine the seven products


    C11 = M1 + M4 - M5 + M7;
    C12 = M3 + M5;
    C21 = M2 + M4;
    C22 = M1 - M2 + M3 + M6;

    disp(C11);

    disp(C12);

    disp(C21);

    disp(C22);

    C = [C11, C12; C21, C22];

    disp(C);

end

