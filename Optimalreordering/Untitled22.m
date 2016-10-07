n = 1024;
m = 2*n + 2;
B_construct = zeros(m);
B_construct(1, 2:n+1) = 1;
B_construct( 2*n+2 - n:2*n+2 -1, 2*n+2 ) = 1;