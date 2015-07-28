clear
% preparing the problem
% trying to find a low approximation to A, an m x n matrix
% where m >= n
m = 2000;
n = 1000;
%// first let's produce example A
A = rand(m,n);
%
% beginning of the algorithm designed to find alow rank matrix of A
% let us define that rank to be equal to k
k = n/10;
% R is an m x l matrix drawn from a N(0,1)
% where l is such that l > c log(n)/ epsilon^2
%
l = 2*k;
% timing the random algorithm
trand =cputime;
R = randn(m,l);
B = 1/sqrt(l)* R' * A;
[a,s,b]=svd(B);
Ak = A*b(:,1:k)*b(:,1:k)';
trandend = cputime-trand
% now timing the normal SVD algorithm
tsvd = cputime;
% doing it the normal SVD way
[U,S,V] = svd(A,0);
Aksvd= U(1:m,1:k)*S(1:k,1:k)*V(1:n,1:k)';
tsvdend = cputime -tsvd