\\run with "gp -q -s "32G" path_to_this_file"
nmax = 100; \\specify nmax to compute the number of connected graphs for arbitrary number of vertices
B(nn, e=2) = {my(v=vector(nn));
for(n=1, nn, v[n] = e^(n*(n-1)) - sum(k=1, n-1, binomial(n, k) * e^((n-1) * (n-k)) * v[k])); v}
Strong(n, e=2) = {my(u=B(n, e), v=vector(n)); v[1]=1; 
for(n=2, #v, v[n] = u[n] + sum(j=1, n-1, binomial(n-1, j-1) * u[n-j] * v[j])); v}
row(n) = { Vecrev(Strong(n, 1+'y)[n]) }
{ for(n=1, nmax, write("num_connected_digraphs.txt", row(n))) }
quit