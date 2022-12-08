\\run with "gp -q -s "32G" path_to_this_file"
nmax = 100; \\specify nmax to compute the number of connected graphs for arbitrary number of vertices
\ps 101; \\ set it to nmax + 1
a(n) = {log( sum(m=0, n, (1 + y)^(m * (m-1)/2) * x^m/m!))}
b(x, n, k) = {polcoeff(polcoeff(x, n), k) * n!}
x = a(nmax);
{for(n=1, nmax, v = vector(binomial(n, 2) + 1); 
for(k = 0, binomial(n, 2), v[k + 1] = b(x, n, k)); 
write("num_connected_graphs.txt", v);)}
quit