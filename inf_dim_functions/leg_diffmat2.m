function D = leg_diffmat2(n)

D = sparse(n,n);

for jj=2:n
    nn = (jj-2):-2:0;
    D((jj-1:-2:1),jj) = 2*nn+1;
end

end