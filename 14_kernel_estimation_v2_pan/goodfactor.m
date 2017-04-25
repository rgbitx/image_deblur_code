function N=goodfactor(N)
    
    f=factor(N);
    
    while(max(f)>7)
      N=N+1;
      f=factor(N);
    end
    