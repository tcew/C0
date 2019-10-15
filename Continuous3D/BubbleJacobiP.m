function phi = BubbleJacobiP(a, alpha, beta, n)
  
  phi = 0.25*(1-a).*(1+a).*SimplifiedJacobiP(a, alpha, beta, n);

