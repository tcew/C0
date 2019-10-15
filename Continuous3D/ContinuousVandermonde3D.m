function V = ContinuousVandermonde3D(N, r, s, t)

  [a,b,c] = rsttoabc(r,s,t);

  cnt = 1;

  Np = (N+1)*(N+2)*(N+3)/6;
  
  V = zeros(length(a),Np);

  %% offset for Jacobi indices
  si = 1;
  
  %% vertices
  V(:,cnt) = (1-a).*(1-b).*(1-c)/8; cnt = cnt+1;
  V(:,cnt) = (1+a).*(1-b).*(1-c)/8; cnt = cnt+1;
  V(:,cnt) = (1+b).*(1-c)/4; cnt = cnt+1; %% ??
  V(:,cnt) = (1+c)/2; cnt = cnt+1;

  %% AB
  for i=0:N-2
    Ba = BubbleJacobiP(a, si, si, i);
    V(:,cnt) = Ba.*((0.5*(1-b)).^(i+2)).*((0.5*(1-c)).^(i+2)); cnt=cnt+1;
  end

  %% AC
  for i=0:N-2
    Bb = BubbleJacobiP(b, si, si, i);
    V(:,cnt) = 0.5*(1-a).*Bb.*((0.5*(1-c)).^(i+2)); cnt=cnt+1;
  end

  %% BC
  for i=0:N-2
    Bb = BubbleJacobiP(b, si, si, i);
    V(:,cnt) = 0.5*(1+a).*Bb.*((0.5*(1-c)).^(i+2)); cnt=cnt+1;
  end
    
  %% AD
  for i=0:N-2
    Bc = BubbleJacobiP(c, si, si, i);
    V(:,cnt) = 0.25*(1-a).*(1-b).*Bc; cnt=cnt+1;
  end

  %% BD
  for i=0:N-2
    Bc = BubbleJacobiP(c, si, si, i);
    V(:,cnt) = 0.25*(1+a).*(1-b).*Bc; cnt=cnt+1;
  end


  %% CD
  for i=0:N-2
    Bc = BubbleJacobiP(c, si, si, i);
    V(:,cnt) = 0.25*(1+b).*Bc; cnt=cnt+1;
  end

  %% ABC
  for i=0:N-3
    for j=0:N-3-i
      Ba = BubbleJacobiP(a, si, si, i);
      Bb = BubbleJacobiP(b, 2*i+2+si, si, j);
      V(:,cnt) = Ba.*((0.5*(1-b)).^(i+1)).*Bb.*((0.5*(1-c)).^(i+j+3)); cnt=cnt+1;
    end
  end

  %% ABD
  for i=0:N-3
    for j=0:N-3-i
      Ba = BubbleJacobiP(a, si, si, i);
      Bc = BubbleJacobiP(c, 2*i+2+si, si, j);
      V(:,cnt) = Ba.*((0.5*(1-b)).^(i+2)).*Bc.*((0.5*(1-c)).^(i+1)); cnt=cnt+1; %% +2 ?
    end
  end
  

  %% ACD
  for i=0:N-3
    for j=0:N-3-i

      Bb = BubbleJacobiP(b, si, si, i);
      Bc = BubbleJacobiP(c, 2*i+2+si, si, j);
      
      V(:,cnt) = 0.5*(1-a).*Bb.*((0.5*(1-c)).^(i+1)).*Bc; cnt=cnt+1; %%  +2 ?
    end
  end


  %% BCD
  for i=0:N-3
    for j=0:N-3-i

      Bb = BubbleJacobiP(b, si, si, i);
      Bc = BubbleJacobiP(c, 2*i+2+si, si, j);

      
      V(:,cnt) = 0.5*(1+a).*Bb.*((0.5*(1-c)).^(i+1)).*Bc; cnt=cnt+1; %%  +2 ?
    end
  end

  %% interior

  for i=0:N-4
    for j=0:N-4-i
      for k=0:N-4-i-j

	Ba = BubbleJacobiP(a, si, si, i);
	Bb = BubbleJacobiP(b, 2*i+3, si, j);
	Bc = BubbleJacobiP(c, 2*i+2*j+5, si, k);
	
	V(:,cnt) = Ba.*((0.5*(1-b)).^(i+1)).*Bb.*((0.5*(1-c)).^(i+j+2)).*Bc; cnt = cnt+1;
      end
    end
  end
      
