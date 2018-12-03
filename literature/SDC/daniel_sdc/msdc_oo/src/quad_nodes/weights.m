
% --- Compute the barycentric FH weights
function w=weights(n,d,x)
for k=0:n
	ji=max(k-d,0);
	jf=min(k,n-d);
	m=1;
	sumcoeff=[];
	product=1;
	for i=ji:jf
		l=1;
		prodterm=[];
		for j=i:i+d
			if(j==k)
				prodterm(l)=1;
			else
				prodterm(l)=(x(k+1)-x(j+1));
			end
			l=l+1;
		end
		product=1/prod(prodterm);
		sumcoeff(m)=((-1)^(i-1))*product;
		m=m+1;
	end
	[Y,I]=sort(abs(sumcoeff));
	Y(:)=sumcoeff(I(:));	
	w(k+1)=sum(Y);
end
w=w(:);
