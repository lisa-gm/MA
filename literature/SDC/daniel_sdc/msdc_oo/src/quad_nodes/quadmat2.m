function [C,w,maxlen] = quadmat2(a,b,x,d,cache)
% [C,w] = quadmat2(a,b,x,d,nrGaussPts)
% Returns the collocation matrix C for cumulative integration at
% the points x on [a,b]. The barycentric blending parameter is d.
% If f is a vector of function values f(x), then g = C*f should
% be an approximation to g(x) with g = cumsum(f).
% Uses cumsum of chebfun.

persistent CACHE

if nargin < 5,
    cache = 0;
end

if cache == 0,
    CACHE = [];
%    disp('quadmat2: cache reset')
end

n = length(x) - 1;

if cache,
    for jc = 1:length(CACHE),
        % check for cache hit:
        if (CACHE{jc}.bminusa - (b-a))<10*eps && norm(CACHE{jc}.xminusa - (x-a))<1e-14 && CACHE{jc}.d == d,
            C = CACHE{jc}.C;
            w = CACHE{jc}.w;
            maxlen = CACHE{jc}.maxlen;
%            disp('quadmat2: cache hit')
            return
        end
    end
end

w = weights(n,d,x);
C = zeros(n+1,n+1);
%wb = waitbar(0,'quadmat2...');

maxlen = 0;

if 0,
    % stefan's version
    for k = 1:(n+1),
        %waitbar(k/(n+1),wb)
        fx = 0*x; fx(k) = 1; % k-th Lagrange function
        c = chebfun(@(xx) bcinterpol(w,x,fx,xx.'),[a,b]); % ! overriding adaptive choice for nr of pts
        maxlen = max(maxlen,length(c));
        c = cumsum(c); c = c(:);
        C(:,k) = c(x) - c(x(1));
    end
    %C(1,:) = 0;
end

if 1,
    % georges version
    if n == 0,
        C = 1;
        %close(wb)
        return;
    end
    % georges version
    for k = 1:ceil((n+1)/2),
        %waitbar(2*k/(n+1),wb)
        fx = 0*x; fx(k) = 1; % k-th Lagrange function
        c = chebfun(@(xx) bcinterpol(w,x,fx,xx.'),[x(1),x(end)]); % ! overriding adaptive choice for nr of pts
        maxlen = max(maxlen,length(c));
        c = cumsum(c); c = c(:);
        C(:,k) = c(x) - c.vals(1);
    end
    for k = ceil((n+1)/2)+1:n+1  % uses symmetry of the Lagrange function, see KB, BIT
        C(:,k) = C(end,n+1-k+1) - C(end:-1:1,n+1-k+1);
    end
end

%close(wb)

if cache,
    jc = length(CACHE) + 1;
    % check for cache hit:
    CACHE{jc}.bminusa = b-a;
    CACHE{jc}.xminusa = x-a;
    CACHE{jc}.d = d;
    CACHE{jc}.C = C;
    CACHE{jc}.w = w;
    CACHE{jc}.maxlen = maxlen;
%    disp('quadmat2: cache write')
end
end