function V = NSvertices()

% Vertices of the no-signaling polytope in the (2,2;2,2) Bell scenario:
% 16 local deterministic vertices + 8 PR-box vertices, 24 in total.
%
% Output: V(lambda, a, b, x, y), with a,b,x,y = 1,2 (i.e., the physical
% values a,b,x,y in {0,1} shifted by +1, matching the convention of
% detPax.m and is_LHV.m), and lambda = 1..24.
%
% Vertices 1..16 : deterministic, V = DA(mu,a,x)*DB(nu,b,y)
% Vertices 17..24: PR boxes, P_{rst}(ab|xy) = 1/2 if
%                  a (+) b = x*y (+) r*x (+) s*y (+) t  (mod 2),
%                  and 0 otherwise, for r,s,t in {0,1}.
%                  (r,s,t) = (0,0,0) is the canonical PR box a(+)b = xy.

DA = detPax(2,2);   % DA(mu, a+1, x+1)
DB = detPax(2,2);

V = zeros(24,2,2,2,2);

% ----- 16 local deterministic vertices -----
k = 0;
for mu = 1:4
    for nu = 1:4
        k = k + 1;
        for a = 0:1
        for b = 0:1
        for x = 0:1
        for y = 0:1
            V(k,a+1,b+1,x+1,y+1) = DA(mu,a+1,x+1)*DB(nu,b+1,y+1);
        end
        end
        end
        end
    end
end

% ----- 8 PR-box vertices -----
for r = 0:1
for s = 0:1
for t = 0:1
    k = k + 1;
    for a = 0:1
    for b = 0:1
    for x = 0:1
    for y = 0:1
        if mod(a+b,2) == mod(x*y + r*x + s*y + t, 2)
            V(k,a+1,b+1,x+1,y+1) = 1/2;
        end
    end
    end
    end
    end
end
end
end

% ----- sanity checks (uncomment on first use) -----
% for k = 1:24
%     for x = 1:2, for y = 1:2
%         assert(abs(sum(sum(V(k,:,:,x,y))) - 1) < 1e-12);          % normalization
%     end, end
%     for a = 1:2, for x = 1:2
%         assert(abs(sum(V(k,a,:,x,1)) - sum(V(k,a,:,x,2))) < 1e-12); % NS: A-marginal
%     end, end
%     for b = 1:2, for y = 1:2
%         assert(abs(sum(V(k,:,b,1,y)) - sum(V(k,:,b,2,y))) < 1e-12); % NS: B-marginal
%     end, end
% end

end
