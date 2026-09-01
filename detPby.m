function y = detPby(ny,nb)

% table of deterministic probability distributions for single side

for lambda = 1:nb^ny
    
    v(lambda,:) = dec2base(lambda-1,nb,ny)-'0';
    
end

for lambda = 1:nb^ny
for b=0:nb-1
for y=0:ny-1

  y(lambda,b+1,y+1) = kronDel(b,v(lambda,y+1));

end
end
end

end