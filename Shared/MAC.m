function mac = MAC(u,v)
[Mu,Nu] = size(u);
[Mv,Nv] = size(v);

mac = [];
if (Nu ~= Nv) || (Mu ~= Mv)
    return;
end;

N = Nu;
M = Mv;

mac = zeros(M,M);
for n = 1:M
    for m = 1:M
        mac(n,m) = (u(n,:) * v(m,:)')^2 /((u(n,:) * u(n,:)') * (v(m,:) * v(m,:)'));
    end
end

% mac = log(mac);
end