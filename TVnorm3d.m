function y = TVnorm3d(x)

N = size(x);

id = [2:N(1),N(1)];

ir = [2:N(2),N(2)];

ib = [2:N(3),N(3)];


y = sum(sum(sum(sqrt((x(:,ir,:) - x).^2 + ...
    (x(id,:,:) - x).^2 + ...
    (x(:,:,ib) - x).^2))));

end