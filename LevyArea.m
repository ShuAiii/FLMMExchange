function Area = LevyArea(tau,n)
subtau=tau/n;

mean = zeros(2, 1);
cov = eye(2);

subbm=mvnrnd(mean, cov, n)';

B = sqrt(subtau) * tril(ones(n));

C = -diag(ones(n - 1, 1), -1);

D = diag(ones(n - 1, 1), 1);

M = C + D;

Area = subbm(1,:) * B' * M * B * subbm(2,:)';



