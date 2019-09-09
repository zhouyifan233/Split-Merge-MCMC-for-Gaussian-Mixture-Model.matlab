mu1 = [-5 5];          % Mean of the 1st component
sigma1 = [2 0.5; 0.5 2]; % Covariance of the 1st component
mu2 = [5 0];        % Mean of the 2nd component
sigma2 = [2 -0.5; -0.5 2];  % Covariance of the 2nd component
mu3 = [-5 0];        % Mean of the 2nd component
sigma3 = [2 0; 0 2];  % Covariance of the 2nd component
mu4 = [5 5];        % Mean of the 2nd component
sigma4 = [2 -0.5; -0.5 2];  % Covariance of the 2nd component

r1 = mvnrnd(mu1,sigma1,180);
r2 = mvnrnd(mu2,sigma2,200);
r3 = mvnrnd(mu3,sigma3,220);
r4 = mvnrnd(mu4,sigma4,200);
label = [ones(1, 180), 2*ones(1, 200), 3*ones(1, 220), 4*ones(1, 200)];
X = [r1; r2; r3; r4]';
