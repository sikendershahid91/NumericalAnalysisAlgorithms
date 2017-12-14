clc; clear; %init

dt = 0.05;
T = 60; 
e_n = [-1, 0];
e_p = [1, 0];
p_s = [1, -0.5]; 

p = zeros(T+1, 2); 
p(1, :) = p_s; 

for t = 2:T+1
	v_p = (p(t-1, :) - e_p) / norm(p(t-1, :) - e_p)^2; 
	v_n = -(p(t-1, :) - e_n) / norm(p(t-1, :) - e_n)^2;
	v = (v_p + v_n); 
	p(t, :) = p(t-1, :) + dt * (v); 
end

scatter(e_n(1), e_n(2), 'O'); 
hold on; 
scatter(e_p(1), e_p(2), 'O'); 

scatter(p(:,1), p(:,2), [], 1:T+1, 'filled')

ylim([-2 2]); 
xlim([-3 3]); 
hold off;
