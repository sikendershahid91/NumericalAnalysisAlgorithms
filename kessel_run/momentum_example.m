dt = 0.01; 
T = 500; 
p_s = [0,0]; 
v =  [-3, 0]; 

G = 1; 

p_g = [-5, -5]; 
m = 100; 

p = zeros(T+1, 2); 
p(1, :) = p_s; 

for t = 2:T+1
	r = p_g - p(t-1,:); 
	r_mag = norm(r); % distance between s and e
	f_mag = G * m/(r_mag^2); 
	f = r/norm(r) * f_mag; 
	v = v+dt *f;
	p(t, :) = p(t-1, :) +dt *v; 

end

scatter(p(:,1), p(:,2), 'O'); 
hold on; 
scatter(-5, -5, 'R');
hold off