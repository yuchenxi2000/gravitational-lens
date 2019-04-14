% solve the photon equation of general relativity
% author : yuchenxi
% @in :  cc        dt = du / sqrt(2 * u^3 - u^2 + cc)
%        theta     array of theta
% @out : r         array of r
%        theta_m   绕行角的一半
%        d         渐近线到极点距离

function orbit = PhotonOrbit(cc, theta)
p = [2, -1, 0, cc];
e = sort(roots(p));
e1 = e(3);
e2 = e(2);
e3 = e(1);

w = sqrt((e1 - e2) / 2.);
m = (e2 - e3) / (e1 - e2);

s = jacobiSN(w * theta, -m);
r = (e2 - (e2 - e3) * (s).^2).^(-1);

m2 = m / (1 + m);
theta0 = acos(sqrt(e2 / (e2 - e3)));

theta_m = 1 / sqrt(1 + m) / w * (ellipticF(pi / 2, m2) - ellipticF(theta0, m2));

% k = e2 / (e2 - e3);

% d = sqrt(2 / (1 - k ^ 2) / (e1 - e2 + (e2 - e3) * k ^ 2)) * 0.5 / sqrt(e2 * (e2 - e3));

d = (e2 - e3) / sqrt(-2.0 * e2 * e3 * (2 * e2 - e3) * ((e1 - e2) * (e2 - e3) + e2 * e2));
% d = 1 / sqrt(-2 * e2 * e3 * (e1 - e2) * (1 + e2 * (e1 - e2) / (e2 - e3)^2));

orbit = {r, theta_m, d};

end

