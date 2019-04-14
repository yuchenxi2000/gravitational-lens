% 用来显示光子轨迹的demo
% cc 建议取 1 / 27.1, 1 / 27.01 这类值
% range 是显示范围，自己调。 10 差不多了。

function hole(cc, range)
if cc > 1 / 27 || cc <= 0 || range <= 0
    return;
end
theta = (0 : 0.001 : 1) * 10 * pi;
orbit = PhotonOrbit(cc, theta);

r = orbit{1};
format long;
disp(['theta_m = ', num2str(orbit{2})])
disp(['d = ', num2str(orbit{3})])

index = 0;

for i = 1:1001
    if r(i) < 0 || r(i) > range
        index = i;
        break;
    end
end

theta = theta(1:index);
r = r(1:index);

polar(theta, r);
hold on;
polar(-theta, r);
% hold on;
% rr = -theta./ theta .* orbit{3} .* (sin(theta - orbit{2})).^(-1);
% 
% polar(theta, rr);

end