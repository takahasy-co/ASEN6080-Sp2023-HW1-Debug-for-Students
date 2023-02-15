function [rho,rhodot] = Range_RangeRate(x_target_N,x_observer_N)
r_target_N   = x_target_N(1:3);
v_target_N   = x_target_N(4:6);
r_observer_N = x_observer_N(1:3);
v_observer_N = x_observer_N(4:6);

rho    = norm(r_target_N-r_observer_N);
rhodot = dot((r_target_N-r_observer_N),(v_target_N-v_observer_N))/rho;
end