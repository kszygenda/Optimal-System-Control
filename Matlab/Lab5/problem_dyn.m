function rozwiazanie = problem_dyn(a)
rozwiazanie=ode45(@model,[0 1],'none');
end

