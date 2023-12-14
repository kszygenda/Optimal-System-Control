clc; clear all;
% Nauka
prob = optimproblem("Description","Factory Location");
x = optimvar("x");
y = optimvar("y");
X = [5 40 70];
Y = [20 50 15];
d=sqrt((x-X).^2+(y-Y).^2);
dTotal=sum(d);
prob.Objective=dTotal;
xvec = linspace(0,75);
yvec = linspace(0,75);
[x,y] = meshgrid(xvec,yvec);
distance = sqrt((x-X(1)).^2 + (y-Y(1)).^2)+...
    sqrt((x-X(2)).^2 + (y-Y(2)).^2)+...
    sqrt((x-X(3)).^2 + (y-Y(3)).^2);
contourf(x,y,distance)
ylabel("Y-Coordinate")
xlabel("X-Coordinate")
colorbar

% jedzenie
load Nutrition.mat 
food

prob = optimproblem("Description","An Optimal Breakfast");
servings = optimvar("servings",16,"LowerBound",0);
C = food.Price .* servings;
prob.Objective = sum(C);
cals=servings .* food.Calories;
totalCals=sum(cals);
prob.Constraints.calories=totalCals == 350;
[sol,optval]=solve(prob);
optServings=sol.servings;
bar(food.Name,optServings)
check=evaluate(totalCals,sol);


