
if Zadanie==1
% -------------------Zadanie 1------------------
x = -20:0.1:20 ; 
y = -5:0.1:20 ; 
figure;
[X,Y] = meshgrid(x,y);
ineq = (2*X-Y<=4) & (Y+X > 3) & (Y+4*X>=-2);
h = pcolor(X,Y,double(ineq)) ;
h.EdgeColor = 'none' ;
xlabel('x')
ylabel('y')

%Definiowanie problemu optymalizacji oraz zmiennych
prob=optimproblem('ObjectiveSense','max');
y=optimvar("y");
x=optimvar("x");
%Definiowanie ograniczeń + Wrzucenie ich do problemu 
ogr1 = 2*x-y<=4;
ogr2 = y+x >=3;
ogr3 = y+4*x >=-2;
prob.Constraints.cos1=ogr1;
prob.Constraints.cos2=ogr2;
prob.Constraints.cos3=ogr3;
prob.Objective=-y;
[sol,optval]=solve(prob);
disp(sol)
show(prob)
end
if Zadanie==2
% ZADANIE 2
X=linspace(0,4,1000);
y=X.^4-4*X.^3-2*X.^2+12*X+9;
figure;
plot(X,y)
[wartosc,probka]=min(y);
prob2=optimproblem('ObjectiveSense','min');
x=optimvar("x",'LowerBound',0,'UpperBound',inf);
y=optimvar("y");
f_celu=x.^4-4*x.^3-2*x.^2+12*x+9;
prob2.Objective=f_celu;
initialGuess.x=4;
[sol2,optval2]=solve(prob2,initialGuess,'Solver','fmincon');
disp(sol2)
show(prob2)
end