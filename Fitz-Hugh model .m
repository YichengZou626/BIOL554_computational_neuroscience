
%  CLEAN UP THE WORKSPACE
clear
close all
clc

%  INITIALIZE FITZ-HUGH-NAGUMO EQUATION PARAMETERS, AND DEFINE A RANGE OF V
%  VALUES TO USE
current = 0;  %  THIS IS THE "I" VALUE IN THE PROBLEM STATEMENT ***
a = 0.3;
epsilon = 0.1;
gam = 0.5;
v_values = -0.6:0.1:1.6;  %  Defines an array of membrane potentials


%  DETERMINE NULLCLINES AND EQUILIBRIA
%  First, determine the w values that give v nullclines
w_for_v = -v_values.^3+(a+1)*v_values.^2-a*v_values+current;  %  REPLACE THE 0 WITH THE EQUATION FOR V NULLCLINES ***

%  Now, determine the w values that give w nullclines
w_for_w = v_values/gam;  %  REPLACE THE 0 WITH THE EQUATION FOR w NULLCLINES ***

%  Determine the equilibrium points by using the "roots" function, since
%  the Fitz-Hugh-Nagumo equations involve polynomials.
v_polynomial = [-1 a+1 -1/gam-a current];   %  REPLACE THE 0 WITH A 1XN ARRAY THAT CONTAINS THE COEFFICIENTS OF THE NECESSARY POLYNOMIAL ***
v_roots = roots(v_polynomial);

%  Get rid of roots with imaginary parts because they are nonphysical
v_roots = remove_imag(v_roots);
w_roots = v_roots/gam;


%  DETERMINE WHETHER THIS PARTICULAR SYSTEM IS STABLE/UNSTABLE, AND WHETHER
%  IT HAS OSCILLATIONS
for counter = 1:length(v_roots)
    %  If we aren't dealing with a real root, continue, otherwise, do
    %  something
    if isnan(v_roots(counter))
        continue
    else
        %  Compute the Jacobian using the included "jac" function
        %  YOU NEED TO MODIFY THIS FUNCTION ***
        J = jac_starter(v_roots(counter), w_roots(counter), a, gam, epsilon);
        
        %  Compute the eigenvalues
        %  YOU NEED TO FILL IN HOW TO COMPUTE EIGENVALUES HERE.  YOUR
        %  EIGENVALUES SHOULD BE PLACED IN THE MATRIX BELOW ***
        tr = J(1,1)+J(2,2);
        det = J(1,1)*J(2,2)-J(2,1)*J(1,2);
        s = [(tr+sqrt(tr.^2-4*det))/2, (tr-sqrt(tr.^2-4*det))/2];
        
        %  Isolate the real and imaginary parts of the eigenvalues
        s1_real = real(s(1));
        s1_imag = imag(s(1));
        s2_real = real(s(2));
        s2_imag = imag(s(2));
        
        %  Display some results.  This will only display anything if the
        %  roots are not nan
        disp('The Jacobian is ');
        disp(J);
        disp(['The eigenvalues are: ', num2str(s)])
        
        
        %  Determine if the equilibria are stable or unstable
        %  NOTE: "||" Means "or", and "&&" means "and"
        if  s1_real>0 || s2_real>0  %  YOU NEED TO ADD IN SOME CONDITIONS HERE ***
            str1 = ['Unstable Equilibrium at (V = ', num2str(v_roots(counter)), ', w = ', num2str(w_roots(counter)), ') --> '];
        elseif  s1_real<0 &&  s2_real<0  %  YOU NEED TO ADD IN SOME CONDITIONS HERE ***
            str1 = ['Stable Equilibrium at (V = ', num2str(v_roots(counter)), ', w = ', num2str(w_roots(counter)), ') --> '];
        end
        
        %  Determine if the equilibria contain oscillations
        if  s1_imag==0 && s2_imag==0   %  YOU NEED TO ADD IN A CONDITION HERE ***
            str2 = 'No Oscillations';
        else
            str2 = 'Oscillations';
        end
        
        disp([str1, str2]);
    end
end

%  ----------------DETERMINE VECTORS TO MAKE A VECTOR FIELD----------------
%  ---------------------YOU DON'T HAVE TO TOUCH THIS-----------------------
%  Define an array of w values (we already defined an array of v values)
w_values = -0.4:0.1:1.4;

%  Use "meshgrid" to make grided versions of v and w so that we have
%  something compatible with quiver
[V_grid, W_grid] = meshgrid(v_values, w_values);

%  Compute grided versions dv/dt and dw/dt so that we have something
%  compatible with the "quiver" function


%  YOU NEED TO ENTER THE FITZ HUGH NAGUMO EQUATIONS HERE, BUT USING THE
%  "V_grid" AND "W_grid" VARIABLES.  I'VE JUST PUT SOMETHING IN AS A 
%  PLACEHOLDER ***
dv = -V_grid.^3+(a+1)*V_grid.^2-a*V_grid-W_grid+current;
dw = epsilon*(V_grid-gam*W_grid);



%  NUMERICALLY SOLVE THE SYSTEM
%  YOU NEED TO NUMERICALLY SOLVE THE EQUATIONS HERE ***
%  For this code to work you will need to use t and y as output variables.
%  Also, t must be N rows x 1 column, and y must be N rows by M columns
%  If you are using python, you might want to check out numpy.multiply for
%  element-wise multiplication
t_span = [0, 200];
init_cond = [0.1, 0.1];
[t, y] = ode45(@(t,y) fitz_hugh_starter(t,y,a, epsilon, gam, current), t_span, init_cond);

% ------------------------- OPTIONAL -------------------------------------
%  Make vectors out of the numerical solution.  These should be the
%  difference between subsequent points, so diff should work here. The idea
%  is that this will plot out vectors ON TOP of the trajectory, so we can
%  see its direction
v_vec = diff(y(:,1));
w_vec = diff(y(:,2));
% -----------------------------END OPTIONAL--------------------------------

%  ----------------------PLOT THE RESULTS OUT------------------------------
%  PLOT THE PHASE PLANE
figure
hold on

%  Plot the nullclines  - REPLACE THE NAN'S HERE ***
plot(v_values, w_for_v, 'linewidth', 2)
plot(v_values, w_for_w, 'linewidth', 2)

%  Plot the trajectory  - REPLACE THE NAN'S HERE ***
plot(y(1:end-1,1), y(1:end-1,2), 'linewidth', 2)

%  Plot arrows on top of the trajectory so we can see what direction it
%  moves in
quiver(y(1:end-1,1), y(1:end-1,2), v_vec, w_vec, 'linewidth', 2)

%  Plot a direction field on the phase plane
quiver(V_grid,W_grid, dv, dw, 'linewidth', 2)

%  Plot the equilibria - REPLACE THE NAN'S HERE ***
plot(v_roots, w_roots, 'ko', 'markersize', 10, 'markerfacecolor', 'k')

%  Set the x and y limits
xlim([-0.55, 1.55]);
ylim([-0.55, 1.55]);

%  Plot out black lines for the "x" and "y" axes
plot(xlim, zeros(size(xlim)), 'k', 'linewidth', 1)
plot(zeros(size(ylim)), ylim, 'k', 'linewidth', 1)

%  Label Everything
xlabel('v')  %   ***
ylabel('w')  %   ***
title('Phase Plane')

%  Set the axis constraints
axis equal


%  PLOT THE NUMERICAL SOLUTION AGAINST TIME
figure
hold on

%  Plot the solutions CURVE ***
%  REPLACE THE NAN'S HERE
plot(t, y, 'linewidth', 2)

%  Label everything
xlabel('t')
ylabel('v,w of t')  %***
title('Numerical Solution')
legend('v', 'w') %***


