function [xCoeff, errorR, numPhi1, numPhi2, timing] = numericSSH_LM(varargin)
    %Numerically computes SSH-type interaction coefficients for a 1-D
    %topological insulator system. NOTE: Expect real eigenvalues
    %
    %To call, try [xCoeff, error, numPhi1, numPhi2, timing] = symbolicSSH_LM(length, tolerance, potential, ...)
    %
    %xCoeff will return the interaction coefficients in the following format:
    %xCoeff = [a_1;b_1;b_2;a_2;b_3;b_4;...] 
    %
    %Format is [onsite, offsite right, offsite left, onsite (next nearest neighbor), offsite right,
    %offsite left,...]
    %
    %error gives two-norm squared difference between continuous
    %bands and numerical bands
    %
    %numPhi1 and numPhi2 are the numerical bands (first and second,
    %respectively)
    %
    %timing gives the time it took for the optimization algorithm to run
    %(from start of symbolic computations to end of while loop. Does
    %not include time taken for user inputs)
    %
    %----------------------------------------------------------------------
    %Inputs are read in the following way:
    %numericSSH_LM(length, tolerance, potential, first_band*, second_band*)
    %
    %length is period of the lattice,
    %
    %tolerance is minimum relative difference in successive iterations of the coefficients to kick out of
    %the optimization algorithm. That is, the algorithm stops if ||xOld - xNew||/||xOld|| < tolerance
    %
    %potential is the potential of the system
    %
    %first_band and second_band are the first and second band data inputs.
    %*This is an optional input.
    %
    %NOTE: At least three inputs are required. If you have a potential, call the function like
    %numericSSH_LM(length, tolerance, potential)
    %
    %If you have band data that you want to fit, call the function like
    %numericSSH_LM(length, tolerance, anything_you_want_here, first_band, second_band)
    %
    %It doesn't matter what you put in the potential input spot in this
    %case, it's ignored. If you input potential and bands, the potential
    %will be ignored and the band data will be used for fitting.
    %
    %Justin T. Cole, Michael J. Nameika, 2024
    

        format long

        %Throws an error if the program is called with an incorrect number of
        %inputs
        if (length(varargin) < 3 || length(varargin) == 4 || length(varargin) > 5)
            error('Must enter 3 or 5 input arguments. Type ''help numericSSH_LM'' for additional details.')
        end

        %Save the length of the lattice and given tolerance
        L = varargin{1};
        tol = varargin{2};


        %corrects length if input was negative
        if (L < 0)
            L = abs(L);
        end
    
        %Number of k points. Sets a default value for the number of k
        %points if working with potential, 
        if length(varargin) < 4
            P = 64;
        else 
            P = length(varargin{4});
        end


       if (length(varargin) == 3)
            %Grabs potential
            potential = varargin{3};

            %If the potential is saved as a row vector, correct it to a
            %column vector (algorithm requires column vector to work).
            if (isrow(potential))
                potential = potential.';
            end

            %grab bands from continuous problem (lambda1 and lambda2)
            [k_int, lambda1, lambda2] = cont_bands(potential, L, P);
            
            
            %Use maximially localized wannier functions to determine if the
            %potential orientation is topological or nontopological and initialize
            %the initial guess respectively.
            wannierCenter = localized_Wannier(potential,L);
    
            %Get the center of mass of the left and right potentials
            [leftCenter, rightCenter] = potentialCenter(potential, L);
    
            %if the wannier Center is between the peaks, nontopological,
            %otherwise, topological
            if (leftCenter < wannierCenter && wannierCenter < rightCenter)
                checkTop = 0;
            else
                checkTop = 1;
            end
        
       end

       %Case where numerical bands are given
       if (length(varargin) == 5)
            checkTop = -1;
            if (length(varargin{4}) ~= length(varargin{5}))
                error('Band data must have the same number of points.')
                
            end

            %Grab first and second bands
            lambda1 = varargin{4};
            lambda2 = varargin{5};
            
            %Changes band data to column vectors if given as row vectors
            if isrow(lambda1)
                lambda1 = lambda1.';
            end
            if isrow(lambda2)
                lambda2 = lambda2.';
            end

            %Ask user if band data is topological or nontopological 
            while (checkTop ~= 0 && checkTop ~= 1)
                promptTop = "Does your band data correspond to a topological or nontopological potential? (Enter 1 for topological and 0 for nontopological): ";
                checkTop = str2double(input(promptTop, 's'));
            end

            %Build k-space vector and remove last point for periodic B.C.s
            k_int = linspace(-pi/L, pi/L, P+1);
            k_int(end) = [];
            k_int = k_int.';
       end
        
       
        %Ask user to input number of coefficients
        prompt = "Enter the number of SSH coefficients (must be divisible by 3): ";
        n = input(prompt);
        
    
        
        %Throws an error if number of coefficients is not divisible by 3
        while (rem(n, 3) ~= 0 || n <= 0)
            disp('ERROR - The number of coefficients must be postive and divisible by 3')
            n = input(prompt);
        end
        
        %Begin program timing. Since user input is done, all time spent
        %from here is in the numerical computations
        tic

        %Prompt the user that the program is starting symbolic
        %computations. These updates will enable the user to see
        %approximately how far along the program is
        fprintf('\n%s\n\n', 'Entering variable initialization and symbolic computations. Please be patient, this may take a while :)')
       
        %initialize initial guess vector
        xVec = zeros(n,1);

        %Initializes counting variable for tracking total number of
        %iterations
        count = 1; 
        
        %initializes errorHold for holding the 2-norm of the residual
        errorHold = [];


        %if the potential orientation is nontopological, start with a
        %nontopological initial guess. Otherwise, use a topological initial
        %guess
        if (checkTop == 0) 
            %nontopological initial guess (c > d)
            xVec(1) = 10;
            xVec(2) = 5;
            xVec(3) = 2;
        else 
            %topological initial guess (c < d)
            xVec(1) = 10;
            xVec(2) = 2;
            xVec(3) = 5;
        end
        
        %Sets the max number of iterations
        max_count = 400;
           
        %Weights for biasing a specific band. w1 weighs first band, w2 weighs
        %second band
        %By default, bands are equally weighted with weight 1.
        w1 = 1;
        w2 = 1;
        
        
        %------------------------------------------------------------------
        %---------------SYMBOLIC VARIABLE INITIALIZATION-------------------

        %Symbolic placeholder variables for interaction coefficients
        A = sym('A', [n 1]);

        %symbolic variable for taking in numerical data
        x = sym('x');
        syms phi(x) phi_2(x)
       
        %---------------END SYMBOLIC VARIABLE INITIALIZATION---------------
        %------------------------------------------------------------------
       
        %------------------------------------------------------------------
        %------------BEGIN OBJECTIVE FUNCTION INITIALIZATION---------------

        %Get objective functions for optimization algorithm
        [phi(x), phi_2(x), ~] = eigenPhiNum(n,L);


        %Symbolically differentiates objective functions for top and 
        %bottom bands; Jac1 and Jac2 are used as placeholders for
        %substituting the values of xVec in for the symbolic variable A
        Jac1 = simplify(gradient(phi, A));
        Jac2 = simplify(gradient(phi_2, A));
   
        %Initialies matrices to hold the numerical values of J1 and J2,
        %respectively
        Jacobian1_double = zeros(P,n);
        Jacobian2_double = zeros(P,n);
        
        %Initializes Jacobian1 and Jacobian2 to hold Jac1 and Jac2 with
        %updated coefficient values, respectively
        Jacobian1(x) = sym(zeros(n,1));
        Jacobian2(x) = sym(zeros(n,1));
        
        %Initializes phiEval1 and phiEval2 to hold phi1 and phi2 with
        %updated coefficient values, respectively
        phiEval1(x) = sym(0);
        phiEval2(x) = sym(0);
        
        %-------------END OBJECTIVE FUNCTION INITIALIZAITON----------------
        %------------------------------------------------------------------

        %Initialize xDiff with arbitrary value to ensure while loop doesn't
        %kick out on first iteration
        xDiff = 1;

        %xHold holds the previous iteration's xVec update for calculating
        %the difference between xVec each iteration (this is used in computing xDiff
        %in the while loop)
        xHold = xVec;

        %rho is the initial Levenberg-Marquardt (LM) parameter
        rho = 5;

        %n by n identity matrix
        I = eye(n,n);

        %LM update parameter
        nu = 1.1;

        errOld = 10;

        %Update the user that variable initializations are complete and the
        %optimization algorithm begins.
        fprintf('%s\n\n', 'Variable initialization and symbolic computation complete.')
        fprintf('%s\n\n', 'Entering optimization algorithm.')

       
        
        while (count < max_count && norm(xDiff) > tol)
           
            %Substitutes the xVec values in for the symbolic variables A in
            %Jacobian1 and Jacobian2
            %NOTE: These are symbolic function vectors
            Jacobian1(x) = subs(Jac1, A, xVec);
            Jacobian2(x) = subs(Jac2, A, xVec);

            
            %Substitutes xVec values in for symbolic variables A for
            %phiEval1 and phiEval2
            phiEval1(x) = subs(phi(x), A, xVec);
            phiEval2(x) = subs(phi_2(x), A, xVec);
            
           
            %passes in k_int from schrodinger_spec_band script and converts to
            %double for the numerics
            numPhi1 = double(phiEval1(k_int));
            numPhi2 = double(phiEval2(k_int));
            

            
            %Defines residual function for first band
            r1 = w1*(numPhi1 - lambda1);

            %Defines residual function for band band
            r2 = w2*(numPhi2 - lambda2);

            %builds total residual by stacking residual for first band on
            %residual for second band
            r = [r1;r2];

            %Here J_frst_band and J_scnd_band are set to Jacobian1(x) and
            %Jacobian2(x). The reason this is done is because indexing on
            %J_frst_band and J_scnd_band will give the elements of the
            %symbolic vector, whereas indexing on Jacobian1 and Jacobian2
            %will substitute the index for x. This is a simple way to
            %circumvent this issue
            J_frst_band = Jacobian1(x);
            J_scnd_band = Jacobian2(x);

            %Initialize placeholder symbolic variables for holding the
            %elements of the Jacobians for the first and second bands
            J1_placehold = sym(0);
            J2_placehold = sym(0);
            
            for j=1:n
                %Set the placeholder variables to be the jth Jacobian
                %element
                J1_placehold(x) = J_frst_band(j);
                J2_placehold(x) = J_scnd_band(j);

                %convert the Jacobians evaluated on the k-space grid to a
                %double data type
                Jacobian1_double(:,j) = double(J1_placehold(k_int));
                Jacobian2_double(:,j) = double(J2_placehold(k_int));
               
            end %end for loop
            
            
            %builds 'total' Jacobian by stacking Jacobian for first band on the
            %Jacobian for the second band
            J = [w1*Jacobian1_double;w2*Jacobian2_double];
            

            %Right Hand Side of least squares problem
            RHS = -J.'*r;
            
           
            %p is the direction vector for updating coefficient values (use
            %MATLAB's built in  '\'  command)
            p = (J.'*J + rho*I)\RHS;
            


            %updates xVec with the direction vector p
            xVec = xVec + p;


            %Update the diff vec
            xDiff = (xVec - xHold)/norm(xVec);
            xHold = xVec;
            
            
            %updates error
            diffErr = 1/2*norm(r,2)^2;
            
            errorHold = [errorHold, diffErr];
            
            %simple update criterion for the adaptive parameter rho. If the
            %error increases, increase rho by a factor of nu. Otherwise,
            %decrease it by a factor of nu
            if (diffErr > errOld)
                rho = nu*rho;
            else
                rho = rho/nu;
            end
            
            %update error
            errOld = diffErr;

            %increments count
            count = count + 1;

        end %end while loop

        
        %Saves final update of xVec to xCoeff 
        xCoeff = xVec;
        
        %Compuptes final least squares error
        errorR = 1/2*norm(r,2)^2;
        % errorR = norm(r,inf);
        
        %append the new computed error to errorHold
        errorHold = [errorHold, errorR];

        %Save time for program to complete
        timing = toc;
        
        %Update user that optimization algorithm has ended
        fprintf('%s\n\n', 'Algorithm complete.')

        %displays number of iterations for algorithm to complete
        count
    

        %phiEval for plotting first band
        phiEval(x) = subs(phi(x), A, xVec);
        numPhi1 = double(phiEval(k_int));

        %Create phi2 for plotting second band
        phi2(x) = phi_2(x);
        phi2Eval(x) = subs(phi2(x), A, xVec);
        numPhi2 = double(phi2Eval(k_int));


        %plots phi1 and phi2
        figure(1)
        plot(k_int, real(lambda1), 'r', k_int, real(lambda2),'b', k_int, real(numPhi1), '--', k_int, real(numPhi2), '--', 'LineWidth',2)
        xlabel('$k$', 'fontsize', 25, 'interpreter', 'latex')
        ylabel('$\lambda(k)$', 'fontsize', 25, 'interpreter', 'latex')
        title(n + " Coefficient Fit", 'fontsize', 16, 'interpreter', 'latex')
        legend('Continuous Bands', 'Continuous Bands', 'Discrete Bands', 'Discrete Bands', 'FontSize', 8, 'interpreter', 'latex')
        axis tight
        
        figure(2)
        plot(k_int, imag(lambda1), 'r', k_int, imag(lambda2),'b', k_int, imag(numPhi1), '--', k_int, imag(numPhi2), '--', 'LineWidth',2)
        xlabel('$k$', 'fontsize', 25, 'interpreter', 'latex')
        ylabel('$\lambda(k)$', 'fontsize', 25, 'interpreter', 'latex')
        title('Imaginary Bands', 'fontsize', 16, 'interpreter', 'latex')
        legend('Continuous Bands', 'Continuous Bands', 'Discrete Bands', 'Discrete Bands', 'FontSize', 8, 'interpreter', 'latex')
        axis tight
        
        figure(3)
        semilogy(errorHold, 'linewidth', 2)
        grid on
        title('Residual vs. Iterations', 'FontSize', 20, 'interpreter', 'latex')
        xlabel('Iteration Count', 'FontSize', 20, 'interpreter', 'latex')
        ylabel('Residual ($\frac{1}{2}||R||_2^2$)', 'FontSize', 20, 'interpreter', 'latex')
        hold off
      
        edgeModeDecision = '';

        %ask user if they want to plot eigenenergies
        while (isempty(edgeModeDecision) || edgeModeDecision(1) ~= 'Y' && edgeModeDecision(1) ~= 'y' && edgeModeDecision(1) ~= 'n' && edgeModeDecision(1) ~= 'N')

            edgeModeDecision = input('Would you like to calculate the edge modes? Y/N: ', 's');

            if (~isempty(edgeModeDecision))
                if (edgeModeDecision(1) == 'Y' || edgeModeDecision(1) == 'y')


                    [~, edgeEigVal] = edgeModeFunction(150, xVec);


                    edgeEigVals = diag(edgeEigVal);

                elseif (edgeModeDecision(1) == 'N' || edgeModeDecision(1) == 'n')

                else
                    disp("Please enter Y/N:")
                end %end nested if
            else
                disp("Please enter Y/N:")
            end %end if

        end %end while

        windingDecision = '';

        %ask user if they want to compute and display the winding number
        while (isempty(windingDecision) || windingDecision(1) ~= 'n' && windingDecision(1) ~= 'N' && windingDecision(1) ~= 'y' && windingDecision(1) ~= 'Y')

            windingDecision = input('Would you like to calculate the winding number? Y/N: ', 's');
            if (~isempty(windingDecision))

                if (windingDecision(1) == 'Y' || windingDecision(1) == 'y')

                    windingNum = discWindingNumCalc(xCoeff, n, L, k_int)

                elseif (windingDecision(1) == 'N' || windingDecision(1) == 'n')

                else
                    disp("Please enter Y/N:")
                end %end nested if

            else 
                disp("Please enter Y/N:")
            end %end if

        end %end while
        
        %prints a nice table with the computed interaction coefficients
        fprintf("\n%s\n", "-----------------------------------------------")
        fprintf("%s\t\t%s\n", "Interaction Term", "     Interaction Value")
        fprintf("%s\n", "-----------------------------------------------")

        xcount = 1;

        for i = 1:n
            if (i == 1)

                fprintf("%s%d%s%f\n", "x", xcount, "(self) . . . . . . . . . . ", xCoeff(i))
                xcount = xcount + 1;

            else

                fprintf("%s%d%s%f\n", "x", xcount, " . . . . . . . . . . . . . ", xCoeff(i))
                xcount = xcount + 1;

            end
        end
        fprintf("\n\n")
    
end %end main

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eig1, eig2, M] = eigenPhiNum(n,L)
%Symbolically builds the objective function for the nonlinear least squares
%process. 
%
%Uses symmetries of the system to efficiently compute the objective
%functions.
%
%Returns the objective function corresponding to the top and bottom band.

    
    A = sym('A', [n 1], 'real');
    x = sym('x', 'real');
    syms Ahold Bhold
    
    Ahold = 0;
    Bhold = 0;

    %populates Ahold for the interaction terms with a-site terms
    Ahold = Ahold + A(1);

    for j=1:n/3-1
        Ahold = Ahold + 2*A(3*j + 1)*cos(L*x*(j));
    end %end for loop

    %populates Bhold for the interaction with b-site terms
    for j = 1:n/3
        Bhold = Bhold + A(3*j - 1)*exp(1i*x*L*(j - 1));
        Bhold = Bhold + A(3*j)*exp(-1i*x*L*(j));
    end %end for loop


    Bsquare = simplify(expand(conj(Bhold)*Bhold), 'Criterion', 'preferReal', 'Steps', 50);

    %saves first band to eig1 and second band to eig2
    eig1 = Ahold + sqrt(Bsquare);
    eig2 = Ahold - sqrt(Bsquare);
    
    M = [Ahold Bhold; conj(Bhold) Ahold];

end %end eigenPhi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [v, d] = edgeModeFunction(N,X)
%edgeModeFunction returns the eigenvalues and eigenfunctions 
%for the (discrete) edge mode boundary value problem (eigenenergies). 
%
%----------Inputs------------
%N is the size of each block matrix; there will be 2N data points
%for evaluating edge modes
%Block matrix has the form [L11 L12]
%                          [L21 L22]
%Where L_{ij} is an NxN matrix 
%
%X is the vector of coefficients from the SSH problem
%
%-----------Outputs----------
%v returns the eigenfunctions, and d returns the corresponding
%eigenvalues in a 2Nx2N matrix


%Saves the total number of coefficients
xLen = length(X);

if (N < xLen)
    N = 2*xLen;
end

%Initializes 
L11 = zeros(N,N);
L12 = zeros(N,N);

%L11 holds the values of X that correspond to the A sites
%The first A site coefficient value goes on the main diagonal,
%all other coefficient values go on the next sub and super diagonal.
for i=1:xLen/3
    
    ind = i-1;
    
    if (i == 1)

        L11 = L11 + X(3*i-2)*diag(ones(N-ind,1),ind);
    else
        L11 = L11 + X(3*i-2)*diag(ones(N-ind,1),ind);
        L11 = L11 + X(3*i-2)*diag(ones(N-ind,1),-ind);
    end %end if
    
end %end for


%L22 holds the values of X that correspond to the B sites.
%Similar to L11, the first B site coefficient goes on the main
%diagonal, following coefficients are placed on the next sub and
%super diagonal
for j=1:2*xLen/3
    
    if (rem(j,2) == 0)
        
        ind = -j/2;
        
        L12 = L12 + X(3*(abs(ind)))*diag(ones(N-abs(ind),1),ind);
        
    else
        
        ind = (j-1)/2;
        L12 = L12 + X(3*ind+2)*diag(ones(N-ind,1),ind);

    end %end if
  
    
end %end for

L22 = L11;
L21 = L12';


A = [L11, L12; L21, L22];

%Saves the eigenfunctions and eigenvalues of A
[v, d] = eig(A);
eValues = diag(d);

%sorts eigenvalues in ascending order
[~, ind] = sort(real(eValues), 'ascend');
eValues = eValues(ind);


%Seperates eigenvalues corresponding to a and b sites
%(First N eigenvalues correspond to a sites, last N 
%eigenvalues correspond to b sites).
a = v(1:N, ind(N));
b = v(N+1:2*N, ind(N));
C = [];

%Saves eigenvalues in a vector in alternating order
%i.e. C = [a(1) b(1) a(2) b(2) ... a(N) b(N)]^T
for i=1:N
   C = [C; a(i); b(i)];
   
end %end for


%plots the eigenvalues 
figure()
plot(real(eValues), 'bo')
grid on
xlabel('Eigenenergy index', 'fontsize', 20, 'interpreter', 'latex');
ylabel('Eigenenergy value', 'fontsize', 20, 'interpreter', 'latex');
title("Eigenenergies for " + xLen + " Coefficients", 'fontsize', 16, 'interpreter', 'latex');


figure()
plot(abs(C), 'linewidth', 1.1)
title('Eigenstate', 'fontsize', 25, 'interpreter', 'latex')


end %end edgeModeFunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W_disc = discWindingNumCalc(X, n, L, k_int)
    %Compute the winding number of the discrete system. A winding number of
    %1 corresponds to a topologically nontrivial system, a winding number
    %of 0 corresponds to a topologically trivial system


    A = sym('A', [n 1], 'real');
    syms x
    
    %Pulls matrix to find the analytic eigenfunctions
    [~,~, Mmat] = eigenPhiNum(n,L);
    
    %Only pull out the real component of X; for large number of coefficients (~12+), small
    %imaginary parts find their way into X.
    X = real(X);
    
    
    %Saves the eigenvectors of the eigenvalue matrix
    Mmat(x) = subs(Mmat, A, X);
    
    M0 = double(Mmat(k_int(1)));
    
    
    [V0, D0] = eig(M0);
    [~, B0] = sort(diag(real(D0)), 'ascend');
    
    gammaK = zeros(length(k_int),1);
    
    for i=1:length(k_int)-1
        M1 = double(Mmat(k_int(i+1)));
        
        [V1, D1] = eig(M1);
        [~, B1] = sort(diag(real(D1)), 'ascend');
        
        gammaK(i) = V0(:,B0(1))'*V1(:,B1(1));
        V0 = V1;
        B0 = B1;
    end %end for

    %Put in the first element at the end for periodicity
    gammaK(end) = gammaK(1);


    %Numerically calculate the Zak Phase
    Z = sum(-imag(log(gammaK)));
    % W_disc = mod(round(Z/pi),2);
    W_disc = mod(Z/pi,2);

end %end function discWindingNumCalc()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function wCntr = localized_Wannier(V,L)
    %Used the algorithm from Vanderbilt (1997) for build a maximally
    %localized Wannier function. This implementation is for one band in the
    %one dimensional case. All vectors in this function are functions of k,
    %that is, each element varies along the k-space.

    %rate constant
    alpha = 0.01;


    %Number of spacial points
    Px = length(V);
    x = linspace(-L/2,L/2,Px+1)';
    x(end) = [];
    dx = x(2) - x(1);

    %number of k points
    Pk = 16;
    k = linspace(-pi/L,pi/L,Pk+1)';
    k(end) = [];

    %k spacing
    deltaK = 2*pi/L/Pk;

    %Differentiation matrices

    %first derivative
    [~, Dx] = fourdif(Px,1);
    Dx = Dx*(2*pi/L);
    
    %second derivative
    [~, Dxx] = fourdif(Px,2);
    Dxx = Dxx*(2*pi/L)^2;
    
    %Rows vary in x, columns in y, and layers correspond to different bands
    uStore = zeros(Px,Pk,2);
    eigStore = zeros(Pk,2);
    
    %Build the eigenfunctions
    for i=1:Pk
    
        %Linear Schrodinger operator
        Lop = Dxx + 2*1i*k(i)*Dx - k(i)^2*eye(Px) + diag(V,0);


        [u,d] = eig(-Lop);
        [a,b] = sort(diag(d), 'ascend');
        
        uStore(:,i,1) = u(:,b(1));
        eigStore(i,1) = a(1);

        uStore(:,i,2) = u(:,b(2));
        eigStore(i,2) = a(2);  

    end %end for

    %Bloch function of the first band
    u1 = uStore(:,:,1);
    u1 = u1/(sqrt(dx*u1(:,1)'*u1(:,1))); %rescaling

    %initialize Mnn from Vanderbilt
    M11 = zeros(Pk,1);

    %Build Mnn for first band - multiply by dx
    for i=1:Pk-1
        M11(i) = dx*u1(:,i)'*u1(:,i+1);
    end %end for

    %multiply the first part by the phase difference
    M11(end) = dx*u1(:,end)'*(exp(-1i*2*pi/L*x).*u1(:,1));
    
    %initialize r vector
    rbar = 0;

    %build r vector
    for i=1:Pk
        rbar = rbar - 1/(Pk)*imag(log(M11(i)));
    end %end for

    %initialize initial guess
    U0 = ones(Pk,1);

    %initialize q
    q = zeros(Pk,1);

    %build q
    for i=1:Pk
        q(i) = imag(log(M11(i)))+rbar;
    end %end for
   
    %for 1D, 1 band, T = q
    T = q;

    %superoperator S(T)
    s = (conj(T) + T)/(2*1i);

    %antiunitary update matrix
    deltaW = alpha*(-s);

    %update U
    U1 = U0.*exp(deltaW);
    U0 = U1;

    %initialize M11 to hold inner product of bloch waves
    M11New = zeros(Pk,1);

    %initialize rSquare
    rSquare = 0;

    %second moment
    for i=1:Pk
        rSquare = rSquare + 1/(Pk*deltaK^2)*(-2*real(log(M11(i))) + (imag(log(M11(i))))^2);
    end %end for

    %spread function
    omega0 = rSquare - (1/deltaK*rbar)^2;

    %initialize diff to be some arbitrary value to start while loop (MATLAB
    %doesn't have a do-while loop)
    diff = 100;

    %tolerance for difference in spread function between iterations
    tol = 1e-10;

    %initialize count
    count = 0;

    %set maximum number of iterations
    maxCount = 1000;

    %iteration loop test
    while (diff > tol && count < maxCount)

        %update M11
        for i=1:Pk-1
            M11New(i) = conj(U0(i))*M11(i)*U0(i+1);
        end %end for

        M11New(end) = conj(U0(end))*M11(end)*U0(1);

        %initialize first and second moments
        rbar = 0;
        rSquare = 0;

        %first and second moment
        for i=1:Pk
            rbar = rbar - 1/(Pk)*imag(log(M11New(i)));
            rSquare = rSquare + 1/(Pk*deltaK^2)*(-2*real(log(M11New(i))) + (imag(log(M11New(i))))^2);
        end %end for

        %update spread function
        omega1 = rSquare - (1/deltaK*rbar)^2;

        %update difference in successive spread functions
        diff = abs(omega1 - omega0);

        %update omega0
        omega0 = omega1;

        %update q
        for i=1:Pk
            q(i) = imag(log(M11New(i)))+rbar;
        end %end for

        %update T and S(T)
        T = q;
        s = (conj(T) + T)/(2*1i);

        %Anti-unitary update matrix
        deltaW = -alpha*s;
    
        %update U
        U1 = U0.*exp(deltaW);
        U0 = U1;

        %increment count
        count = count + 1;

    end %end while
    
    %initialize vector for holding bloch wave u corresponding to maximally
    %localized Wannier function
    u1Smooth = zeros(Px,Pk);

    %build u
    for i=1:Px
        u1Smooth(i,:) = U1.*u1(i,:).';
        
    end %end for

    %renormalize
    u1Smooth = u1Smooth/(sqrt(dx*u1Smooth(:,1)'*u1Smooth(:,1)));


    product = zeros(Pk,1);

    for i=1:Pk-1
        product(i) = dx*u1Smooth(:,i)'*u1Smooth(:,i+1);
        
    end %end for

    product(end) = dx*u1Smooth(:,end)'*(exp(-1i*2*pi/L*x).*u1Smooth(:,1));
    wCntr = -L/(2*pi)*imag(log(prod(product)));

end %end function localized_Wannier()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [k_int, lambda1, lambda2] = cont_bands(V, L, P)

    N = length(V); % # of spatial points
    
    % spectral differentiation matrix
    % 1st derivative
    [~,Dx] = fourdif(N,1);
    Dx = Dx*(2*pi/L)^1;
    
    % 2nd derivative
    [~,Dxx] = fourdif(N,2);
    Dxx = Dxx*(2*pi/L)^2;
    
    % quasimomentum values in interval [-pi/L, pi/L]
    k_int = linspace(-pi/L,pi/L,P+1)';
    k_int(end) = [];
    
    % number of bands to store
    num = 2;
    eig_store = zeros(num);
    
    for ii = 1:P
        
        % choose and fix the quasimomentum value    
        k = k_int(ii);    
    
        % construct the Schrodinger operator (after applying Bloch theory)
        L_mat = -Dxx - 2*1i*k*Dx + k^2*eye(N,N) - diag(V,0);  
    
        eigmat = eig(-L_mat);
        [~, I] = sort(real(eigmat), 'descend');
        % compute eigenvalues and sort
        eig_values = eigmat(I);
        eig_store(:, ii) = eig_values(1:num);
    
    end %end for
    
    %bands
    lambda1 = eig_store(1,:).';
    lambda2 = eig_store(2,:).';


end %end function cont_bands()

function [lftCntr, rghtCntr] = potentialCenter(V,L)
    %numerically computes the center of mass of the left and right
    %potentials
    
    x = linspace(-L/2, L/2, length(V)+1).';
    x(end) = [];
    
    %left potential + left half of domain
    vLeft = V(1:floor(length(V)/2) + 1);
    xLeft = x(1:floor(length(x)/2) + 1);

    %right potential + right half of domain
    vRight = V(floor(length(V)/2) + 1:end);
    xRight = x(floor(length(x)/2) + 1:end);
    
    %compute center of masses
    lftCntr = sum(vLeft.^10.*xLeft)/sum(vLeft.^10);
    rghtCntr = sum(vRight.^10.*xRight)/sum(vRight.^10);

end %end potentialCenter()


function [x, DM] = fourdif(N,m)
%
% The function [x, DM] = fourdif(N,m) computes the m'th derivative Fourier 
% spectral differentiation matrix on grid with N equispaced points in [0,2pi)
% 
%  Input:
%  N:        Size of differentiation matrix.
%  M:        Derivative required (non-negative integer)
%
%  Output:
%  x:        Equispaced points 0, 2pi/N, 4pi/N, ... , (N-1)2pi/N
%  DM:       m'th order differentiation matrix
%
% 
%  Explicit formulas are used to compute the matrices for m=1 and 2. 
%  A discrete Fouier approach is employed for m>2. The program 
%  computes the first column and first row and then uses the 
%  toeplitz command to create the matrix.

%  For m=1 and 2 the code implements a "flipping trick" to
%  improve accuracy suggested by W. Don and A. Solomonoff in 
%  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
%  The flipping trick is necesary since sin t can be computed to high
%  relative precision when t is small whereas sin (pi-t) cannot.
%
%  S.C. Reddy, J.A.C. Weideman 1998.  Corrected for MATLAB R13 
%  by JACW, April 2003.
 

    x=2*pi*(0:N-1)'/N;                       % gridpoints
    h=2*pi/N;                                % grid spacing
    zi=sqrt(-1);
    kk=(1:N-1)';
    n1=floor((N-1)/2); n2=ceil((N-1)/2);
    if m==0                                % compute first column
      col1=[1; zeros(N-1,1)];                % of zeroth derivative
      row1=col1;                             % matrix, which is identity

    elseif m==1                            % compute first column
      if rem(N,2)==0                         % of 1st derivative matrix
	topc=cot((1:n2)'*h/2);
        col1=[0; 0.5*((-1).^kk).*[topc; -flipud(topc(1:n1))]]; 
      else
	topc=csc((1:n2)'*h/2);
        col1=[0; 0.5*((-1).^kk).*[topc; flipud(topc(1:n1))]];
      end
      row1=-col1;                            % first row

    elseif m==2                             % compute first column  
      if rem(N,2)==0                         % of 2nd derivative matrix
	topc=csc((1:n2)'*h/2).^2;
        col1=[-pi^2/3/h^2-1/6; -0.5*((-1).^kk).*[topc; flipud(topc(1:n1))]];
      else
	topc=csc((1:n2)'*h/2).*cot((1:n2)'*h/2);
        col1=[-pi^2/3/h^2+1/12; -0.5*((-1).^kk).*[topc; -flipud(topc(1:n1))]];
      end
      row1=col1;                             % first row 

    else                                     % employ FFT to compute
      N1=floor((N-1)/2);                     % 1st column of matrix for m>2
      N2 = (-N/2)*rem(m+1,2)*ones(rem(N+1,2)); 
      mwave=zi*[(0:N1) N2 (-N1:-1)];
      col1=real(ifft((mwave.^m).*fft([1 zeros(1,N-1)])));
      if rem(m,2)==0
	row1=col1;                           % first row even derivative
      else
	col1=[0 col1(2:N)]'; 
	row1=-col1;                          % first row odd derivative
      end
    end
    DM=toeplitz(col1,row1);  
end

%------------------------END SUB FUNCTIONS---------------------------------