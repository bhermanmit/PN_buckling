% compute diffusion coefficient
load '../data/e_grid.mat'
load '../data/H1.mat'

% allocate vars
ng = 70;
fluxx = zeros(ng,1);
enep = ene(1:ng);
a = 0.988;
b = 2.249;
chi = @(E) exp(-E/a).*sinh(sqrt(b*E));
M = zeros(ng);
D = zeros(ng);
Q = zeros(ng,1);
S = zeros(ng,1);
B2 = 0.0001;

% estimate initial diffusion coefficients
diff = 1./(3*sigt);

% begin iteration loop
for i = 1:5
    
    M(:,:) = 0;
    D(:,:) = 0;
    Q(:) = 0;
    S(:) = 0;
    
    % create flux matrices
    for g = 1:ng
        
        for h = 1:ng
       
            M(g,h) = -sigs0(h,g);
            if (g == h)
                M(g,h) = M(g,h) + sigt(g) + diff(g)*B2;
            end
            
        end
        Q(g) = quad(chi,Elow(g),Ehi(g));
    end
    
    % solve for flux
    fluxx = M\Q;
    
    % create diff matrices
    for g = 1:ng
        
        for h = 1:ng
       
            D(g,h) = -sigs1(h,g)/sigt(g)*fluxx(h)/fluxx(g);
            if (g == h)
                D(g,h) = 1 + D(g,h);
            end
            
        end
        S(g) = 1/(3*sigt(g));
    end
    
    % solve for diff coefs
    diff = D\S;    
   
end

trans = 1./(3*diff);
rat = trans./sigt;
semilogx(enep,rat,'.')