% compute diffusion coefficient
load '../data/e_grid'
H1 = load('../data/H1.mat');
O16 = load('../data/O16.mat');
%O16 = load('studsvik.mat');

% allocate vars
ng = 70;
fluxx = zeros(ng,1);
J = zeros(ng,1);
enep = ene(1:ng);
a = 0.988;
b = 2.249;
chi = @(E) exp(-E/a).*sinh(sqrt(b*E));
M = zeros(ng);
D = zeros(ng);
Q = zeros(ng,1);
S = zeros(ng,1);
B2 = 0.0001;

% specify number densities 
NH1 = 4.9457e-02;
NO16 = NH1/2;
NH1=0;

% create macro xs
sigt = NH1*H1.sigt + NO16*O16.sigt;
sigs0 = NH1*H1.sigs0 + NO16*O16.sigs0;
sigs1 = NH1*H1.sigs1 + NO16*O16.sigs1;

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
                M(g,h) = M(g,h) + sigt(g);
            end
            
        end
        Q(g) = quad(chi,Elow(g),Ehi(g)) + sqrt(B2)*J(g);
    end
    
    % solve for flux
    fluxx = M\Q;
    
    % create diff matrices
    for g = 1:ng
        
        for h = 1:ng
       
            D(g,h) = -sigs1(h,g);
            if (g == h)
                D(g,h) = D(g,h) + sigt(g);
            end
            
        end
        S(g) = B2*fluxx(g)/3;
    end
    
    % solve for diff coefs
    J = D\S;    
   
end

for g = 1:ng
    
    inscatter = 0;
    
    for h = 1:ng
        inscatter = inscatter + sigs1(h,g)*J(h);
    end
    
    diff(g) = 1/(3*(sigt(g) - inscatter/J(g)));
end

trans = 1./(3*diff);
rat = trans./sigt;
semilogx(enep,rat,'b.-')