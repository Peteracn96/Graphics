import atoms;

settings.prc       = false;
settings.render    = 10; //quality
settings.outformat = "png"; //output 

unitsize(0.04cm);

size(14cm,0);

// Colors for the atoms and bonds
pen Mcolor = heavyred;
pen Xcolor = blue;
pen bond = gray(0.4);

// Atom surfaces
surface metal = scale3(0.01)*unitsphere; // 70 pm
surface calchogenide = scale3(0.08)*unitsphere; // 25 pm

// Lattice parameters
real a =  3.153*3; // a parameter in pm
real bond_length =  a; // in pm
real c = 1.0; // off-plane height

// Number of unit cells
int Nx = 5;
int Nxplus = Nx + 1;

int Ny = 6;
int Nyminus = Ny - 1;

// Bravais lattice and motif vectors
triple R_1 = (sqrt(3)*a/2,1.5a,0);
triple R_2 = (-sqrt(3)*a/2,1.5a,0);


triple delta_1 = (sqrt(3)*bond_length/2,0.5bond_length,0*c);
triple delta_2 = (-sqrt(3)*bond_length/2,0.5bond_length,0*c);
triple delta_3 = (0,bond_length,c);

triple delta_4 = (sqrt(3)*bond_length/2,0.5bond_length,-c);
triple delta_5 = (-sqrt(3)*bond_length/2,0.5bond_length,-c);
triple delta_6 = (0,bond_length,-0c);

real bond_radius = 1.0;
real radius_scale = 0.2;

Atom[] M_1_array;
Atom[] X_1_array;

Atom[] M_2_array;
Atom[] X_2_array;

// Generates and stores the lattice sites

for(int x = 0; x < Nx; ++x) {
    for(int y = 0; y < Nyminus; ++y) {

        triple positionM = x*(R_1 + R_2) + y*(R_1 - R_2);
        triple positionXtop = positionM + delta_2 - 0*delta_6;
        triple positionXbottom = positionXtop - 2*c*Z;

        int index = x*Nyminus + y;

        Atom M_1 = Atom("Mo", positionM);

        M_1.color = Mcolor;

        M_1_array[index] = M_1;

        Atom X_2_top = Atom("S", positionXtop + R_1);

        X_2_array[index] = X_2_top;
    }
}


for(int x = 0; x < Nxplus; ++x) {
    for(int y = 0; y < Ny; ++y) {

        triple positionM = x*(R_1 + R_2) + y*(R_1 - R_2);
        triple positionMshifted = positionM + R_2;

        triple positionXtop = positionM + delta_2 - 0*delta_6;
        triple positionXbottom = positionXtop - 2*c*Z;

        int index = x*Ny + y;

        Atom X_1_top = Atom("S", positionXtop);

        X_1_array[index] = X_1_top;


        Atom M_2 = Atom("Mo", positionMshifted);

        M_2_array[index] = M_2;

        M_2.color = Mcolor;

    }
}


// Draws the atoms

for(int x = 0; x < Nx; ++x) {
    for(int y = 0; y < Nyminus; ++y) {
        
        int index = x*Nyminus + y;

        M_1_array[index].draw();
    
        X_2_array[index].draw();
    }
}

for(int x = 0; x < Nxplus; ++x) {
    for(int y = 0; y < Nyminus; ++y) {
        
        int index = x*Nyminus + y;

        M_2_array[index].draw();

        X_1_array[index].draw();
    }
}

// Draws the bonds

int jump = 0;
for(int y = 0; y < Nyminus; ++y) {
    for(int x = 0; x < Nx; x = x + 1) {

        int index_1 = x + y*Nx;

        int index_2 = index_1 + jump;
    
        Bond(M_1_array[index_1],X_1_array[index_2]).draw(radius=bond_radius);
    }

    ++jump;
}


//draw axes for reference
draw(O--10*a*X, black, arrow = Arrow3(TeXHead2, black), L=Label("$x$", position = EndPoint, align=W)); //x-axis
draw(O--20*a*Y, black, arrow = Arrow3(TeXHead2, black), L=Label("$y$", position = EndPoint, align=W)); //y-axis
draw(O--2*a*Z, black, arrow = Arrow3(TeXHead2, black), L=Label("$z$", position = EndPoint, align=W)); //z-axis

currentprojection = perspective(12*a*X + 2.5*a*Y + 4*a*Z);
