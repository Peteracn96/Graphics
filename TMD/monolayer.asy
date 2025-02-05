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
real scale = 1.6; //scales distance between atoms, instead of changing the radius of the atoms
real alpha = scale*3.16; // Å, distance between neighbouring atoms of the same element = lattice parameter
real a =  alpha/sqrt(3); // side of the hexagon
real bond_length = scale*2.41; // in pm
real c = scale*1.58; // off-plane height

real bond_radius = 3.0; // radius of the cylindric bonds
real radius_scale = 0.2;

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

Atom[] M_1_array;
Atom[] X_1_array;
Atom[] X_1_bottom_array;

Atom[] M_2_array;
Atom[] X_2_array;
Atom[] X_2_bottom_array;

// Generates and stores the lattice sites

for(int x = 0; x < Nx; ++x) {
    for(int y = 0; y < Nyminus; ++y) {

        triple positionM = x*(R_1 + R_2) + y*(R_1 - R_2);
        triple positionXtop = positionM + delta_2 + c*Z;
        triple positionXbottom = positionXtop - 2*c*Z;

        int index = x*Nyminus + y;

        Atom M_1 = Atom("S", positionM);

        M_1.color = Mcolor;

        M_1_array[index] = M_1;

        Atom X_2_top = Atom("S", positionXtop + R_1);
        Atom X_2_bottom = Atom("S", positionXbottom + R_1);

        X_2_top.color = Xcolor;
        X_2_bottom.color = Xcolor;

        X_2_array[index] = X_2_top;
        X_2_bottom_array[index] = X_2_bottom;
    }
}


for(int x = 0; x < Nxplus; ++x) {
    for(int y = 0; y < Ny; ++y) {

        triple positionM = x*(R_1 + R_2) + y*(R_1 - R_2);
        triple positionMshifted = positionM + R_2;

        triple positionXtop = positionM + delta_2 + c*Z;
        triple positionXbottom = positionXtop - 2*c*Z;

        int index = x*Ny + y;

        Atom X_1_top = Atom("S", positionXtop);
        Atom X_1_bottom = Atom("S", positionXbottom);

        X_1_top.color = Xcolor;
        X_1_bottom.color = Xcolor;

        X_1_array[index] = X_1_top;
        X_1_bottom_array[index] = X_1_bottom;

        Atom M_2 = Atom("S", positionMshifted);

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
        X_2_bottom_array[index].draw();
    }
}

for(int x = 0; x < Nxplus; ++x) {
    for(int y = 0; y < Nyminus; ++y) {
        
        int index = x*Nyminus + y;

        M_2_array[index].draw();

        X_1_array[index].draw();
        X_1_bottom_array[index].draw();
    }
}

// Draws the bonds

int jump = 0;
for(int y = 0; y < Nyminus; ++y) {
    for(int x = 0; x < Nx; x = x + 1) {

        int index_1 = x + y*Nx;

        int index_2 = index_1 + jump;
        int index_3 = index_2 + 1;
    
        Bond(M_1_array[index_1],X_1_array[index_2]).draw(radius=bond_radius);
        Bond(M_1_array[index_1],X_1_array[index_3]).draw(radius=bond_radius);

        Bond(M_1_array[index_1],X_1_bottom_array[index_2]).draw(radius=bond_radius);
        Bond(M_1_array[index_1],X_1_bottom_array[index_3]).draw(radius=bond_radius);
    }

    ++jump;
}

jump = 0;
for(int y = 0; y < Nyminus; ++y) {
    for(int x = 0; x < Nx; x = x + 1) {

        int index_1 = x + y*Nx;

        int index_2 = index_1 + jump;

        int index_3 = index_2 + 1;
    
        Bond(M_2_array[index_2],X_2_array[index_1]).draw(radius=bond_radius);
        Bond(M_2_array[index_3],X_2_array[index_1]).draw(radius=bond_radius);

        Bond(M_2_array[index_2],X_2_bottom_array[index_1]).draw(radius=bond_radius);
        Bond(M_2_array[index_3],X_2_bottom_array[index_1]).draw(radius=bond_radius);
    }

    ++jump;
}

jump = 0;
for(int x = 0; x < Nxplus; ++x) {
    for(int y = 0; y < Nyminus; ++y) {
        
        int index = x*Nyminus + y;
        int index_2 = index + Nx;

        Bond(X_1_array[index],M_2_array[index]).draw(radius=bond_radius);
        Bond(X_1_bottom_array[index],M_2_array[index]).draw(radius=bond_radius);
    }
}

jump = 0;
for(int y = 0; y < Nyminus - 1; ++y) {
    for(int x = 0; x < Nx; ++x) {

        
        int index = x + y*Nx;
        int index_2 = index + Nx;

        Bond(M_1_array[index_2],X_2_array[index]).draw(radius=bond_radius);
        Bond(M_1_array[index_2],X_2_bottom_array[index]).draw(radius=bond_radius);
    }
}

//draw axes for reference
draw(O--10*a*X, black, arrow = Arrow3(TeXHead2, black), L=Label("$x$", position = EndPoint, align=W)); //x-axis
draw(O--20*a*Y, black, arrow = Arrow3(TeXHead2, black), L=Label("$y$", position = EndPoint, align=W)); //y-axis
draw(O--2*a*Z, black, arrow = Arrow3(TeXHead2, black), L=Label("$z$", position = EndPoint, align=W)); //z-axis

currentprojection = perspective(12*a*X + 3.0*a*Y + 2*a*Z);
