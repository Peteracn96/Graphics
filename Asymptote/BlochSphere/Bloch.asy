import three;
import math;

settings.render    = 10; //quality
settings.outformat = "pdf"; //output

size(200);

draw(unitsphere,gray+opacity(0.4));
currentprojection=orthographic((4,2,2));

real r=1.4;
pen p=black;
draw(Label("$x$",1),O--2.0*X,p,Arrow3);
draw(Label("$y$",1),O--r*Y,p,Arrow3);
draw(Label("$z$",1),O--r*Z,p,Arrow3);
label("$\rm O$",(0,0,0),-1.5Y-X);

real theta = 0.7;
real phi = 0.95;

real x = cos(phi)*sin(theta);
real y = sin(phi)*sin(theta);
real z = cos(theta);

triple P = X + Y + Z;

triple Psi = (x,y,z);

draw(Label("$|\psi\rangle$",1),O--Psi,p);
dot(Psi,black);
dot(Label("$|0\rangle$",1,black),(0,0,1),align=(0.,0.8,0.),white);
dot(Label("$|1\rangle$",1,black),(0,0,-1),white);

draw("$\theta$",arc(O,0.4*Z,0.4*Psi),align=1.0*dir(theta/2,phi),p);
draw("$\phi$",arc(O,0.3*X,0.3*(x,y,0)),align=-0.5*dir(theta,phi/2),p);

draw(O--(x,y,0),dashed+black);
draw(Psi--(x,y,0),dashed+black);