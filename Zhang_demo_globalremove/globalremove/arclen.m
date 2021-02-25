function l=arclen(XI,XJ)
s=referenceSphere('Unit Sphere');
l=distance('gc',XI*180/pi,XJ*180/pi,s);
