
>>Prism Folder

1) prism_case_1.txt: example of an input file containing the coordinates of the prismatic source
2) topology_prism.txt: the common source topology file for all considered 28 special cases
3) Cases_Prism_Coords.txt: includes the prism's coordinates for all considered 28 cases
4) Results_Prism_Cases.txt: output from function GPolyhedron.m for the aforementioned data

Example for running case 1 for the prism:
[V,Vx,Vy,Vz,Vij]=GPolyhedron(6.67259e-11,2670,'topology_prism.txt','prism_case_1.txt');



>>Asteroid Eros Folder

1) Eros_Coords.txt: coordinates of asteroid Eros in [km], with respect to its center of mass
2) eros_case_1.txt, eros_case_2.txt, etc.: input files containing the coordinates of asteroid 
   Eros in [m] for the 15 examined cases, with the origin of the coordinate system located always 
   at the computation point
3) topology_eros.txt: the common source topology file for all 15 cases
4) Results_Eros_Cases.txt: output from function GPolyhedron.m for the aforementioned data

Example for running case 1 for asteroid Eros:
[V,Vx,Vy,Vz,Vij]=GPolyhedron(6.67259e-11,2670,'topology_eros.txt','eros_case_1.txt');