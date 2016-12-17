Main file: C3D8.m

The code requires the following:

- Node info (in csv file ): Node ID, X-Coordinate, Y-Coordinate, Z-Coordinate
- Element info (in csv file ): Element ID, Material ID, Element connectivity matrix
- Material info (in txt file ): Material ID, Elastic Modulus, Poisson’s Ratio
- Boundary Conditions (To be edited in the code):
	– Dirichlet BC: Node ID, DOF, Value
	– Neumann BC: Element ID, Nodes, DOF, Value

With these information, the code should be able to handle most 3D elastostatic problems usinng C3D8 elements.
