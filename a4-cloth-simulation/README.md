## Background

### Generalized Coordinates and Velocities

First, choose basis, or shape functions with which to approxiamte functions on triangular mesh. A trianle as 3 nodes the approximations become:

$$
f(\mathbf Y) = \sum_{i=0}^2 f_i \phi_i(\mathbf Y)
$$

where $\phi_i$ are the `barycentric coordinates` for a 2D triangle and $\mathbf Y \in \R^2$ is the 2D coordinate in the undeformed space.

However, cloth is really a thin volumetric construct, of which our triangle only represents a small part. We might need information lying slightly off the triangle surface. To account for this we will need to modify our FEM model a bit. First, let's assume our triangle is actually embedded in a 3D undeformed space $\mathbf X \in \R^3$ . Let's try and and build and appropriate mapping from this space to the world space.

Given any point $\mathbf X$ in the undeformed space, we can compute the barycentric coordinates of the nearest point on our triangle by solving

$$
\begin{bmatrix}
\phi_1(\mathbf X) \\
\phi_2(\mathbf X)
\end{bmatrix} = (T^TT)^{-1}T^T(\mathbf X - \mathbf X_0)
$$

where $T = \begin{bmatrix}(\mathbf X_1 - \mathbf X_0)& (\mathbf X_2 - \mathbf X_0)]\end{bmatrix}$ is a matrix of edge vectors. We use the constraint $\phi_0 + \phi_1 + \phi_2 = 1$ to reconstruct $\phi_0$. This equation finds the barycentric coordinates of the nearest point on the triangle to $\mathbf X$ in a least squares fashion.

The error in this least squares solve will be orthogonal to the column space of $T$, our triangle. For any point we can work out its offset from the triangle by computing $(\mathbf X - \mathbf X_0)^TN$, where $\mathbf X_0$ is the first vertex of our triangle and $N$ is the undeformed surface normal of the triangle. Because our triangle has a constant normal, we don't need to worry about where we compute it, which makes this all very convenient.

Let's assume that our point $\mathbf X$ maintains a constant offset from the triangle when deformed. This implies we can reconstruct the world space position by offsetting our point the same distance along the world space normal $n$. This gives us the following mapping from reference space to world space:

$$
x(\mathbf X) = \sum_{i=0}^2x_i\phi_i(\mathbf X) + (\mathbf X - \mathbf X_0)^TN \cdot n(\mathbf x_0, \mathbf x_1, \mathbf x_2)
$$

Now we can choose the generalized coordinates ($q \in \mathbb R^9$) to be the stacked vector of vertex positions, which defines the generalized velocities as the stacked $9D$ vector of per-vertex velocities.

### Deformation Gradient

In this assignment we will be able to avoid these more complicated solution due to our particular choice of undeformed to world space mapping which allows us to directly compute a $3 \times 3$ deformation gradient as:

$$
F = \frac{\partial {\mathbf x}}{\partial {\mathbf X}} = \begin{pmatrix}
 \mathbf x_0 &  \mathbf x_1 &  \mathbf x_2  &\mathbf n
\end{pmatrix} \begin{pmatrix}
-1^T(T^TT)^{-1}T^T \\
(T^TT)^{-1}T^T\\
N^T
\end{pmatrix}
$$

### Kinetic Energy

Armed with the generalized velocities, the formula for the per-triangle kinetic energy is eerily similar to that of assignment 3. It's an integral of the local kinetic energy over the entire triangle, multiplied by the thickness of the cloth, $h$. For this assignment you are free to assume the thickness of the cloth is $1$.

$$
T_{triangle} = \frac{1}{2}\dot q^T(h \begin{pmatrix}
\phi_0\phi_0 I & \phi_0\phi_1 I & \phi_0\phi_2 I \\
\phi_1\phi_0 I & \phi_1\phi_1 I & \phi_1\phi_2 I \\
\phi_2\phi_0 I & \phi_2\phi_1 I & \phi_2\phi_2 I
\end{pmatrix}) \dot q
$$

and can be compute analytically using a symbolic math package. The per-element mass matrices for every cloth triangle can then be assembled into the mass matrix for the entire mesh.

### Potential Energy

For this assignment we will use a different type of material model to describe the elastic behaviour of the cloth. This is motivated by the fact that cloth is typically very resistant to stretching. For these materials, a linear stress-strain relationship is often desirable. Unfortunately, cloth triangles also rotate a lot (every time they fold-over for instance). Rotations are **NOT** linear and so a purely linear relationship will suffer from severe artifacts. To avoid this we will build a material model that only measures the in plane deformation of the cloth via its principal stretches.

#### Principal Stretches

Recall that in the previous assignment we used the right Cauchy strain tensor ($F^TF$) to measure deformation and the rationale for using this was that it measures the squared deformed length of an arbitrary, infinitesimal line of material, $d\mathbf X$. In other words,$|d\mathbf x|^2 = d\mathbf X^T(F^TF)d\mathbf X$ . Because $F$ is symmetric and positive definite, we can perform an eigendecomposition such that $F^TF=UVU^T$ where $V$ is the orthogonal matrix of eigenvectors and is the diagonal matrix of eigenvalues. This means we can think of this squared length as $|d\mathbf x|^2 = \hat{d\mathbf X}^T\hat{d\mathbf X}$ where $\hat{d\mathbf X} = U^Td\mathbf X$. In other words, if we transform $d\mathbf X$ just right, its deformation is completely characterized by $V$.

$V$ are the eigenvalues of $F^TF$ and also the squared singular values of $F$. We call these singular values of $F$ the principal stretches. They measure deformation independently of the orientation (or rotation/reflection) of the finite element.

#### Linear Elasticity without the Pesky Rotations

Now we can formulate a linear elastic model using the principal stretches which "filters out" any rotational components. Much like the Neohookean material model, this model will have one energy term which measures deformation and one energy term that tries to preserve volume (well area in the case of cloth). We already know we can measure deformation using the principal stretches. We also know that the determinant of $F$ measures the change in volume of a 3D object. In the volumetric case this determinant is just the product of the principal stretches.

$$
\psi(s_0, s_1, s_2) = \mu\sum_{i=0}^2(s_i - 1)^2 + \frac{\lambda}{2}(s_0 + s_1 + s_2 - 3)^2
$$

where $\lambda$ and $\mu$ are the material properties for the cloth. The first term in this model attempts to keep $s_0$ and $s_1$ close to one (limiting deformation) while the second term is attempting to preserve the volume of the deformed triangle (it's a linearization of the determinant). This model is called **co-rotational linear elasticity** because it is linear in the principal stretches but rotates with each finite element. When we use energy models to measure the in-plane stretching of the cloth (or membrane), we often refer to them as membrane energies.

#### The Gradient of Principal Stretch Models

The strain energy density for principal stretch models, like the one above, are relatively easy to implement and understand. This is a big reason we like them in graphics. We'll also see that the gradient of this model (needed for force computation) is also pretty easy to compute.
Really, the derivative we need to understand how to compute is $\frac{\partial \psi}{\partial F}$. Once we have this we can use $\frac{\partial \psi}{\partial q}$ to compute the gradient wrt to the generalized coordinates. Conveniently, we have the following for principal stretch models.

$$
\frac{\partial \psi}{\partial F} = U \underbrace{\begin{bmatrix}
\frac{\partial \psi}{\partial s_0} & 0 & 0 \\
0 & \frac{\partial \psi}{\partial s_1} & 0 \\
0 & 0 & \frac{\partial \psi}{\partial s_2}
\end{bmatrix}}_{dS} V^T
$$

where $F = USV^T$ is the singular value decomposition.

#### The Hessian of Principal Stretch Models

Unfortunately, the gradient of the principal stretch energy is not enough. That's because our favourite implicit integrators require second order information to provide the stability and performance we crave in computer graphics. This is where things get messy. The good news is that, if we can just compute $\frac{\partial \psi}{\partial F \partial F}$ then we can use $\frac{\partial F}{\partial q}$ to compute our Hessian wrt to the generalized coordinates (this follows from the linearity of the FEM formulation wrt to the generalized coordinates). This formula is going to get ugly so, in an attempt to make it somewhat clear, we are going to consider derivatives wrt to single entries of $F_{ij}$, denoted . In this context we are trying to compute

$$
\frac{\partial}{\partial F_{ij}} \frac{\partial \psi}{\partial F} = \frac{\partial U}{\partial F_{ij}}dSV^T + Udiag(ds_{ij})V^T + UdS\frac{\partial V}{\partial F_{ij}}^T
$$

Here $diag()$ takes a $3\times 1$ vector as input and converts it into a diagonal matrix, with the entries of the matrix on the diagonal. In our case, we define $ds$ as:

$$
ds_{ij} = \begin{bmatrix}
\frac{\partial^2 \psi}{\partial s_0^2} & \frac{\partial^2 \psi}{\partial s_0 \partial s_1} & \frac{\partial^2 \psi}{\partial s_0 \partial s_2} \\
\frac{\partial^2 \psi}{\partial s_0\partial s_1} & \frac{\partial^2 \psi}{\partial s_1^2} & \frac{\partial^2 \psi}{\partial s_1 \partial s_2}  \\
\frac{\partial^2 \psi}{\partial s_0\partial s_2} & \frac{\partial^2 \psi}{\partial s_1 \partial s_2} & \frac{\partial^2 \psi}{s_2^2}
\end{bmatrix} \begin{bmatrix}
\frac{\partial s_0}{\partial F_{ij}} \\
\frac{\partial s_1}{\partial F_{ij}} \\
\frac{\partial s_2}{\partial F_{ij}}
\end{bmatrix}
$$

If we define the $svd$ of a matrix as $F=USV^T$, this code returns $\frac{\partial U}{\partial F} \in \mathcal R^{3\times3\times3\times3}$, $\frac{\partial V}{\partial F} \in \mathcal R^{3\times3\times3\times3}$ and $\frac{\partial U}{\partial S} \in \mathcal R^{3\times3\times3\times3}$.

Yes this code returns 3 and four dimensional tensors storing this quantities, yes I said never to do this in class, consider this the exception that makes the rule. The latter two indices on each tensor are the and indices used in the formula above.

The hardest part of implementing this gradient correctly is handling the $SVD$ terms. These gradients have a different form based on whether your $F$ matrix is square or rectangular. This is one big reason that the $3\times3$ deformation gradient we use in this assignment is desirable. It allows one to use the same singular value decomposition code for volumetric and cloth models.

### Collision Detection with Sphere

To make this assignment a little more visually interesting, you will implement simple collision detection and resolution with an analytical sphere. **Collision Detection** is the process of detecting contact between two or more objects in the scene and **Collision Resolution** is the process of modifying the motion of the object in response to these detected collisions.

For this assignment we will implement per-vertex collision detection with the sphere. This is as simple as detecting if the distance from the center of the sphere to any vertex in your mesh is less than the radius of the sphere. If you detect such a collision, you need to store an **index** to the colliding vertex, along with the outward facing **contact normal($n$)**. In this case, the outward facing contact normal is the sphere normal at the point of contact.

### Collision Resolution

The minimal goal of any collision resolution algorithm is to prevent collisions from getting worse locally. To do this we will implement a simple velocity filter approach. Velocity filters are so named because the "filter out" components of an objects velocity that will increase the severity of a collision. Given a vertex that is colliding with our sphere, the only way that the collision can get worse locally is if that vertex moves into the sphere. One way we can check if this is happening is to compute the projection of the vertex velocity onto the outward facing contact normal ($n^T\dot q_i$, $i$ selects the $i^{th}$ contacting vertex). If this number is $>0$ we are OK, the vertex is moving away from the sphere. If this number is $<0$ we better do something.

The thing we will do is project out, or filter out, the component of the velocity moving in the negative, normal direction like so:

$$
\dot q_i^{filtered} = \dot q_i - nn^T\dot q_i
$$

This "fixes" the collision. This approach to collision resolution is fast but for more complicated scenes it is fraught with peril. For instance it doesn't take into account how it is deforming the simulated object which can cause big headaches when objects are stiff or rigid. We'll see a cleaner mathematical approach to content in the final assignment of the course.

## Appendix

- Heron's formula for area
  If $a$, $b$, and $c$ are the lengths of the three sides of a triangle, and $s$ is the semi-perimeter (half of the perimeter), then Heron's formula for the area ($A$) is given by:
  $$
  A= \sqrt{s⋅(s−a)⋅(s−b)⋅(s−c)}
  $$
  where, the semi-perimeter ($s$) is calculated as:
  $$
  s = \frac{a+b+c}{2}
  $$
