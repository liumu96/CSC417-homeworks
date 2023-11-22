# Finite Elements Methods

## Background

## Assignement Implementation

[`phi_linear_tetrahedron.cpp`](./src/phi_linear_tetrahedron.cpp)

Evaluate the linear shape functions for a tetrahedron. This function returns a 4D vector which contains the values of the shape functions for each vertex at the world space point x (assumed to be inside the element).

```cpp

phi.setZero();

Eigen::Vector3d X0 = V.row(element(0));
Eigen::Vector3d X1 = V.row(element(1));
Eigen::Vector3d X2 = V.row(element(2));
Eigen::Vector3d X3 = V.row(element(3));

Eigen::Matrix3d T;
T << (X1 - X0), (X2 - X0), (X3 - X0);

phi.tail<3>() = T.inverse() * (x - X0);
phi(0) = 1 - phi.tail<3>().sum();

```

[`dphi_linear_tetrahedron_dX.cpp`](./src/dphi_linear_tetrahedron_dX.cpp)

Piecewise constant gradient matrix for linear shape functions. Row $i$ of the returned matrix contains the gradient of the $i^{th}$ shape function.

$$
X = \sum_{i = 0}^3 X_i\phi_i(X)
$$

$$
\begin{pmatrix}
X \\ Y \\ X \end{pmatrix} = \sum_{i = 0}^3\begin{pmatrix}
X_i \\ Y_i \\ X_i \end{pmatrix} \phi_i(X)
$$

$$
\underbrace {\begin{pmatrix}
X \\ Y \\ X \end{pmatrix}}_b = \underbrace {\begin{pmatrix}
X_1 - X_0 & X_2 - X_0 & X_3 - X_0 \\
Y_1 - Y_0 & Y_2 - Y_0 & Y_3 - Y_0 \\
Z_1 - Z_0 & Z_2 - Z_0 & Z_3 - Z_0 \\
\end{pmatrix}}_A \underbrace{\begin{pmatrix}
\phi_0(X)  \\
\phi_1(X) \\
\phi_2(X) \\
\phi_3(X)
\end{pmatrix}}_{\mathbf x}
$$

$$
\phi_0(\mathbf X) = 1 - \phi_1(\mathbf X) - \phi_2(\mathbf X) - \phi_3(\mathbf X)
$$

$$
\begin{pmatrix}
X - X_0  \\
Y - Y_0  \\
Z - Z_0  \\
\end{pmatrix} = \begin{pmatrix}
X_1 - X_0 & X_2 - X_0 & X_3 - X_0 \\
Y_1 - Y_0 & Y_2 - Y_0 & Y_3 - Y_0 \\
Z_1 - Z_0 & Z_2 - Z_0 & Z_3 - Z_0 \\
\end{pmatrix} \begin{pmatrix}
\phi_0(X)  \\
\phi_1(X) \\
\phi_2(X) \\
\end{pmatrix} = \underbrace {\begin{pmatrix}
\Delta X_1 & \Delta X_2 & \Delta X_3 \\
\Delta Y_1 & \Delta Y_2 & \Delta Y_3 \\
\Delta Z_1 & \Delta Z_2 & \Delta Z_3 \\
\end{pmatrix}}_T \begin{pmatrix}
\phi_0(X)  \\
\phi_1(X) \\
\phi_2(X) \\
\end{pmatrix}
$$

Barycentric Coordinates

$$
\begin{pmatrix}
\phi_0(X)  \\
\phi_1(X) \\
\phi_2(X) \\
\end{pmatrix} = T^{-1}(X - X_0) \\
\phi_0(\mathbf X) = 1 - \phi_1(\mathbf X) - \phi_2(\mathbf X) - \phi_3(\mathbf X)
$$

$$
\frac{d \phi}{d X} = \begin{pmatrix} -1^TT^{-1} \\ T^{-1} \end{pmatrix} \text{where } 1^T = \begin{pmatrix} 1 & 1 & 1 & 1\end{pmatrix}
$$

```cpp
dphi.setZero();

Eigen::Vector3d X0 = V.row(element(0));
Eigen::Vector3d X1 = V.row(element(1));
Eigen::Vector3d X2 = V.row(element(2));
Eigen::Vector3d X3 = V.row(element(3));

Eigen::Matrix3d T;
T << (X1 - X0), (X2 - X0), (X3 - X0);

dphi.block<3, 3>(1, 0) = T.inverse();
dphi.block<1, 3>(0, 0) = -1.0 * Eigen::RowVector3d::Ones() * dphi.block<3, 3>(1, 0);
// dphi.block<1, 3>(0, 0) = dphi.block<3, 3>(1, 0).colwise().sum() * -1.0;
```

[`psi_neo_hookean.cpp`](./src/psi_neo_hookean.cpp)

Compute the Neo-Hookean strain energy density.

$$
J = det(F) \text{ }\text{ }
$$

$$
\Psi(F) =  \frac{\mu}{2} (J^{-2/3}tr(F^T F) - 3) + \frac{\lambda}{2}(J - 1)^2
$$

```cpp
double J = F.determinant();
psi = C * (pow(J, -2.0 / 3.0) * (F.transpose() * F).trace() - 3) + D * pow(J - 1, 2.0);
```

[`dpsi_neo_hookean_dF.cpp`](./src/dpsi_neo_hookean_dF.cpp)

Compute the first Piola-Kirchoff (the gradient of the strain energy density with respect to the deformation gradient). You can use a symbolic math package to do this (but don't just look it up on the internet please).

$$
\frac{d(tr(F^TF))}{dF} = 2 F
$$

$$
\frac{d(J)}{dF} = adj(F) = F^{-T}, J = det(F)
$$

so:

$$
\Psi(F) =  C (J^{-2/3}tr(F^T F) - 3) + D(J - 1)^2, C = \frac{\mu}{2}, D = \frac{\lambda}{2}
$$

$$
\frac{d\Psi}{dF} = \frac{2CF}{J^{2/3}} - \frac{2C tr(F^T F)}{3}\frac{1}{J^{5/3}}F^{-T} + 2D(J-1)F^{-T}
$$

```cpp
Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

// Compute the Green deformation tensor
Eigen::Matrix3d FTF = F.transpose() * F;

// Compute the determinant of the deformation gradient
double J = F.determinant();

// Compute the trace of the Green deformation tensor
double trFTF = FTF.trace();

// Compute the partial derivative of the trace term
Eigen::Matrix3d d_trace_term = -2.0 / 3.0 * C * std::pow(J, -5.0 / 3.0) * trFTF * F + 2.0 * C * std::pow(J, -2.0 / 3.0) * F;

// Compute the partial derivative of the determinant term
Eigen::Matrix3d d_det_term = 2.0 * D * (J - 1) * F.inverse().transpose();

// Combine the partial derivatives and multiply by C and D
Eigen::Matrix3d term1 = C * d_trace_term;
Eigen::Matrix3d term2 = 2.0 * D * (J - 1) * F.inverse().transpose();

psi = term1 + term2;
```

[`d2psi_neo_hooken_dF2.cpp`](./src/d2psi_neo_hooken_dF2.cpp)

Compute the hessian of the strain energy density with respect to the deformation gradient. You can use a symbolic math package to do this (but don't just look it up on the internet please).

$$
\Psi(F) =  C (J^{-2/3}tr(F^T F) - 3) + D(J - 1)^2, C = \frac{\mu}{2}, D = \frac{\lambda}{2}
$$

$$
\frac{d\Psi}{dF} = \frac{2CF}{J^{2/3}} - \frac{2C tr(F^T F)}{3}\frac{1}{J^{5/3}}F^{-T} + 2D(J-1)F^{-T}
$$

For $ g = \frac{2CF}{J^{2/3}}$

$$
\frac{\partial g}{\partial F_{kl}} = \frac{2C}{J^{2/3}}|_{i=k, j=l} - \frac{4CF_{ij}F^{-T}_{kl}}{3J^{5/3}}
$$

For $h = - \frac{2C tr(F^T F)}{3}\frac{F^{-T}}{J^{5/3}}$

$$
\frac{\partial h}{\partial F_{kl}} = - \frac{4C F_{k,l}}{3}\frac{F^{-T}_{i,j}}{J^{5/3}} + - \frac{2C tr(F^T F)}{3}[-\frac{5F^{-T}_{i,j} F^{-T}_{k,l}}{3J^{8/3}} + \frac{[-(F^{-T}_{jl})(F^{-T}_{li})]}{J^{5/3}}]
$$

For $k = 2D(J-1)F^{-T}$

$$
\frac{\partial k}{\partial F_{kl}} = 2DF^{-T}_{i,j}F^{-T}_{k,l} + 2D(J-1)[-(F^{-T}_{jl})(F^{-T}_{li})]
$$

```cpp

```

[`T_linear_tetrahedron.cpp`](./src/T_linear_tetrahedron.cpp)

Compute the kinetic energy of a single tetrahedron.

```cpp

```

[`quadrature_single_point.h`](./include/quadrature_single_point.h)

Single point quadrature for a constant strain tetrahedron (CST).

```cpp
// Call the integrand function to compute the value at the given coordinates
    integrand(integrated, q, element);

// Multiply the result by the volume of the tetrahedron
integrand *= volume;
```

[`V_linear_tetrahedron.cpp`](./src/V_linear_tetrahedron.cpp)

Compute the potential energy of a single tetrahedron. Note: you will need both `psi_neo_hookean.cpp` and `quadrature_single_point.h` to do this.

$$
\Delta \mathbf x = \mathbf x(X_0 + \Delta X) - \mathbf x(X_0)
$$

$$
\Delta \mathbf x  \approx \mathbf x(X_0) + \frac{\partial \mathbf x}{\partial X}\Delta X - \mathbf x(X_0) = \underbrace{\frac{\partial \mathbf x}{\partial X}}_F\Delta X
$$

Strain: $l^2 - l_0^2$,

$l$: deformed length squared, $l_0$: rest length squared

$$
l^2 = \Delta \mathbf x^T \Delta \mathbf x \text{ }\text { } l_0^2 = \Delta X^T \Delta X
$$

$$
\Delta X^T F^T F \Delta X - \Delta X^T \Delta X = \Delta X^T(F^TF - I) \Delta X
$$

where $F^TF$ is Right Caychy Green Deformation. $F^TF - I$ is Green Lagrange Strain.

From deformation to potential energy:  
strain energy density: $\Psi$; $d\Omega$ is infinitesimal volume;

$$
\Psi(F(X))d\Omega
$$

$$
\int_{\Omega}\Psi(F(X))d\Omega
$$

$$
\mathbf x(X) = \underbrace{\begin{pmatrix}
\phi_0I & \phi_1I & \phi_2I & \phi_3I
 \end{pmatrix}}_{N(X)} \underbrace{\begin{pmatrix}
 \mathbf x_0 \\
 \mathbf x_1 \\
 \mathbf x_2 \\
 \mathbf x_3 \\
 \end{pmatrix}}_{q(t)}
$$

$$
\mathbf x(X) = \mathbf x_0 + \begin{pmatrix}
 \mathbf x_0 &
 \mathbf x_1 &
 \mathbf x_2 &
 \mathbf x_3 &
 \end{pmatrix} \begin{pmatrix}
 -1^T T^{-1} \\
 T^{-1}
 \end{pmatrix} (X - X_0)
$$

$$
F = \frac{\partial \mathbf x}{\partial X} = \begin{pmatrix}
 \mathbf x_0 &
 \mathbf x_1 &
 \mathbf x_2 &
 \mathbf x_3 &
 \end{pmatrix} \begin{pmatrix}
 -1^T T^{-1} \\
 T^{-1}
 \end{pmatrix}
$$

where $\begin{pmatrix}
 -1^T T^{-1} \\
 T^{-1}
 \end{pmatrix} = \frac{\partial \phi}{\partial X}$, $1^T = \begin{pmatrix} 1 & 1 & 1 & 1 \end{pmatrix} $

$$
V = vol \cdot \psi(F_0) \\
V_j = vol_j \cdot \psi(F_j(q_j)) \\
V = \sum_{j = 0}^{m-1} vol_j \cdot \psi(F_j(E_j q))
$$

```cpp

```

[`dV_linear_tetrahedron_dq.cpp`](./src/dV_linear_tetrahedron_dq.cpp)

Compute the gradient of the potential energy of a single tetrahedron. Note: you will need both `dpsi_neo_hookean_dq.cpp` and `quadrature_single_point.h` to do this.

```cpp
auto neohookean_linear_tet = [&](
    Eigen::Vector12d &dV,
    Eigen::Ref<const Eigen::VectorXd> q,
    Eigen::Ref<const Eigen::RowVectorXi> element)
{
    Eigen::Vector3d X = V.row(element(0)).transpose();

    Eigen::Matrix3d x;
    x << q.segment<3>(element(0) * 3),
        q.segment<3>(element(1) * 3),
        q.segment<3>(element(2) * 3),
        q.segment<3>(element(3) * 3);

    Eigen::Matrix43d dphi;
    dphi_linear_tetrahedron_dX(dphi, V, element, x);

    Eigen::Vector9d dpsi;
    dpsi_neo_hookean_dF(dpsi, x * dphi, C, D);

    Eigen::MatrixXd B(9, 12);
    B.setZero();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            B.block<3, 3>(i * 3, j * 3) = Eigen::Matrix3d::Identity() * dphi(j, i);
        }
    }

    dV = B.transpose() * dpsi;
};

quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);

```

[`d2V_linear_tetrahedron_dq2.cpp`](./src/d2V_linear_tetrahedron_dq2.cpp)

Compute the hessian of the potential energy of a single tetrahedron. Note: you will need both `d2psi_neo_hookean_dq2.cpp` and `quadrature_single_point.h` to do this.

```cpp

```

[`V_spring_particle_particle.cpp`](./src/V_spring_particle_particle.cpp)

The potential energy of a non-zero rest length spring attached to two vertices of the mesh. Use your code from the last assignment.

$$
V = \frac{1}{2}k(||q_1 - q_0|| - l_0)^2
$$

```cpp
V = 1. / 2. * stiffness * std::pow((q1 - q0).norm() - l0, 2);
```

[`dV_spring_particle_particle_dq.cpp`](./src/dV_spring_particle_particle_dq.cpp)

The gradient of the spring potential energy. Use your code from the last assignment.

```cpp

```

[`mass_matrix_linear_tetrahedron.cpp`](./src/mass_matrix_linear_tetrahedron.cpp)

Compute the dense mass matrix for a single tetrahedron.

$$
T = \frac{1}{2}\dot q^T(\underbrace{\int _{\Omega}\rho N(\mathbf X)^TN(\mathbf X)d\Omega}_{M_0})\dot q
$$

$$
M_0 = \int _{\Omega}\rho N(\mathbf X)^TN(\mathbf X)d\Omega = \int _{\Omega}\rho \begin{pmatrix}
\phi_0\phi_0 \mathbf I & \phi_0\phi_1 \mathbf I & \phi_0\phi_2 \mathbf I & \phi_0\phi_3 \mathbf I  \\
\phi_1\phi_0 \mathbf I & \phi_1\phi_1 \mathbf I & \phi_1\phi_2 \mathbf I & \phi_1\phi_3 \mathbf I  \\
\phi_2\phi_0 \mathbf I & \phi_2\phi_1 \mathbf I & \phi_2\phi_2 \mathbf I & \phi_2\phi_3 \mathbf I  \\
\phi_3\phi_0 \mathbf I & \phi_3\phi_1 \mathbf I & \phi_3\phi_2 \mathbf I & \phi_3\phi_3 \mathbf I  \\
 \end{pmatrix}d\Omega
$$

Evaluate each term separately:

$$
\rho \int _{\Omega}\phi_r(\mathbf X)\phi_s(\mathbf X) d\Omega \mathbf I
$$

integration using barycentric coordinates:

$$
6\rho \cdot vol \cdot \int_0^1 \int_0^{1 - \phi_1} \int_0^{1 - \phi_1 - \phi_2}(\phi_r\phi_s)d\phi_3d\phi_2d\phi_1
$$

$$
\phi_0(\mathbf X) = 1 - \phi_1(\mathbf X) - \phi_2(\mathbf X) - \phi_3(\mathbf X)
$$

For Diagonal Elements:

$$
6\rho \cdot vol \cdot \int_0^1 \int_0^{1 - \phi_1} \int_0^{1 - \phi_1 - \phi_2}(\phi_1\phi_1)d\phi_3d\phi_2d\phi_1 = 6\rho \cdot vol \cdot \frac{1}{60} = \frac{\rho \cdot vol}{10}
$$

The same answer for [0, 0], [2, 2], [3, 3].

For non-diagonale Elements:

$$
6\rho \cdot vol \cdot \int_0^1 \int_0^{1 - \phi_1} \int_0^{1 - \phi_1 - \phi_2}(\phi_1\phi_2)d\phi_3d\phi_2d\phi_1 = 6\rho \cdot vol \cdot \frac{1}{120} = \frac{\rho \cdot vol}{20}
$$

so $M_0$ can be written as:

$$
M_0 = \rho \cdot vol \cdot \begin{bmatrix}
1/10 & 1/20 & 1/20 & 1/20 \\
1/20 & 1/10 & 1/20 & 1/20 \\
1/20 & 1/20 & 1/10 & 1/20 \\
1/20 & 1/20 & 1/20 & 1/10 \\
\end{bmatrix}
$$

```cpp
M.setZero();

Eigen::Matrix4d M_0;

M_0 << 1.0 / 10.0, 1.0 / 20.0, 1.0 / 20.0, 1.0 / 20.0,
    1.0 / 20.0, 1.0 / 10.0, 1.0 / 20.0, 1.0 / 20.0,
    1.0 / 20.0, 1.0 / 20.0, 1.0 / 10.0, 1.0 / 20.0,
    1.0 / 20.0, 1.0 / 20.0, 1.0 / 20.0, 1.0 / 10.0;
M_0 *= density * volume;

for (int i = 0; i < 4; i++)
{
    for (int j = 0; j < 4; j++)
    {
        M.block<3, 3>(i * 3, j * 3) = Eigen::Matrix3d::Identity() * M_0(i, j);
    }
}
```

[`mass_matrix_mesh.cpp`](./src/mass_matrix_mesh.cpp)

Assemble the full mass matrix for the entire tetrahedral mesh.

```cpp
M.resize(qdot.size(), qdot.size());
M.setZero();
std::vector<Eigen::Triplet<double>> M_entries;

for (int i = 0; i < T.rows(); i++)
{
    Eigen::RowVector4i current_tetrahedron = T.row(i);

    Eigen::Matrix1212d current_tetrahedron_M;
    mass_matrix_linear_tetrahedron(current_tetrahedron_M, qdot, T.row(i), density, v0(i));

    // Iterate to populate 16 total d2V/d(corner_i)(corner_i) blocks
    for (int phi_i = 0; phi_i < 4; phi_i++)
    {
        for (int phi_j = 0; phi_j < 4; phi_j++)
        {
            // Fill up diagonal entry of each block
            for (int qdot_i = 0; qdot_i < 3; qdot_i++)
            {
                M_entries.push_back(
                    Eigen::Triplet<double>(
                        current_tetrahedron(phi_i) * 3 + qdot_i,
                        current_tetrahedron(phi_j) * 3 + qdot_i,
                        current_tetrahedron_M(phi_i * 3 + qdot_i, phi_j * 3 + qdot_i)));
            }
        }
    }
}

M.setFromTriplets(M_entries.begin(), M_entries.end());
```

[`assemble_forces.cpp`](./src/assemble_forces.cpp)

Assemble the global force vector for the finite element mesh.

```cpp
f.resize(q.size());
f.setZero();

for (int i = 0; i < T.rows(); i++)
{
    Eigen::RowVector4i current_tetrahedron = T.row(i);

    Eigen::Vector12d dV;

    dV_linear_tetrahedron_dq(dV, q, V, current_tetrahedron, v0(i), C, D);

    for (int vertex_i = 0; vertex_i < 4; vertex_i++)
    {
        f.segment<3>(current_tetrahedron(vertex_i) * 3) -= dV.segment<3>(vertex_i * 3);
    }
}
```

[`assemble_stiffness.cpp`](./src/assemble_stiffness.cpp)

Assemble the global stiffness matrix for the finite element mesh.

```cpp

```

[`build_skinning_matrix.cpp`](./src/build_skinning_matrix.cpp)

Build the skinning matrix that maps position from the coarse simulation mesh to the high resolution rendering mesh.

```cpp

```

[`fixed_point_constraints.cpp`](./src/fixed_point_constraints.cpp)

Use your code from the last assignment

```cpp
double min_vertex = V(0, 1);

for (unsigned int vi = 0; vi < V.rows(); vi++)
{
    min_vertex = (V(vi, 1) < min_vertex ? V(vi, 1) : min_vertex);
}

std::cout << min_vertex << std::endl;

for (unsigned int vi = 0; vi < V.rows(); vi++)
{
    if (std::abs(V(vi, 1) - min_vertex) <= tol)
    {
        indices.push_back(vi);
    }
}
```

[`pick_nearest_vertices.cpp`](./src/pick_nearest_vertices.cpp)

Use your code from the last assignment

```cpp

```

[`linearly_implicit_euler.h`](./include/linearly_implicit_euler.h)

Use your code from the last assignment

```cpp

```

[`newtons_method.h`](./include/newtons_method.h)

Implement Newton's method with backtracking line search. Use the following parameter values: alpha (initial step length) = 1, p (scaling factor) = 0.5, c (ensure sufficient decrease) = 1e-8.

```cpp

```

[`implicit_euler.h`](./include/implicit_euler.h)

Using your Newton's method, implement a fully implicit solver. To ensure reasonable performance, use a maximum of five (5) iterations.

```cpp


```
