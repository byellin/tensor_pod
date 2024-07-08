# # #This code should take a particular value of mu
# # # #and compute the solution to the Poisseuille flow
# #
# #
# # from fenics import *
# # import numpy as np
# # import matplotlib.pyplot as plt
# #
# # # Define domain and mesh
# # L = 2.0
# # mesh = RectangleMesh(Point(-1, -1), Point(1, 1), 8, 8)
# #
# # # Define function spaces
# # V = VectorFunctionSpace(mesh, 'P', 2)
# # Q = FunctionSpace(mesh, 'P', 1)
# #
# # # Define boundary conditions
# # inflow = Expression(('1.0*(1 - x[1]*x[1])', '0'), degree=2)
# # noslip = Constant((0, 0))
# #
# # bc_inflow = DirichletBC(V, inflow, 'near(x[0], -1)')
# # bc_walls = DirichletBC(V, noslip, 'near(x[1], -1) || near(x[1], 1)')
# # bc_outflow = DirichletBC(Q, Constant(0), 'near(x[0], 1)')
# #
# # bcu = [bc_inflow, bc_walls]
# # bcp = [bc_outflow]
# #
# # # Define trial and test functions
# # u = TrialFunction(V)
# # v = TestFunction(V)
# # p = TrialFunction(Q)
# # q = TestFunction(Q)
# #
# # # Define functions
# # u_n = Function(V)
# # u_ = Function(V)
# # p_n = Function(Q)
# # p_ = Function(Q)
# #
# # # Define coefficients and parameter range
# # rho = 1.0
# # mu_vals = [10**k for k in range(-10,11)]
# # dt = 0.01
# # T = 1
# #
# # # Initialize storage for velocities and pressures
# # num_steps = int(T / dt)
# # velocities_tensor = np.zeros((V.dim(),num_steps, len(mu_vals)))
# # pressures_tensor = np.zeros((Q.dim(),num_steps, len(mu_vals)))
# # print("velocities shape", velocities_tensor.shape)
# # print("pressures shape", pressures_tensor.shape)
# # # Loop over viscosity values
# # for mu_count, mu in enumerate(mu_vals):
# #
# #     # Define expressions used in variational forms
# #     U = 0.5 * (u_n + u)
# #     n = FacetNormal(mesh)
# #     f = Constant((0, 0))
# #     k = Constant(dt)
# #     mu = Constant(mu)
# #     rho = Constant(rho)
# #
# #     # Define the stress tensor
# #     def epsilon(u):
# #         return sym(nabla_grad(u))
# #
# #     def sigma(u, p):
# #         return 2*mu*epsilon(u) - p*Identity(len(u))
# #
# #     # Define variational problem for step 1
# #     F1 = rho * dot((u - u_n) / k, v) * dx + rho * dot(dot(u_n, nabla_grad(u_n)), v) * dx \
# #          + inner(sigma(U, p_n), epsilon(v)) * dx + dot(p_n * n, v) * ds - dot(mu * nabla_grad(U) * n, v) * ds - dot(f,
# #                                                                                                                     v) * dx
# #     a1 = lhs(F1)
# #     L1 = rhs(F1)
# #
# #     # Define variational problem for step 2
# #     a2 = dot(nabla_grad(p), nabla_grad(q)) * dx
# #     L2 = dot(nabla_grad(p_n), nabla_grad(q)) * dx - (1 / k) * div(u_) * q * dx
# #
# #     # Define variational problem for step 3
# #     a3 = dot(u, v) * dx
# #     L3 = dot(u_, v) * dx - k * dot(nabla_grad(p_ - p_n), v) * dx
# #
# #     # Assemble matrices
# #     A1 = assemble(a1)
# #     A2 = assemble(a2)
# #     A3 = assemble(a3)
# #
# #     # Time-stepping
# #     t = 0
# #     step = 0
# #     while t < T:
# #         t += dt
# #
# #         # Step 1: Tentative velocity step
# #         b1 = assemble(L1)
# #         [bc.apply(A1, b1) for bc in bcu]
# #         solve(A1, u_.vector(), b1)
# #
# #         # Step 2: Pressure correction step
# #         b2 = assemble(L2)
# #         [bc.apply(A2, b2) for bc in bcp]
# #         solve(A2, p_.vector(), b2)
# #
# #         # Step 3: Velocity correction step
# #         b3 = assemble(L3)
# #         solve(A3, u_.vector(), b3)
# #
# #         # Update previous solution
# #         u_n.assign(u_)
# #         p_n.assign(p_)
# #
# #         # Store the results
# #         velocities_tensor[:,step, mu_count] = u_.vector().get_local()
# #         pressures_tensor[:, step, mu_count] = p_.vector().get_local()
# #
# #         step += 1
# #
# #         # Print progress
# #         print(f'Time step {t:.2f}/{T:.2f} completed')
# #
# #     # Save solution in VTK format
# #     file_u = File(f'velocity_mu{float(mu):.3f}.pvd')
# #     file_u << u_
# #
# #     file_p = File(f'pressure_mu{float(mu):.3f}.pvd')
# #     file_p << p_
# #
# #     # Plot solution (optional)
# #     # plt.figure()
# #     # # plt.imshow(np.asarray(u_))
# #     # plot(u_)
# #     # plt.title(f'Velocity (mu={float(mu):.3f})')
# #     # plt.xlabel('x')
# #     # plt.ylabel('y')
# #     # # plt.colorbar()
# #     # plt.savefig(f'velocity_mu{float(mu):.3f}.png')
# #     #
# #     # plt.figure()
# #     # # plt.imshow(np.asarray(p_))
# #     # plot(p_)
# #     # plt.title(f'Pressure (mu={float(mu):.3f})')
# #     # plt.xlabel('x')
# #     # plt.ylabel('y')
# #     # # plt.colorbar()
# #     # plt.savefig(f'pressure_mu{float(mu):.3f}.png')
# #
# # # Save tensors to file
# # np.save('velocities_tensor.npy', velocities_tensor)
# # np.save('pressures_tensor.npy', pressures_tensor)
# #
# # # print("velocities: ", velocities_tensor)
# # # print("pressures: ", pressures_tensor)


from fenics import *
import numpy as np
# from mshr import *

# Create a square domain mesh
nx, ny = 64, 64
mesh = UnitSquareMesh(nx, ny)

# Define function spaces (P2 for velocity, P1 for pressure)
V = VectorFunctionSpace(mesh, 'P', 2)
Q = FunctionSpace(mesh, 'P', 1)

# Define boundary conditions
inflow_profile = Expression(('4.0*x[1]*(1.0 - x[1])', '0'), degree=2)

bc_inflow = DirichletBC(V, inflow_profile, 'near(x[0], 0)')
bc_outflow = DirichletBC(Q, Constant(0), 'near(x[0], 1)')
bc_walls = [DirichletBC(V, Constant((0, 0)), 'near(x[1], 0)'),
            DirichletBC(V, Constant((0, 0)), 'near(x[1], 1)')]

bcs = [bc_inflow, bc_walls[0], bc_walls[1]]

# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Define functions for solutions at previous and current time steps
u_n = Function(V)
u_ = Function(V)
p_n = Function(Q)
p_ = Function(Q)

# Define constants
rho = 1.0  # density
# mu = 1.0  # dynamic viscosity
mu_vals = [10**-k for k in range(-10,2)]

T=2
dt=0.1
num_steps = int(T / dt)
velocities_tensor = np.zeros((V.dim(),num_steps, len(mu_vals)))
pressures_tensor = np.zeros((Q.dim(),num_steps, len(mu_vals)))




for mu_count,mu in enumerate(mu_vals):

    print("mu is ", mu)
    # Define source term
    f = Constant((0, 0))

    # Define variational problem for step 1
    F1 = (rho * dot((u - u_n) / dt, v) * dx
          + rho * dot(dot(u_n, nabla_grad(u_n)), v) * dx
          + mu * inner(nabla_grad(u), nabla_grad(v)) * dx
          - p_n * div(v) * dx
          - dot(f, v) * dx)
    a1, L1 = lhs(F1), rhs(F1)

    # Define variational problem for step 2
    a2 = dot(nabla_grad(p), nabla_grad(q)) * dx
    L2 = dot(nabla_grad(p_n), nabla_grad(q)) * dx - (1 / dt) * div(u_) * q * dx

    # Define variational problem for step 3
    a3 = dot(u, v) * dx
    L3 = dot(u_, v) * dx - dt * dot(nabla_grad(p_ - p_n), v) * dx

    # Assemble matrices
    A1 = assemble(a1)
    A2 = assemble(a2)
    A3 = assemble(a3)

    print("right before t loop")

    # Time-stepping
    t = 0.0
    step = 0

    while t < T:
        t += dt
        print('t', t)
        # Step 1: Tentative velocity step
        b1 = assemble(L1)
        [bc.apply(A1, b1) for bc in bcs]
        solve(A1, u_.vector(), b1)

        # Step 2: Pressure correction step
        b2 = assemble(L2)
        [bc.apply(A2, b2) for bc in [bc_outflow]]
        solve(A2, p_.vector(), b2)

        # Step 3: Velocity correction step
        b3 = assemble(L3)
        solve(A3, u_.vector(), b3)

        # Update previous solution
        u_n.assign(u_)
        p_n.assign(p_)

        # Store the results
        velocities_tensor[:,step, mu_count] = u_.vector().get_local()
        pressures_tensor[:, step, mu_count] = p_.vector().get_local()

        step += 1

    velocities_has_nan = np.isnan(velocities_tensor).any()
    pressures_has_nan = np.isnan(pressures_tensor).any()

    print("Velocities has nan values: ", velocities_has_nan)
    print("Pressures has nan values: ", pressures_has_nan)

    # Save solution in VTK format
    File("velocity.pvd") << u_
    File("pressure.pvd") << p_
