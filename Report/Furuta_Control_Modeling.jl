using LinearAlgebra
using DifferentialEquations
using Plots
import Polynomials


J0hat = .0025
J2hat = .0001
L1 = .5
l2 = .3
m2 = .2
b1 = .0001
b2 = .0001
g = 9.8

A31 = 0
A32 = g*m2^2*l2^2*L1/(J0hat*J2hat-m2^2*L1^2*l2^2)
A33 = -b1*J2hat/(J0hat*J2hat-m2^2*L1^2*l2^2)
A34 = -b2*m2*l2*L1/(J0hat*J2hat-m2^2*L1^2*l2^2)
A41 = 0
A42 = g*m2*l2*J0hat/(J0hat*J2hat-m2^2*L1^2*l2^2)
A43 = -b1*m2*l2*L1/(J0hat*J2hat-m2^2*L1^2*l2^2)
A44 = -b2*J0hat/(J0hat*J2hat-m2^2*L1^2*l2^2)

B31 = J2hat/(J0hat*J2hat-m2^2*L1^2*l2^2)
B32 = m2*L1*l2/(J0hat*J2hat-m2^2*L1^2*l2^2)
B41 = m2*L1*l2/(J0hat*J2hat-m2^2*L1^2*l2^2)
B42 = J0hat/(J0hat*J2hat-m2^2*L1^2*l2^2)

A = [0 0 1 0;0 0 0 1;A31 A32 A33 A34;A41 A42 A43 A44]
# assuming that tau2 = 0
B = [0,0,B31,B41]


OS = 5/100
Ts = .3
zeta = sqrt(log(OS)^2/(pi^2+log(OS)^2))
omega = 4/Ts/zeta
Ts2 = 5*Ts
omega2 = 4/Ts2/zeta

# p1 = Polynomials.poly([-1,-2,-3,-4])
s1 = -zeta*omega+1im*omega*sqrt(1 - zeta^2)
s2 = -zeta*omega2+1im*omega2*sqrt(1 - zeta^2)

p1 = Polynomials.poly([s1, conj(s1), s2, conj(s2)  ])
Λ = zeros(4,4)
for i = 0:4
        global Λ += p1.a[i+1]*A^i
end
Λ = real(Λ)

e_c = hcat(B, A*B, A^2*B, A^3*B)
K = (e_c\Λ)[end,:]
eigvals(A - B*K')


# K = (e_c\lambda)[:,end]
u0 = 0
R(t) = 0
poles = eigvals(A-B*K')


#this funtion takes the expressions for both accelerations and takes τ1 to be u and takes τ2 to be 0.
function ode2(dz,z,p,t)
        J0hat,J2hat,L1,l2,m2,b1,b2,g, R, K = p


        x = z - [0. π 0. 0.]
        u = u0 + R(t) - dot(K,x)

        dz[1] = z[3]
        dz[2] = z[4]
        dz[3] = (1/(J0hat*J2hat+J2hat^2*sin(z[2])*sin(z[2])
        -m2^2*L1^2*l2^2*cos(z[2])*cos(z[2])))*(-J2hat*b1*z[3]
        +m2*L1*l2*cos(z[2])*b2*z[4]-J2hat^2*sin(2*z[2])*z[3]*z[4]
        -0.5*J2hat*m2*L1*l2*cos(z[2])*sin(2*z[2])*z[3]^2
        +J2hat*m2*L1*l2*sin(z[2])*z[4]^2+J2hat*u+0.5*m2^2*l2^2*L1*sin(2*z[2])*g)

        dz[4] = (m2*L1*l2*cos(z[2])*b1*z[3]-b2*(J0hat+J2hat*sin(z[2])^2)*z[4]+m2*
        L1*l2*J2hat*cos(z[2])*sin(2*z[2])*z[3]*z[4]
        -0.5*sin(2*z[2])*(J0hat*J2hat+J2hat^2*sin(z[2])^2)*z[3]^2
        -0.5*m2^2*L1^2*l2^2*sin(2*z[2])*z[4]^2-m2*L1*l2*cos(z[2])*u
        -m2*l2*sin(z[2])*(J0hat+J2hat*sin(z[2])^2)*g)/
        (J0hat*J2hat+J2hat^2*sin(z[2])^2-m2^2*L1^2*l2^2*cos(z[2])^2)


end

tspan = (0.0,5)
θ0 = 10 # the angle offset in DEGREES from vertical
z0 = [0 (π - θ0*π/180) 0 0.]
p = [J0hat,J2hat,L1,l2,m2,b1,b2,g, R, K]
prob = ODEProblem(ode2,z0,tspan,p)
sol = DifferentialEquations.solve(prob)
plot(sol, vars = (0,2))
plot(sol, vars = (0,3))

### NOTE this one does both arms
let sol = sol
        @gif for i = 1:size(sol[1,:],1)
                x = 1
                y = .5
                z = 1
                pyplot()
                # r1 is the actuated arm
                r1(θ1) = [cos(θ1),sin(θ1), 0. ]
                # where θ1 is sol[1,:] and sol has 27 elements
                r2(θ1,θ2) = [r1(θ1)[1],r1(θ1)[2]+sin(θ2),
                r1(θ1)[3]-cos(θ2)]
                path3d([0., r1(sol[1,i])[1]], [0., r1(sol[1,i])[2]],
                [0., r1(sol[1,i])[3]],xlims=(-x,x),
                ylims=(-y,y),zlims=(-0,z))
                path3d!([r1(sol[1,i])[1], r1(sol[1,i])[1]],
                 [r1(sol[1,i])[2], r2(sol[1,i],sol[2,i])[2]],
                 [r1(sol[1,i])[3],r2(sol[1,i],sol[2,i])[3]])
        end every 1
end
