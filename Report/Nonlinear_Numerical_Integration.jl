


f(x) = [x[3] x[4] -(L1.*l2.*m2.*(0.5*J2.*x[3].^2*sin(2.0*x[2]) - 0.5*J2.*x[3].^2*sin(2.0*x[2]) + b2.*x[4] + g.*l2.*m2.*sin(x[2]) - 0.5*l2.^2*m2.*x[3].^2*sin(2.0*x[2]) - 0).*cos(x[2]) + (J2 + l2.^2*m2).*(J2.*x[3].*x[4].*sin(2.0*x[2]) - J2.*x[3].*x[4].*sin(2.0*x[2]) + L1.*l2.*m2.*x[4].^2*sin(x[2]) - b1.*x[3] - l2.^2*m2.*x[3].*x[4].*sin(2.0*x[2]) + (u0 + R(t) - dot(K,x))))./(L1.^2*l2.^2*m2.^2*cos(x[2]).^2 - (J2 + l2.^2*m2).*(J1 - J2.*sin(x[2]).^2 + J2 + J2.*sin(x[2]).^2 + L1.^2*m2 + l1.^2*m1 + l2.^2*m2.*sin(x[2]).^2)) (L1.*l2.*m2.*(J2.*x[3].*x[4].*sin(2.0*x[2]) - J2.*x[3].*x[4].*sin(2.0*x[2]) + L1.*l2.*m2.*x[4].^2*sin(x[2]) - b1.*x[3] - l2.^2*m2.*x[3].*x[4].*sin(2.0*x[2]) + (u0 + R(t) - dot(K,x))).*cos(x[2]) + (0.5*J2.*x[3].^2*sin(2.0*x[2]) - 0.5*J2.*x[3].^2*sin(2.0*x[2]) + b2.*x[4] + g.*l2.*m2.*sin(x[2]) - 0.5*l2.^2*m2.*x[3].^2*sin(2.0*x[2]) - 0).*(J1 - J2.*sin(x[2]).^2 + J2 + J2.*sin(x[2]).^2 + L1.^2*m2 + l1.^2*m1 + (l2.^2).*m2.*sin(x[2]).^2))./(L1.^2*l2.^2*m2.^2*cos(x[2]).^2 + (-J2 - l2.^2*m2).*(J1 - J2.*sin(x[2]).^2 + J2 + J2.*sin(x[2]).^2 + L1.^2*m2 + l1.^2*m1 + l2.^2*m2.*sin(x[2]).^2))]


f(z) = [z[3] z[4] (1/(J0hat*J2hat+J2hat^2*sin(z[2])*sin(z[2])-m2^2*L1^2*l2^2*cos(z[2])*cos(z[2])))*(-J2hat*b1*z[3]+m2*L1*l2*cos(z[2])*b2*z[4]-J2hat^2*sin(2*z[2])*z[3]*z[4]-0.5*J2hat*m2*L1*l2*cos(z[2])*sin(2*z[2])*z[3]^2+J2hat*m2*L1*l2*sin(z[2])*z[4]^2+J2hat*(u0 + R(t) - dot(K,([z[1] z[2] z[3] z[4]] - [0. π 0. 0.])))+0.5*m2^2*l2^2*L1*sin(2*z[2])*g) (m2*L1*l2*cos(z[2])*b1*z[3]-b2*(J0hat+J2hat*sin(z[2])^2)*z[4]+m2*
    L1*l2*J2hat*cos(z[2])*sin(2*z[2])*z[3]*z[4]-0.5*sin(2*z[2])*(J0hat*J2hat+J2hat^2*sin(z[2])^2)*z[3]^2-0.5*m2^2*L1^2*l2^2*sin(2*z[2])*z[4]^2-m2*L1*l2*cos(z[2])*(u0 + R(t) - dot(K,([z[1] z[2] z[3] z[4]] - [0. π 0. 0.])))-m2*l2*sin(z[2])*(J0hat+J2hat*sin(z[2])^2)*g)/(J0hat*J2hat+J2hat^2*sin(z[2])^2-m2^2*L1^2*l2^2*cos(z[2])^2)]






x0 = z0 - [0. π 0. 0.]
tf = 1.
h =.01
f(x0)

include("../GL_RK.jl")
x_GL,t = GL_RK.gl_rk(f,tf,h,x0)
plot(t,x_GL[:,1],label = "theta")
    plot!(t,x_GL[:,2],label = "thetadot")
    xlabel!("time")
    ylabel!("Angle (Radians)")
    title!("Gauss Legendre Runge Kutta")



include("../../COMPUTATIONAL-DYNAMICS/Homework/HW05-ODEs_Part2/BackwardEuler.jl")

x_BE,t = BackwardEuler.beuler(f,tf,h,z0)
    plot(t,x_BE[:,1],label = "qdot")
        plot!(t,x_BE[:,2],label = "q")
        xlabel!("Time")
        ylabel!("State")
        title!("Backwards Euler")














###
function odeNL(dx,x,p,t)
        J1,J2,L1,l2,m2,b1,b2,g, R, K = p


        u = u0 + R(t) - dot(K,x)

        dx[1] = x[3]
        dx[2] = x[4]
        dx[3] = -(L1.*l2.*m2.*(0.5*J2.*x[3].^2*sin(2.0*x[2]) - 0.5*J2.*x[3].^2*sin(2.0*x[2]) + b2.*x[4] + g.*l2.*m2.*sin(x[2]) - 0.5*l2.^2*m2.*x[3].^2*sin(2.0*x[2]) - 0).*cos(x[2]) + (J2 + l2.^2*m2).*(J2.*x[3].*x[4].*sin(2.0*x[2]) - J2.*x[3].*x[4].*sin(2.0*x[2]) + L1.*l2.*m2.*x[4].^2*sin(x[2]) - b1.*x[3] - l2.^2*m2.*x[3].*x[4].*sin(2.0*x[2]) + u))./(L1.^2*l2.^2*m2.^2*cos(x[2]).^2 - (J2 + l2.^2*m2).*(J1 - J2.*sin(x[2]).^2 + J2 + J2.*sin(x[2]).^2 + L1.^2*m2 + l1.^2*m1 + l2.^2*m2.*sin(x[2]).^2))

        dx[4] = (L1.*l2.*m2.*(J2.*x[3].*x[4].*sin(2.0*x[2]) - J2.*x[3].*x[4].*sin(2.0*x[2]) + L1.*l2.*m2.*x[4].^2*sin(x[2]) - b1.*x[3] - l2.^2*m2.*x[3].*x[4].*sin(2.0*x[2]) + u).*cos(x[2]) + (0.5*J2.*x[3].^2*sin(2.0*x[2]) - 0.5*J2.*x[3].^2*sin(2.0*x[2]) + b2.*x[4] + g.*l2.*m2.*sin(x[2]) - 0.5*l2.^2*m2.*x[3].^2*sin(2.0*x[2]) - 0).*(J1 - J2.*sin(x[2]).^2 + J2 + J2.*sin(x[2]).^2 + L1.^2*m2 + l1.^2*m1 + (l2.^2).*m2.*sin(x[2]).^2))./(L1.^2*l2.^2*m2.^2*cos(x[2]).^2 + (-J2 - l2.^2*m2).*(J1 - J2.*sin(x[2]).^2 + J2 + J2.*sin(x[2]).^2 + L1.^2*m2 + l1.^2*m1 + l2.^2*m2.*sin(x[2]).^2))


end

tspan = (0.0,1)
θ0 = 10 # the angle offset in DEGREES from vertical
z0 = [0 (π - θ0*π/180) 0 0.]
x0 = z0 - [0 π 0 0.]
p = [J1,J2,L1,l2,m2,b1,b2,g, R, K]
prob = ODEProblem(odeNL,z0,tspan,p)
sol = DifferentialEquations.solve(prob,BS3())
plot(sol, vars = (0,2))
plot(sol, vars = (0,3))
