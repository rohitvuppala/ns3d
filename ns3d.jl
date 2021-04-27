#Import useful stuff
using Plots
using OffsetArrays
using WriteVTK

module consts
    γ = 1.4
    Ma= 0.08
end

#Create the grid
function grid_init(nx,ny,nz)
    #Intialise the grid
    lx = 2*pi
    ly = 2*pi
    lz = 2*pi

    dx = lx/nx
    dy = ly/ny
    dz = lz/nz

    x = zeros(nx+7,ny+7,nz+7)
    y = zeros(nx+7,ny+7,nz+7)
    z = zeros(nx+7,ny+7,nz+7)

    #Offset arrays to reflect ghost points
    x = OffsetArray(x,-3:nx+3,-3:ny+3,-3:nz+3)
    y = OffsetArray(y,-3:nx+3,-3:ny+3,-3:nz+3)
    z = OffsetArray(z,-3:nx+3,-3:ny+3,-3:nz+3)

    #Set the values for x,y,z and create grid
    for k in -3:nz+3
        for j in -3:ny+3
            for i in -3:nx+3
                x[i,j,k] = 0.0 + dx*i
                y[i,j,k] = 0.0 + dy*j
                z[i,j,k] = 0.0 + dz*k
            end
        end
    end

    return x,y,z,dx,dy,dz
end

#Initialise
function init_3d(Ma,x,y,z,nx,ny,nz)

    γ = consts.γ

    # To store cons ρ,ρu,ρv,ρw,ρe
    q = zeros(5,nx+7,ny+7,nz+7)
    q = OffsetArray(q,1:5,-3:nx+3,-3:ny+3,-3:nz+3)

    # To store ρ,u,v,w,p,e
    pq= zeros(6,nx+7,ny+7,nz+7)
    pq= OffsetArray(pq,1:6,-3:nx+3,-3:ny+3,-3:nz+3)

    #Temp arrays
    cons = zeros(5)
    xyz  = zeros(3)

    #Initialise
    for k in 0:nz
        for j in 0:ny
            for i in 0:nx

                xyz = [x[i,j,k],y[i,j,k],z[i,j,k]]

                #Primitives
                ρ = 1.0
                u = sin(xyz[1])*cos(xyz[2])*cos(xyz[3])
                v =-cos(xyz[1])*sin(xyz[2])*cos(xyz[3])
                w = 0.0
                p = 1/(γ*Ma^2)
                  + (cos(2*xyz[1])+cos(2*xyz[2]))*(cos(2*xyz[3])+2.0)/16.0
                e = p/(ρ*(γ-1))
                  + 0.5*(u^2 + v^2 + w^2)

                pq[1:6,i,j,k] = [ρ,u,v,w,p,e]

                #Conservatives
                cons[1] = ρ
                cons[2] = ρ*u
                cons[3] = ρ*v
                cons[4] = ρ*w
                cons[5] = ρ*e

                q[1:5,i,j,k] = cons[1:5]
            end
        end
    end

    return q,pq
end

#Calculate time step
function calc_dt(cfl,γ,q,nx,ny,nz,dx,dy,dz)
    a = 0.0
    a = maximum([a,0.0])
    for k in 0:nz
        for j in 0:ny
            for i in 0:nx
                ρ,ρu,ρv,ρw,ρe = q[:,i,j,k]
                u,v,w,e       = ρu/ρ, ρv/ρ, ρw/ρ, ρe/ρ
                p = ρ*(γ-1)*(e-0.5*(u^2+v^2+w^2))
                c = sqrt(γ*p/ρ)
                a = maximum([a,abs(u),abs(u+c),abs(u-c)
                             ,abs(v),abs(v+c),abs(v-c)
                             ,abs(w),abs(w+c),abs(w-c)])

            end
        end
    end

    dt = cfl* minimum([dx,dy,dz])/a

    return dt
end

#Boundary Conditions
function expbc!(q,nx,ny,nz)
    #Periodic Boundary Conditions
    # In x-direction
    q[:,-1,:,:]   = q[:,nx,:,:]
    q[:,-2,:,:]   = q[:,nx-1,:,:]
    q[:,-3,:,:]   = q[:,nx-2,:,:]
    q[:,nx+1,:,:] = q[:,0,:,:]
    q[:,nx+2,:,:] = q[:,1,:,:]
    q[:,nx+3,:,:] = q[:,2,:,:]

    # In y-direction
    q[:,:,-1,:]   = q[:,ny,:,:]
    q[:,:,-2,:]   = q[:,ny-1,:,:]
    q[:,:,-3,:]   = q[:,ny-2,:,:]
    q[:,:,ny+1,:] = q[:,0,:,:]
    q[:,:,ny+2,:] = q[:,1,:,:]
    q[:,:,ny+3,:] = q[:,2,:,:]

    # In z-direction
    q[:,:,:,-1]   = q[:,:,:,nz]
    q[:,:,:,-2]   = q[:,:,:,nz-1]
    q[:,:,:,-3]   = q[:,:,:,nz-2]
    q[:,:,:,nz+1] = q[:,:,:,0]
    q[:,:,:,nz+2] = q[:,:,:,1]
    q[:,:,:,nz+3] = q[:,:,:,2]

end

#Create function for writing output
function output_data(q,x,y,z,nx,ny,nz)
    γ = consts.γ
    ρ = q[1,0:nx,0:ny,0:nz]
    u = q[2,0:nx,0:ny,0:nz]./ρ
    v = q[3,0:nx,0:ny,0:nz]./ρ
    w = q[4,0:nx,0:ny,0:nz]./ρ
    e = q[5,0:nx,0:ny,0:nz]./ρ
    p = ρ.*(γ-1).*(e-0.5.*(u.^2+v.^2+w.^2))
    c = sqrt.(γ.*p./ρ)

    vtkfile = vtk_grid("output",x[0:nx,0:ny,0:nz],y[0:nx,0:ny,0:nz],z[0:nx,0:ny,0:nz])
    vtkfile["Velocity"] = (u,v,w)
    vtkfile["Density" ] = ρ
    vtkfile["Pressure"] = p
    vtkfile["Speed of Sound"] = c

    vtk_save(vtkfile)
    return nothing
end

#Create a function for the ns3d run
function ns3d(cfl=0.5,nx=32,ny=32,nz=32,nitermax=10000,tend=1.0)

    #Setup required
    γ  = consts.γ
    Ma = 0.08
    time = 0.0

    #Create the grid
    x,y,z,dx,dy,dz = grid_init(nx,ny,nz)

    #Initialise
    q,pq_init = init_3d(Ma,x,y,z,nx,ny,nz)
    qnew = zeros(5,nx+7,ny+7,nz+7)
    qnew = OffsetArray(q,1:5,-3:nx+3,-3:ny+3,-3:nz+3)

    #Calc_dt
    dt = calc_dt(cfl,γ,q,nx,ny,nz,dx,dy,dz)

    #Boundary Conditions
    expbc!(q,nx,ny,nz)

    for niter in 1:nitermax
        dt = calc_dt(cfl,γ,q,nx,ny,nz,dx,dy,dz)
        if time+dt > tend
            dt = tend-time
        end

        expbc!(q,nx,ny,nz)

        qnew = tvdrk3(nx,ny,nz,dx,dy,dz,q,dt)

        expbc!(qnew,nx,ny,nz)

        q = copy(qnew)
        time = time + dt

        @info niter,time,dt
        if (time >= tend)
            break
        end
    end


    #Output data
    output_data(q,x,y,z,nx,ny,nz)

    return nothing
end

#Flux calculation
function flux(nx,ny,nz,q)
    F = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-2:nx+2,-2:ny+2,-2:nz+2)

    γ = 1.4
    ρ = q[1,:]
    u = q[2,:]./ρ
    v = q[3,:]./ρ
    w = q[4,:]./ρ
    e = q[5,:]./ρ
    p = (γ-1)*(q[5,:] - 0.5*ρ.*(u.^2+v.^2+w.^2))

    F[1,:] = ρ.*u
    F[2,:] = ρ.* (u.^2) + p
    F[3,:] = u.* (ρ.*e + p)

    return F
end

#RHS calculation
function rhs(nx,ny,nz,dx,dy,dz,q)

    r = OffsetArray(zeros(3,nx+5),1:3,-2:nx+2)

    qL,qR = weno5(nx,q)
    FL = flux(nx,qL)
    FR = flux(nx,qR)

    cs = cs_weno(q,nx)

    #Compute flux with Rusanov
    F = OffsetArray(zeros(3,nx+5),1:3,-2:nx+2)

    for i in 0:nx-1
        F[:,i] = 0.5*(FL[:,i]+FR[:,i]) + 0.5*cs[i]*(qL[:,i]-qR[:,i])
    end

    for i in 1:nx-1
        r[:,i] = -(F[:,i]-F[:,i-1])./dx
    end

    return r
end

#Propagation speed calculation
function cs_weno(q,nx,ny,nz)
    γ = 1.4
    ρ = q[1,:]
    u = q[2,:]./ρ
    e = q[3,:]./ρ
    p = (γ-1)*(q[3,:] - 0.5*(ρ.*u.^2))

    a  = sqrt.(γ*p./ρ)

    cs = OffsetArray(zeros(nx+1),0:nx)
    r  = OffsetArray(zeros(nx+5),-2:nx+2)

    for i in -2:nx+2
        r[i] = maximum([abs(u[i]),abs(u[i]-a[i]),abs(u[i]+a[i])])
    end
    for i in 0:nx-1
        cs[i] = maximum([abs(r[i-2]),abs(r[i-1]),abs(r[i]),
                         abs(r[i+1]),abs(r[i+2]),abs(r[i+3])])
    end
    return cs
end

#Time stepping RK3
function tvdrk3(nx,ny,nz,dx,dy,dz,q,dt)
    qq = copy(q)
    qn = copy(q)

    #First step
    r = rhs(nx,dx,q)
    for i in 1:nx-1
        qq[:,i] = q[:,i] + dt*r[:,i]
    end

    #Second step
    r = rhs(nx,dx,qq)
    for i in 1:nx-1
        qq[:,i] = 0.75*q[:,i] + 0.25*qq[:,i] + 0.25*dt*r[:,i]
    end

    #Third Step
    r = rhs(nx,dx,qq)
    for i in 1:nx-1
        qn[:,i] = 1/3*q[:,i] + 2/3*qq[:,i] + 2/3*dt*r[:,i]
    end
    
    return qn
end
