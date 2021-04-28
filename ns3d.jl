#Import useful stuff
using Plots
using OffsetArrays
using WriteVTK

module consts
    γ = 1.4
    Ma= 0.08
end

#Create the grid
function grid_init(nx,ny,nz,idbg)
    #Intialise the grid
    lx = 2*pi
    ly = 2*pi
    lz = 2*pi

    if (idbg==1)
        lx = 1.0
        ly = 1.0
        lz = 1.0
    end
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
function init_3d(Ma,x,y,z,nx,ny,nz,idbg)

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

                if (idbg==1)
                    #Sod Shock Testing
                    if (x[i,j,k] <= 0.5)
                        ρ = 1.0
                        u = 0.0
                        v = 0.0
                        w = 0.0
                        p = 1.0
                        e = p/(ρ*(γ-1)) +
                            0.5*(u^2 + v^2 + w^2)
                    else
                        ρ = 0.125
                        u = 0.0
                        v = 0.0
                        w = 0.0
                        p = 0.1
                        e = p/(ρ*(γ-1)) +
                            0.5*(u^2 + v^2 + w^2)
                    end
                    pq[1:6,i,j,k] = [ρ,u,v,w,p,e]

                else

                    #Primitives
                    ρ = 1.0
                    u = sin(xyz[1])*cos(xyz[2])*cos(xyz[3])
                    v =-cos(xyz[1])*sin(xyz[2])*cos(xyz[3])
                    w = 0.0
                    p = 1/(γ*Ma^2) +
                        (cos(2*xyz[1])+cos(2*xyz[2]))*(cos(2*xyz[3])+2.0)/16.0
                    e = p/(ρ*(γ-1)) +
                        0.5*(u^2 + v^2 + w^2)

                    pq[1:6,i,j,k] = [ρ,u,v,w,p,e]

                end
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
function expbc!(q,nx,ny,nz,idbg)
    if (idbg==1)
        # In x-direction
        q[:,-1,:,:]   = q[:,0,:,:]
        q[:,-2,:,:]   = q[:,0,:,:]
        q[:,-3,:,:]   = q[:,0,:,:]
        q[:,nx+1,:,:] = q[:,nx,:,:]
        q[:,nx+2,:,:] = q[:,nx,:,:]
        q[:,nx+3,:,:] = q[:,nx,:,:]


        # In y-direction
        q[:,:,-1,:]   = q[:,:,0,:]
        q[:,:,-2,:]   = q[:,:,0,:]
        q[:,:,-3,:]   = q[:,:,0,:]
        q[:,:,ny+1,:] = q[:,:,ny,:]
        q[:,:,ny+2,:] = q[:,:,ny,:]
        q[:,:,ny+3,:] = q[:,:,ny,:]

        # In z-direction
        q[:,:,:,-1]   = q[:,:,:,0]
        q[:,:,:,-2]   = q[:,:,:,0]
        q[:,:,:,-3]   = q[:,:,:,0]
        q[:,:,:,nz+1] = q[:,:,:,nz]
        q[:,:,:,nz+2] = q[:,:,:,nz]
        q[:,:,:,nz+3] = q[:,:,:,nz]

    else
        #Periodic Boundary Conditions
        # In x-direction
        q[:,-1,:,:]   = q[:,nx,:,:]
        q[:,-2,:,:]   = q[:,nx-1,:,:]
        q[:,-3,:,:]   = q[:,nx-2,:,:]
        q[:,nx+1,:,:] = q[:,0,:,:]
        q[:,nx+2,:,:] = q[:,1,:,:]
        q[:,nx+3,:,:] = q[:,2,:,:]

        # In y-direction
        q[:,:,-1,:]   = q[:,:,ny,:]
        q[:,:,-2,:]   = q[:,:,ny-1,:]
        q[:,:,-3,:]   = q[:,:,ny-2,:]
        q[:,:,ny+1,:] = q[:,:,0,:]
        q[:,:,ny+2,:] = q[:,:,1,:]
        q[:,:,ny+3,:] = q[:,:,2,:]

        # In z-direction
        q[:,:,:,-1]   = q[:,:,:,nz]
        q[:,:,:,-2]   = q[:,:,:,nz-1]
        q[:,:,:,-3]   = q[:,:,:,nz-2]
        q[:,:,:,nz+1] = q[:,:,:,0]
        q[:,:,:,nz+2] = q[:,:,:,1]
        q[:,:,:,nz+3] = q[:,:,:,2]
    end
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

#3D Weno function
function weno5(nx,ny,nz,q,axis)

    #Swap axes as required NOTE: First index stores variables
    if (axis==1)
        #q = copy(q)
        n1,n2,n3 = nx,ny,nz
    elseif (axis==2)
        q = permutedims(q,[1,3,2,4])
        n1,n2,n3 = ny,nx,nz
    elseif (axis==3)
        q = permutedims(q,[1,4,3,2])
        n1,n2,n3 = nz,ny,nx
    else
        println("Error at Axes Weno")
    end

    #qL and qR
    qL = OffsetArray(zeros(5,n1+7,n2+7,n3+7),1:5,-3:n1+3,-3:n2+3,-3:n3+3)
    qR = OffsetArray(zeros(5,n1+7,n2+7,n3+7),1:5,-3:n1+3,-3:n2+3,-3:n3+3)

    eps  = 1e-6
    pweno= 2

    c0 = 1/6
    c1 = 13/12
    c2 = 1/4

    d0 = 1/10
    d1 = 3/5
    d2 = 3/10

    #Compute smoothness
    β0 = OffsetArray(zeros(5,n1+7,n2+7,n3+7),1:5,-3:n1+3,-3:n2+3,-3:n3+3)
    β1 = OffsetArray(zeros(5,n1+7,n2+7,n3+7),1:5,-3:n1+3,-3:n2+3,-3:n3+3)
    β2 = OffsetArray(zeros(5,n1+7,n2+7,n3+7),1:5,-3:n1+3,-3:n2+3,-3:n3+3)

    for i in -1:n1+1
        β0[:,i,:,:] = c1*(q[:,i-2,:,:]-2*q[:,i-1,:,:]+q[:,i,:,:]).^2 +
                      c2*(q[:,i-2,:,:]-4*q[:,i-1,:,:]+3*q[:,i,:,:]).^2
        β1[:,i,:,:] = c1*(q[:,i-1,:,:]-2*q[:,i,:,:]+q[:,i+1,:,:]).^2 +
                      c2*(q[:,i-1,:,:]-q[:,i+1,:,:]).^2
        β2[:,i,:,:] = c1*(q[:,i,:,:]-2*q[:,i+1,:,:]+q[:,i+2,:,:]).^2 +
                      c2*(3*q[:,i,:,:]-4*q[:,i+1,:,:]+q[:,i+2,:,:]).^2
    end

    for i in -1:n1
        #Positive reconstruction
        α0 = d0./(β0[:,i,:,:].+eps).^pweno
        α1 = d1./(β1[:,i,:,:].+eps).^pweno
        α2 = d2./(β2[:,i,:,:].+eps).^pweno

        w0 = α0./(α0+α1+α2)
        w1 = α1./(α0+α1+α2)
        w2 = α2./(α0+α1+α2)

        q0 = c0.*(2.0.*q[:,i-2,:,:].-7.0.*q[:,i-1,:,:].+11.0.*q[:,i,:,:])
        q1 = c0.*(-q[:,i-1,:,:]+5.0.*q[:,i,:,:].+2.0.*q[:,i+1,:,:])
        q2 = c0.*(2.0.*q[:,i,:,:].+5.0.*q[:,i+1,:,:].-q[:,i+2,:,:])

        qL[:,i,:,:] = w0.*q0 + w1.*q1 + w2.*q2


        #Negative reconstruction
        α0 = d0./(β2[:,i+1,:,:].+eps).^pweno
        α1 = d1./(β1[:,i+1,:,:].+eps).^pweno
        α2 = d2./(β0[:,i+1,:,:].+eps).^pweno

        w0 = α0./(α0+α1+α2)
        w1 = α1./(α0+α1+α2)
        w2 = α2./(α0+α1+α2)

        q0 = c0.*(2.0.*q[:,i+3,:,:].-7.0.*q[:,i+2,:,:].+11.0.*q[:,i+1,:,:])
        q1 = c0.*(-q[:,i+2,:,:].+5.0.*q[:,i+1,:,:].+2.0.*q[:,i,:,:])
        q2 = c0.*(2.0.*q[:,i+1,:,:].+5.0.*q[:,i,:,:].-q[:,i-1,:,:])

        qR[:,i,:,:] = w0.*q0 + w1.*q1 + w2.*q2

    end

    #Swap axes as required NOTE: First index stores variables
    if (axis==1)
        #Nothing required
    elseif (axis==2)
        qL = permutedims(qL,[1,3,2,4])
        qR = permutedims(qR,[1,3,2,4])
    elseif (axis==3)
        qL = permutedims(qL,[1,4,3,2])
        qR = permutedims(qR,[1,4,3,2])
    else
        println("Error at Axes Weno")
    end

    return qL,qR
end

#Flux calculation
function flux(nx,ny,nz,q,axis)

    γ = consts.γ
    ρ = q[1,:,:,:]
    u = q[2,:,:,:]./ρ
    v = q[3,:,:,:]./ρ
    w = q[4,:,:,:]./ρ
    e = q[5,:,:,:]./ρ
    p = (γ-1)*(q[5,:,:,:] - 0.5*ρ.*(u.^2+v.^2+w.^2))
    h = e + p./ρ

    if (axis==1)
        F = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)
        F[1,:,:,:] = ρ.*u
        F[2,:,:,:] = ρ.*(u.^2) + p
        F[3,:,:,:] = ρ.*u.*v
        F[4,:,:,:] = ρ.*u.*w
        F[5,:,:,:] = ρ.*u.*h

        return F

    elseif (axis==2)
        G = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)
        G[1,:,:,:] = ρ.*v
        G[2,:,:,:] = ρ.*(u.*v)
        G[3,:,:,:] = ρ.*(v.^2) + p
        G[4,:,:,:] = ρ.*v.*w
        G[5,:,:,:] = ρ.*v.*h

        return G
    elseif (axis==3)
        H = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)
        H[1,:,:,:] = ρ.*w
        H[2,:,:,:] = ρ.*(u.*w)
        H[3,:,:,:] = ρ.*w.*v
        H[4,:,:,:] = ρ.*(w.^2) + p
        H[5,:,:,:] = ρ.*w.*h

        return H
    else
        println("Error in Axes index at Flux")
    end
end

#Propagation speed calculation
function cs_weno(q,nx,ny,nz,axis)
    γ = consts.γ
    ρ = q[1,:,:,:]
    u = q[2,:,:,:]./ρ
    v = q[3,:,:,:]./ρ
    w = q[4,:,:,:]./ρ
    e = q[5,:,:,:]./ρ
    p = (γ-1)*(q[5,:,:,:] - 0.5*ρ.*(u.^2+v.^2+w.^2))

    #a  = sqrt.(γ*p./ρ)


    a = OffsetArray(zeros(nx+7,ny+7,nz+7),-3:nx+3,-3:ny+3,-3:nz+3)
    for k in -3:nz+3
    for j in -3:ny+3
    for i in -3:nx+3
        if (ρ[i,j,k]<0 || p[i,j,k]<0)
        @info i,j,k,p[i,j,k],ρ[i,j,k]
        end
        a[i,j,k] = sqrt(γ*p[i,j,k]/ρ[i,j,k])
    end
    end
    end


    r  = OffsetArray(zeros(nx+7,ny+7,nz+7),-3:nx+3,-3:ny+3,-3:nz+3)
    cs = OffsetArray(zeros(nx+7,ny+7,nz+7),-3:nx+3,-3:ny+3,-3:nz+3)

    if (axis==1)
        r = max.(abs.(u),abs.(u-a),abs.(u+a))

        for i in -1:nx
            cs[i,:,:]=  max.(abs.(r[i-2,:,:]),abs.(r[i-1,:,:]),abs.(r[i,:,:]),
                             abs.(r[i+1,:,:]),abs.(r[i+2,:,:]),abs.(r[i+3,:,:]))
        end

    elseif (axis==2)
        r = max.(abs.(u),abs.(u-a),abs.(u+a))

        for j in -1:ny
            cs[:,j,:]=  max.(abs.(r[:,j-2,:]),abs.(r[:,j-1,:]),abs.(r[:,j,:]),
                             abs.(r[:,j+1,:]),abs.(r[:,j+2,:]),abs.(r[:,j+3,:]))
        end

    elseif (axis==3)
        r= max.(abs.(u),abs.(u-a),abs.(u+a))

        for k in -1:nz
            cs[:,:,k]=  max.(abs.(r[:,:,k-2]),abs.(r[:,:,k-1]),abs.(r[:,:,k]),
                             abs.(r[:,:,k+1]),abs.(r[:,:,k+2]),abs.(r[:,:,k+3]))
        end
    else
        println("Error Axes Index at cs_weno")
    end

    return cs
end

#Rusonov Flux Calc
function rusanov_3d(q,qL,fluxL,qR,fluxR,nx,ny,nz,axis)
    cs = cs_weno(q,nx,ny,nz,axis)

    if (axis==1)
        FL = copy(fluxL)
        FR = copy(fluxR)

        F = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)

        for i in -1:nx
            for n in 1:5
                F[n,i,:,:] = 0.5*(FL[n,i,:,:]+FR[n,i,:,:]) + 0.5*cs[i,:,:].*(qL[n,i,:,:]-qR[n,i,:,:])
            end
        end

        return F
    elseif (axis==2)
        GL = copy(fluxL)
        GR = copy(fluxR)

        G = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)

        for j in -1:ny
            for n in 1:5
                G[n,:,j,:] = 0.5*(GL[n,:,j,:]+GR[n,:,j,:]) + 0.5*cs[:,j,:].*(qL[n,:,j,:]-qR[n,:,j,:])
            end
        end

        return G

    elseif (axis==3)
        HL = copy(fluxL)
        HR = copy(fluxR)

        H = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)

        for k in -1:nz
            for n in 1:5
                H[n,:,:,k] = 0.5*(HL[n,:,:,k]+HR[n,:,:,k]) + 0.5*cs[:,:,k].*(qL[n,:,:,k]-qR[n,:,:,k])
            end
        end

        return H
    else
        println("Axes Index Error at rusanov_3d")
    end
end

#RHS calculation
function rhsInv(nx,ny,nz,dx,dy,dz,q)

    r = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)

    #x-direction
    qLx,qRx = weno5(nx,ny,nz,q,1)
    FLx = flux(nx,ny,nz,qLx,1)
    FRx = flux(nx,ny,nz,qRx,1)
    Fx  = rusanov_3d(q,qLx,FLx,qRx,FRx,nx,ny,nz,1)

    #y-direction
    qLy,qRy = weno5(nx,ny,nz,q,2)
    FLy = flux(nx,ny,nz,qLy,2)
    FRy = flux(nx,ny,nz,qRy,2)
    Fy  = rusanov_3d(q,qLy,FLy,qRy,FRy,nx,ny,nz,2)

    #z-direction
    qLz,qRz = weno5(nx,ny,nz,q,3)
    FLz = flux(nx,ny,nz,qLz,3)
    FRz = flux(nx,ny,nz,qRz,3)
    Fz  = rusanov_3d(q,qLz,FLz,qRz,FRz,nz,ny,nz,3)


    for i in 0:nx
        r[:,i,:,:] = -(Fx[:,i,:,:]-Fx[:,i-1,:,:])./dx
    end

    for j in 0:ny
        r[:,:,j,:] = r[:,:,j,:] -(Fy[:,:,j,:]-Fy[:,:,j-1,:])./dy
    end

    for k in 0:nz
        r[:,:,:,k] = r[:,:,:,k] -(Fz[:,:,:,k]-Fz[:,:,:,k-1])./dz
    end

    return r
end

#RHS
function rhs(nx,ny,nz,dx,dy,dz,q)
    ri = rhsInv(nx,ny,nz,dx,dy,dz,q)
    #rv = rhsVis(nx,ny,nz,dx,dy,dz,q)
    r  = ri #+rv
    return r
end

#Time stepping RK3
function tvdrk3(nx,ny,nz,dx,dy,dz,q,dt,idbg)
    qq = copy(q)
    qn = copy(q)

    #First step
    expbc!(q,nx,ny,nz,idbg)
    r  = rhs(nx,ny,nz,dx,dy,dz,q)
    qq = q + dt*r

    #Second step
    expbc!(qq,nx,ny,nz,idbg)
    r  = rhs(nx,ny,nz,dx,dy,dz,qq)
    qq = 0.75*q + 0.25*qq + 0.25*dt*r

    #Third Step
    expbc!(qq,nx,ny,nz,idbg)
    r  = rhs(nx,ny,nz,dx,dy,dz,qq)
    qn = 1/3*q + 2/3*qq + 2/3*dt*r

    return qn
end


#Create a function for the ns3d run
function ns3d(cfl=0.5,nx=16,ny=16,nz=16,nitermax=10000,tend=1.0,idbg=0)

    #Setup required
    γ  = consts.γ
    Ma = 0.08
    time = 0.0

    #Create the grid
    x,y,z,dx,dy,dz = grid_init(nx,ny,nz,idbg)

    #Initialise
    q,pq_init = init_3d(Ma,x,y,z,nx,ny,nz,idbg)
    qnew = zeros(5,nx+7,ny+7,nz+7)
    qnew = OffsetArray(qnew,1:5,-3:nx+3,-3:ny+3,-3:nz+3)

    #Boundary Conditions
    expbc!(q,nx,ny,nz,idbg)

    #Calc_dt
    dt = calc_dt(cfl,γ,q,nx,ny,nz,dx,dy,dz)

    for niter in 1:nitermax
        dt = calc_dt(cfl,γ,q,nx,ny,nz,dx,dy,dz)
        if time+dt > tend
            dt = tend-time
        end

        expbc!(q,nx,ny,nz,idbg)

        qnew = tvdrk3(nx,ny,nz,dx,dy,dz,q,dt,idbg)

        expbc!(qnew,nx,ny,nz,idbg)

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
