#Import useful stuff
using Plots
using OffsetArrays
using WriteVTK
using NPZ
using YAML
using Printf

module consts
    γ = 1.4
    Ma= 0.08
    Pr= 0.72
    function calc_nu(T,iflag=1)
    if (iflag==1)
        Tref = 300
        νref = 1.716e-5
        S    = 110.4
        ν    = νref*(T.^1.5)/(Tref^1.5)*(Tref+S)./(T.+S)
        return ν
    elseif (iflag==2)
        Re = T
        ν  = 1.0/Re
    end
    end
end

#Create the grid
function grid_init(nx,ny,nz)
    #Intialise the grid
    lx = 2*pi
    ly = 2*pi
    lz = 2*pi

    dx = lx/(nx+1)
    dy = ly/(ny+1)
    dz = lz/(nz+1)


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
                x[i,j,k] = dx/2 + dx*i
                y[i,j,k] = dy/2 + dy*j
                z[i,j,k] = dz/2 + dz*k
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
                p = 1/(γ*Ma^2) +
                    (cos(2*xyz[1])+cos(2*xyz[2]))*(cos(2*xyz[3])+2.0)/16.0
                e = p/(ρ*(γ-1)) +
                    0.5*(u^2 + v^2 + w^2)

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

#Create function for writing output
function output_data(q,x,y,z,nx,ny,nz,fname)

    fname = "output/"*fname

    γ = consts.γ
    ρ = q[1,0:nx,0:ny,0:nz]
    u = q[2,0:nx,0:ny,0:nz]./ρ
    v = q[3,0:nx,0:ny,0:nz]./ρ
    w = q[4,0:nx,0:ny,0:nz]./ρ
    e = q[5,0:nx,0:ny,0:nz]./ρ
    p = ρ.*(γ-1).*(e-0.5.*(u.^2+v.^2+w.^2))
    c = sqrt.(γ.*p./ρ)

    vtkfile = vtk_grid(fname,x[0:nx,0:ny,0:nz],y[0:nx,0:ny,0:nz],z[0:nx,0:ny,0:nz])
    vtkfile["Velocity"] = (u,v,w)
    vtkfile["Density" ] = ρ
    vtkfile["Pressure"] = p
    vtkfile["Speed of Sound"] = c

    vtk_save(vtkfile)
    return vtkfile
end

#3D Fifth-Order Upwind
function upwind5(nx,ny,nz,q,axis)

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
        println("Error at Axes Upwind")
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

    w0 = 1/30
    w1 =-13/60
    w2 = 47/60
    w3 = 27/60
    w4 =- 1/20

    for i in -1:n1

        #Positive reconstruction
        qL[:,i,:,:] = w0.*q[:,i-2,:,:] +
                      w1.*q[:,i-1,:,:] +
                      w2.*q[:,i,:,:]   +
                      w3.*q[:,i+1,:,:] +
                      w4.*q[:,i+2,:,:]


        #Negative reconstruction
        qR[:,i,:,:] = w0.*q[:,i+3,:,:] +
                      w1.*q[:,i+2,:,:] +
                      w2.*q[:,i+1,:,:] +
                      w3.*q[:,i,:,:]   +
                      w4.*q[:,i-1,:,:]

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
        println("Error at Axes Upwind")
    end

    return qL,qR
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
function rhsInv(nx,ny,nz,dx,dy,dz,q,iflx)
    expbc!(q,nx,ny,nz)
    r = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)

    #x-direction
    if (iflx==1)
        qLx,qRx = weno5(nx,ny,nz,q,1)
    elseif (iflx==2)
        qLx,qRx = upwind5(nx,ny,nz,q,1)
    end
    FLx = flux(nx,ny,nz,qLx,1)
    FRx = flux(nx,ny,nz,qRx,1)
    Fx  = rusanov_3d(q,qLx,FLx,qRx,FRx,nx,ny,nz,1)

    #y-direction
    if (iflx==1)
        qLy,qRy = weno5(nx,ny,nz,q,2)
    elseif (iflx==2)
        qLy,qRy = upwind5(nx,ny,nz,q,2)
    end
    FLy = flux(nx,ny,nz,qLy,2)
    FRy = flux(nx,ny,nz,qRy,2)
    Fy  = rusanov_3d(q,qLy,FLy,qRy,FRy,nx,ny,nz,2)

    #z-direction
    if (iflx==1)
        qLz,qRz = weno5(nx,ny,nz,q,3)
    elseif (iflx==2)
        qLz,qRz = upwind5(nx,ny,nz,q,3)
    end
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
#Function find factors for derv
function dfac(a,fac)
    return 1.0./(a.*float(fac))
end

#Function to calc Viscous Flux
function rhsVis(nx,ny,nz,dx,dy,dz,q,Re)

    ρ = q[1,:,:,:]
    u = q[2,:,:,:]./ρ
    v = q[3,:,:,:]./ρ
    w = q[4,:,:,:]./ρ
    e = q[5,:,:,:]./ρ

    #Interpolation coefficients
    g = 1.0/180.0
    a = 245.0
    b =-75.0
    c = 10.0

    dx1,dx2,dx3,dx4,dx5,dx6 = dfac(dx,[1.0,2.0,3.0,4.0,5.0,6.0])
    dy1,dy2,dy3,dy4,dy5,dy6 = dfac(dy,[1.0,2.0,3.0,4.0,5.0,6.0])
    dz1,dz2,dz3,dz4,dz5,dz6 = dfac(dz,[1.0,2.0,3.0,4.0,5.0,6.0])

    #For interpolation (FV)
    gi = 1.0/60.0
    ai = 37.0
    bi =-8.0
    ci = 1.0

    #for cross derivatives
    gd = 1.0/10.0
    ad = 15.0
    bd =-6.0
    cd = 1.0

    Tref = 300.0
    Pr   = consts.Pr
    Ma   = consts.Ma
    γ    = consts.γ

    g1 = 1.0/Re
    g2 = (2.0/3.0)/Re
    g3 = -1.0/(Re*Pr*(Ma^2)*(γ-1.0))
    g4 = (γ-1)*γ*(Ma^2)

    cc = 110.4/Tref

    #Compute Temp and Viscocity
    te = g4*(e - 0.5*(u.^2+v.^2+w.^2))
    mu = (te.^(1.5))*(1.0+cc)./(te.+cc)


    #Allocate for the flux vector
    r = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)

    """
    #X-direction
    """
    #uy
    cd1 = calc_coeff_cd(u,gd,ad,bd,cd,dy2,dy4,dy6,nx,ny,nz,2,1)
    #uz
    cd2 = calc_coeff_cd(u,gd,ad,bd,cd,dz2,dz4,dz6,nx,ny,nz,3,1)
    #vy
    cd3 = calc_coeff_cd(v,gd,ad,bd,cd,dy2,dy4,dy6,nx,ny,nz,2,1)
    #wz
    cd4 = calc_coeff_cd(w,gd,ad,bd,cd,dz2,dz4,dz6,nx,ny,nz,3,1)

    #Compute vis flux at interfaces
    #ux,vx,wx,tx
    ux = calc_dd(u,g,a,b,c,dx1,dx3,dx5,nx,ny,nz,1)
    vx = calc_dd(v,g,a,b,c,dx1,dx3,dx5,nx,ny,nz,1)
    wx = calc_dd(w,g,a,b,c,dx1,dx3,dx5,nx,ny,nz,1)
    tx = calc_dd(te,g,a,b,c,dx1,dx3,dx5,nx,ny,nz,1)

    #uy,uz,vy,wz
    uy = calc_ddi(cd1,gi,ai,bi,ci,nx,ny,nz,1)
    uz = calc_ddi(cd2,gi,ai,bi,ci,nx,ny,nz,1)
    vy = calc_ddi(cd3,gi,ai,bi,ci,nx,ny,nz,1)
    wz = calc_ddi(cd4,gi,ai,bi,ci,nx,ny,nz,1)

    #uu,vv,ww,mux
    uu = calc_ddi(u,gi,ai,bi,ci,nx,ny,nz,1)
    vv = calc_ddi(v,gi,ai,bi,ci,nx,ny,nz,1)
    ww = calc_ddi(w,gi,ai,bi,ci,nx,ny,nz,1)
    mux= calc_ddi(mu,gi,ai,bi,ci,nx,ny,nz,1)

    txx= g2*mux.*(2.0*ux - vy - wz)
    txy= g1*mux.*(uy + vx)
    txz= g1*mux.*(uz + wx)
    qx = g3*mux.*tx

    vf = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)
    """
    Doubtful check indexes if doesnt work
    """
    vf[1,:,:,:] .= 0.0
    vf[2,:,:,:]  = txx
    vf[3,:,:,:]  = txy
    vf[4,:,:,:]  = txz
    vf[5,:,:,:]  = uu.*txx + vv.*txy + ww.*txz - qx

    #Compute RHS contribution (central and FV)
    for i in 0:nx
        r[:,i,:,:] = (vf[:,i,:,:]-vf[:,i-1,:,:])./dx
    end
    #Deallocate
    cd1,cd2,cd3,cd4 = nothing,nothing,nothing,nothing
    ux,vx,wx,tx     = nothing,nothing,nothing,nothing
    uy,uz,vy,wz     = nothing,nothing,nothing,nothing
    uu,vv,ww,mux    = nothing,nothing,nothing,nothing
    txx,txy,txz,qx  = nothing,nothing,nothing,nothing
    vf              = nothing

    """
    #Y-direction
    """

    #vx
    cd1 = calc_coeff_cd(v,gd,ad,bd,cd,dx2,dx4,dx6,nx,ny,nz,1,2)
    #vz
    cd2 = calc_coeff_cd(v,gd,ad,bd,cd,dz2,dz4,dz6,nx,ny,nz,3,2)
    #ux
    cd3 = calc_coeff_cd(u,gd,ad,bd,cd,dx2,dx4,dx6,nx,ny,nz,1,2)
    #wz
    cd4 = calc_coeff_cd(w,gd,ad,bd,cd,dz2,dz4,dz6,nx,ny,nz,3,2)

    #Compute vis flux at interfaces
    #uy,vy,wy,ty
    uy = calc_dd(u,g,a,b,c,dy1,dy3,dy5,nx,ny,nz,2)
    vy = calc_dd(v,g,a,b,c,dy1,dy3,dy5,nx,ny,nz,2)
    wy = calc_dd(w,g,a,b,c,dy1,dy3,dy5,nx,ny,nz,2)
    ty = calc_dd(te,g,a,b,c,dy1,dy3,dy5,nx,ny,nz,2)

    #vx,vz,ux,wz
    vx = calc_ddi(cd1,gi,ai,bi,ci,nx,ny,nz,2)
    vz = calc_ddi(cd2,gi,ai,bi,ci,nx,ny,nz,2)
    ux = calc_ddi(cd3,gi,ai,bi,ci,nx,ny,nz,2)
    wz = calc_ddi(cd4,gi,ai,bi,ci,nx,ny,nz,2)

    #uu,vv,ww,muy
    uu = calc_ddi(u,gi,ai,bi,ci,nx,ny,nz,2)
    vv = calc_ddi(v,gi,ai,bi,ci,nx,ny,nz,2)
    ww = calc_ddi(w,gi,ai,bi,ci,nx,ny,nz,2)
    muy= calc_ddi(mu,gi,ai,bi,ci,nx,ny,nz,2)

    #txy,tyy,tyz,qy
    txy= g1*muy.*(uy+vx)
    tyy= g2*muy.*(2.0*vy - ux - wz)
    tyz= g1*muy.*(vz + wy)
    qy = g3*muy.*ty

    vg = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)
    """
    Doubtful check indexes if doesnt work
    """
    vg[1,:,:,:] .= 0.0
    vg[2,:,:,:]  = txy
    vg[3,:,:,:]  = tyy
    vg[4,:,:,:]  = tyz
    vg[5,:,:,:]  = uu.*txy + vv.*tyy + ww.*tyz - qy

    #Compute RHS contribution (central and FV)
    for j in 0:ny
        r[:,:,j,:] = r[:,:,j,:] + (vg[:,:,j,:]-vg[:,:,j-1,:])./dy
    end
    #Deallocate
    cd1,cd2,cd3,cd4 = nothing,nothing,nothing,nothing
    uy,vy,wy,ty     = nothing,nothing,nothing,nothing
    vx,vz,ux,wz     = nothing,nothing,nothing,nothing
    uu,vv,ww,muy   = nothing,nothing,nothing,nothing
    txy,tyy,tyz,qy  = nothing,nothing,nothing,nothing
    vg              = nothing

    """
    #Z-direction
    """

    #wx
    cd1 = calc_coeff_cd(w,gd,ad,bd,cd,dx2,dx4,dx6,nx,ny,nz,1,3)
    #wy
    cd2 = calc_coeff_cd(w,gd,ad,bd,cd,dz2,dz4,dz6,nx,ny,nz,2,3)
    #ux
    cd3 = calc_coeff_cd(u,gd,ad,bd,cd,dx2,dx4,dx6,nx,ny,nz,1,3)
    #vy
    cd4 = calc_coeff_cd(v,gd,ad,bd,cd,dz2,dz4,dz6,nx,ny,nz,2,3)

    #Compute vis flux at interfaces
    #uz,vz,wz,tz
    uz = calc_dd(u,g,a,b,c,dy1,dy3,dy5,nx,ny,nz,3)
    vz = calc_dd(v,g,a,b,c,dy1,dy3,dy5,nx,ny,nz,3)
    wz = calc_dd(w,g,a,b,c,dy1,dy3,dy5,nx,ny,nz,3)
    tz = calc_dd(te,g,a,b,c,dy1,dy3,dy5,nx,ny,nz,3)

    #wx,wy,ux,vy
    wx = calc_ddi(cd1,gi,ai,bi,ci,nx,ny,nz,3)
    wy = calc_ddi(cd2,gi,ai,bi,ci,nx,ny,nz,3)
    ux = calc_ddi(cd3,gi,ai,bi,ci,nx,ny,nz,3)
    vy = calc_ddi(cd4,gi,ai,bi,ci,nx,ny,nz,3)

    #uu,vv,ww,muz
    uu = calc_ddi(u,gi,ai,bi,ci,nx,ny,nz,3)
    vv = calc_ddi(v,gi,ai,bi,ci,nx,ny,nz,3)
    ww = calc_ddi(w,gi,ai,bi,ci,nx,ny,nz,3)
    muz= calc_ddi(mu,gi,ai,bi,ci,nx,ny,nz,3)

    #txz,tyz,tzz,qz
    txz= g1*muz.*(uz+wx)
    tyz= g1*muz.*(vz+wy)
    tzz= g2*muz.*(2.0*wz - ux - vy)
    qz = g3*muz.*tz

    vh = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)
    """
    Doubtful check indexes if doesnt work
    """
    vh[1,:,:,:] .= 0.0
    vh[2,:,:,:]  = txz
    vh[3,:,:,:]  = tyz
    vh[4,:,:,:]  = tzz
    vh[5,:,:,:]  = uu.*txz + vv.*tyz + ww.*tzz - qz

    #Compute RHS contribution (central and FV)
    for k in 0:nz
        r[:,:,:,k] = r[:,:,:,k] + (vh[:,:,:,k]-vh[:,:,:,k-1])./dz
    end
    #Deallocate
    cd1,cd2,cd3,cd4 = nothing,nothing,nothing,nothing
    uz,vz,wz,tz     = nothing,nothing,nothing,nothing
    wx,wy,ux,vy     = nothing,nothing,nothing,nothing
    uu,vv,ww,muz    = nothing,nothing,nothing,nothing
    txz,tyz,tzz,qz  = nothing,nothing,nothing,nothing
    vh              = nothing

    return r
end

#Function to calculate value at interface
function calc_ddi(u,g,a,b,c,nx,ny,nz,axis)
    ddi = OffsetArray(zeros(nx+7,ny+7,nz+7),-3:nx+3,-3:ny+3,-3:nz+3)
    if (axis==1)
        for k in 0:nz
        for j in 0:ny
        for i in -1:nx
            ddi[i,j,k] =  g*(a*( u[i+1,j,k] + u[i,j,k])   +
                             b*( u[i+2,j,k] + u[i-1,j,k]) +
                             c*( u[i+3,j,k] + u[i-2,j,k]) )
        end
        end
        end
    elseif (axis==2)
        for k in 0:nz
        for j in -1:ny
        for i in 0:nx
            ddi[i,j,k] =  g*(a*( u[i,j+1,k] + u[i,j,k])   +
                             b*( u[i,j+2,k] + u[i,j-1,k]) +
                             c*( u[i,j+3,k] + u[i,j-2,k]) )
        end
        end
        end
    elseif (axis==3)
        for k in -1:nz
        for j in 0:ny
        for i in 0:nx
            ddi[i,j,k] =  g*(a*( u[i,j,k+1] + u[i,j,k])   +
                             b*( u[i,j,k+2] + u[i,j,k-1]) +
                             c*( u[i,j,k+3] + u[i,j,k-2]) )
        end
        end
        end
    else
        print("Error in Axes Index at calc_cd")
    end

    return ddi
end

#Function to calculate the derivates
function calc_dd(u,g,a,b,c,d1,d2,d3,nx,ny,nz,axis)
    dd = OffsetArray(zeros(nx+7,ny+7,nz+7),-3:nx+3,-3:ny+3,-3:nz+3)
    if (axis==1)
        for k in 0:nz
        for j in 0:ny
        for i in -1:nx
            dd[i,j,k] =  g*(a*d1*( u[i+1,j,k] - u[i,j,k]) +
                            b*d2*( u[i+2,j,k] - u[i-1,j,k]) +
                            c*d3*( u[i+3,j,k] - u[i-2,j,k]) )
        end
        end
        end
    elseif (axis==2)
        for k in 0:nz
        for j in -1:ny
        for i in 0:nx
            dd[i,j,k] =  g*(a*d1*( u[i,j+1,k] - u[i,j,k]) +
                            b*d2*( u[i,j+2,k] - u[i,j-1,k]) +
                            c*d3*( u[i,j+3,k] - u[i,j-2,k]) )
        end
        end
        end
    elseif (axis==3)
        for k in -1:nz
        for j in 0:ny
        for i in 0:nx
            dd[i,j,k] =  g*(a*d1*( u[i,j,k+1] - u[i,j,k]) +
                            b*d2*( u[i,j,k+2] - u[i,j,k-1]) +
                            c*d3*( u[i,j,k+3] - u[i,j,k-2]) )
        end
        end
        end
    else
        print("Error in Axes Index at calc_cd")
    end

    return dd
end

#Function to calculate the coefficients for cross derivates
function calc_coeff_cd(u,g,a,b,c,d1,d2,d3,nx,ny,nz,axis,flxdir)
    cd = OffsetArray(zeros(nx+7,ny+7,nz+7),-3:nx+3,-3:ny+3,-3:nz+3)
    if (flxdir==1)
        nxb,nxe,nyb,nye,nzb,nze = -3,nx+3,0,ny,0,nz
    elseif (flxdir==2)
        nxb,nxe,nyb,nye,nzb,nze = 0,nx,-3,ny+3,0,nz
    elseif (flxdir==3)
        nxb,nxe,nyb,nye,nzb,nze = 0,nx,0,ny,-3,nz+3
    end

    if (axis==1)
        for k in nzb:nze
        for j in nyb:nye
        for i in nxb:nxe
            cd[i,j,k] =  g*(a*d1*( u[i+1,j,k] - u[i-1,j,k]) +
                            b*d2*( u[i+2,j,k] - u[i-2,j,k]) +
                            c*d3*( u[i+3,j,k] - u[i-3,j,k]) )
        end
        end
        end
    elseif (axis==2)
        for k in nzb:nze
        for j in nyb:nye
        for i in nxb:nxe
            cd[i,j,k] =  g*(a*d1*( u[i,j+1,k] - u[i,j-1,k]) +
                            b*d2*( u[i,j+2,k] - u[i,j-2,k]) +
                            c*d3*( u[i,j+3,k] - u[i,j-3,k]) )
        end
        end
        end
    elseif (axis==3)
        for k in nzb:nze
        for j in nyb:nye
        for i in nxb:nxe
            cd[i,j,k] =  g*(a*d1*( u[i,j,k+1] - u[i,j,k-1]) +
                            b*d2*( u[i,j,k+2] - u[i,j,k-2]) +
                            c*d3*( u[i,j,k+3] - u[i,j,k-3]) )
        end
        end
        end
    else
        print("Error in Axes Index at calc_cd")
    end

    return cd
end

#RHS
function rhs(nx,ny,nz,dx,dy,dz,q,ivis,iflx,Re)
    ri = rhsInv(nx,ny,nz,dx,dy,dz,q,iflx)
    if (ivis==1)
        rv = rhsVis(nx,ny,nz,dx,dy,dz,q,Re)
        r  = ri + rv
    else
        r  = ri
    end

    return r
end

#Time stepping RK3
function tvdrk3(nx,ny,nz,dx,dy,dz,q,dt,ivis,iflx,Re)
    qq = copy(q)
    qn = copy(q)

    #First step
    expbc!(q,nx,ny,nz)
    r  = rhs(nx,ny,nz,dx,dy,dz,q,ivis,iflx,Re)
    qq = q + dt*r

    #Second step
    expbc!(qq,nx,ny,nz)
    r  = rhs(nx,ny,nz,dx,dy,dz,qq,ivis,iflx,Re)
    qq = 0.75*q + 0.25*qq + 0.25*dt*r

    #Third Step
    expbc!(qq,nx,ny,nz)
    r  = rhs(nx,ny,nz,dx,dy,dz,qq,ivis,iflx,Re)
    qn = 1/3*q + 2/3*qq + 2/3*dt*r

    return qn
end

#Compute Turbulent Kinetic Energy
function calc_tke(q,nx,ny,nz)
    ρ = q[1,0:nx,0:ny,0:nz]
    u = q[2,0:nx,0:ny,0:nz]./ρ
    v = q[3,0:nx,0:ny,0:nz]./ρ
    w = q[4,0:nx,0:ny,0:nz]./ρ

    tke = sum(0.5*(u.^2 + v.^2 + w.^2))./((nx+1)*(ny+1)*(nz+1))
    return tke
end


#Create a function for the ns3d run
function ns3d(cfl=0.5,nx=16,ny=16,nz=16,nitermax=10000,tend=1.0,nout=10,ivis=0,iflx=1,Re=1600)

    #Setup required
    γ  = consts.γ
    Ma = 0.08
    time = 0.0

    nx = nx - 1
    ny = ny - 1
    nz = nz - 1

    #Create the grid
    x,y,z,dx,dy,dz = grid_init(nx,ny,nz)

    #Initialise
    q,pq_init = init_3d(Ma,x,y,z,nx,ny,nz)
    @info calc_tke(q,nx,ny,nz)
    qnew = zeros(5,nx+7,ny+7,nz+7)
    qnew = OffsetArray(qnew,1:5,-3:nx+3,-3:ny+3,-3:nz+3)

    #Boundary Conditions
    expbc!(q,nx,ny,nz)

    #Calc_dt
    dt = calc_dt(cfl,γ,q,nx,ny,nz,dx,dy,dz)

    tke_old = calc_tke(q,nx,ny,nz)
    tkelist = append!(zeros(0),tke_old)
    dElist  = append!(zeros(0),0.0)
    tlist   = append!(zeros(0),0.0)

    #Initialise the PVD File
    if (!ispath("output"))
        mkdir("output")
    else
        rm("output/",recursive=true)
        mkdir("output")
    end
    pvd = paraview_collection("output/output_all")
    fname = "output_initial"
    vtkfile = output_data(q,x,y,z,nx,ny,nz,fname)
    pvd[time] = vtkfile

    #Open log file for residual
    io1 = open("residual.log","w")
    @printf(io1,"   iter \t dt \t time \n")
    for niter in 1:nitermax
        dt = calc_dt(cfl,γ,q,nx,ny,nz,dx,dy,dz)
        if time+dt > tend
            dt = tend-time
        end
        #tke_new = calc_tke(q,nx,ny,nz)
        expbc!(q,nx,ny,nz)

        qnew = tvdrk3(nx,ny,nz,dx,dy,dz,q,dt,ivis,iflx,Re)

        expbc!(qnew,nx,ny,nz)

        q = copy(qnew)
        time = time + dt

        @info niter,dt,time
        @printf(io1," %8i  %.3e %.3e \n",niter,dt,time)
        flush(io1)

        tke_new = calc_tke(q,nx,ny,nz)
        dEdt    = (tke_new - tke_old)/dt

        tkelist  = append!(tkelist,tke_new)
        dElist   = append!(dElist,dEdt)
        tlistc   = append!(tlist,time)
        if (time >= tend)
            break
        end
        tke_old = copy(tke_new)
        if (niter%nout ==0)
            fname = "output_"*string(Int(floor(niter/nout)))
            vtkfile = output_data(q,x,y,z,nx,ny,nz,fname)
            pvd[time] = vtkfile
        end
    end

    #Output data
    fname    = "output_final"
    vtkfile  = output_data(q,x,y,z,nx,ny,nz,fname)
    pvd[time]= vtkfile
    vtk_save(pvd)

    #Close any opened files
    close(io1)

    #Save NPZ file fpr the fields also
    fname = "data_fields.npz"

    γ = consts.γ
    ρ = q[1,:,:,:]
    u = q[2,:,:,:]./ρ
    v = q[3,:,:,:]./ρ
    w = q[4,:,:,:]./ρ
    e = q[5,:,:,:]./ρ
    p = (γ-1)*(q[5,:,:,:] - 0.5*ρ.*(u.^2+v.^2+w.^2))
    c = sqrt.(γ*p./ρ)

    npzwrite(fname, Dict("u" => u,"v" => v,"w" => w, "Density" => ρ, "Pressure" => p, "c" => c))

    return tkelist,dElist,tlist
end
#%%
inp_data = YAML.load_file("input.yaml")

#Read the input file and parameters
#Runtime parameters
nitermax= inp_data["nitermax"]
nx      = inp_data["nx"]
ny      = inp_data["ny"]
nz      = inp_data["nz"]
cfl     = inp_data["cfl"]
Re      = inp_data["Re"]
tend    = inp_data["tend"]
nplot   = inp_data["nplot"]

#Specific Flags
ivis = inp_data["ivis"]
iflx = inp_data["iflx"]

ihpc = inp_data["ihpc"]
tke,dEdt,tlist = ns3d(cfl,nx,nx,nx,nitermax,tend,nplot,ivis,iflx,Re)
npzwrite("data.npz", Dict("tkelist" => tke, "dEdt" => -dEdt, "tlist" => tlist))

#%%
# Operations for plotting only on local PC
if (ihpc==0)
    vars = npzread("data.npz")
    tkelist = vars["tkelist"]
    dEdt = vars["dEdt"]
    tlist = vars["tlist"]
    p1 = plot(tlist,tke)
    p2 = plot(tlist[2:lastindex(tlist)],dEdt[2:lastindex(dEdt)])
end
