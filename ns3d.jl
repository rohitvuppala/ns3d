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

function weno5(nx,ny,nz,q)
	#qL and qR
	qL = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)
	qR = OffsetArray(zeros(5,nx+7,ny+7,nz+7),1:5,-3:nx+3,-3:ny+3,-3:nz+3)

	eps  = 1e-6
	pweno= 2

	c0 = 1/6
	c1 = 13/12
	c2 = 1/4

	d0 = 1/10
	d1 = 3/5
	d2 = 3/10

	#Compute smoothness
	β0 = OffsetArray(zeros(5,nx+7),1:5,-3:nx+3)
	β1 = OffsetArray(zeros(5,nx+7),1:5,-3:nx+3)
	β2 = OffsetArray(zeros(5,nx+7),1:5,-3:nx+3)

	for i in 0:nx
		β0[:,i] = c1*(q[:,i-2]-2*q[:,i-1]+q[:,i]).^2
		        + c2*(q[:,i-2]-4*q[:,i-1]+3*q[:,i]).^2
		β1[:,i] = c1*(q[:,i-1]-2*q[:,i]+q[:,i+1]).^2
		        + c2*(q[:,i-1]-q[:,i+1]).^2
		β2[:,i] = c1*(q[:,i]-2*q[:,i+1]+q[:,i+2]).^2
		        + c2*(3*q[:,i]-4*q[:,i+1]+q[:,i+2]).^2
	end

	for i in 0:nx-1
		#Positive reconstruction
		α0 = d0./(β0[:,i].+eps).^pweno
		α1 = d1./(β1[:,i].+eps).^pweno
		α2 = d2./(β2[:,i].+eps).^pweno

		w0 = α0./(α0+α1+α2)
		w1 = α1./(α0+α1+α2)
		w2 = α2./(α0+α1+α2)

		q0 = c0.*(2.0.*q[:,i-2].-7.0.*q[:,i-1].+11.0.*q[:,i])
		q1 = c0.*(-q[:,i-1]+5.0.*q[:,i].+2.0.*q[:,i+1])
		q2 = c0.*(2.0.*q[:,i].+5.0.*q[:,i+1].-q[:,i+2])

		qL[:,i] = w0.*q0 + w1.*q1 + w2.*q2
		#@info qL[:,i]

		#Negative reconstruction
		α0 = d0./(β2[:,i+1].+eps).^pweno
		α1 = d1./(β1[:,i+1].+eps).^pweno
		α2 = d2./(β0[:,i+1].+eps).^pweno

		w0 = α0./(α0+α1+α2)
		w1 = α1./(α0+α1+α2)
		w2 = α2./(α0+α1+α2)

		q0 = c0.*(2.0.*q[:,i+3].-7.0.*q[:,i+2].+11.0.*q[:,i+1])
		q1 = c0.*(-q[:,i+2].+5.0.*q[:,i+1].+2.0.*q[:,i])
		q2 = c0.*(2.0.*q[:,i+1].+5.0.*q[:,i].-q[:,i-1])

		qR[:,i] = w0.*q0 + w1.*q1 + w2.*q2

	end

	return qL,qR
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
function ns3d(cfl=0.5,nx=32,ny=32,nz=32)

    #Setup required
    γ  = consts.γ
    Ma = 0.08
    time = 0.0

    #Create the grid
    x,y,z,dx,dy,dz = grid_init(nx,ny,nz)

    #Initialise
    q,pq_init = init_3d(Ma,x,y,z,nx,ny,nz)

    #Calc_dt
    dt = calc_dt(cfl,γ,q,nx,ny,nz,dx,dy,dz)

    #Boundary Conditions
    expbc!(q,nx,ny,nz)

    #Output data
    output_data(q,x,y,z,nx,ny,nz)

    return nothing
end
