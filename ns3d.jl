#Import useful stuff
using Plots
using OffsetArrays

#Create the grid
function grid_init(nx,ny,nz)
    #Intialise the grid
    lx = 2*pi
    ly = 2*pi
    lz = 2*pi

    dx = lx/nx
    dy = ly/ny
    dz = lz/nz

    x = zeros(nx+5,ny+5,nz+5)
    y = zeros(nx+5,ny+5,nz+5)
    z = zeros(nx+5,ny+5,nz+5)

    #Offset arrays to reflect ghost points
    x = OffsetArray(x,-2:nx+2,-2:ny+2,-2:nz+2)
    y = OffsetArray(y,-2:nx+2,-2:ny+2,-2:nz+2)
    z = OffsetArray(z,-2:nx+2,-2:ny+2,-2:nz+2)

    #Set the values for x,y,z and create grid
    for k in -2:nz+2
        for j in -2:ny+2
            for i in -2:nx+2
                x[i,j,k] = 0.0 + dx*i
                y[i,j,k] = 0.0 + dy*j
                z[i,j,k] = 0.0 + dz*k
            end
        end
    end

    return x,y,z
end

#Initialise
function init_3d(Ma,γ,x,y,z,nx,ny,nz)

    # To store cons ρ,ρu,ρv,ρw,ρe
    q = zeros(5,nx+5,ny+5,nz+5)
    q = OffsetArray(q,1:5,-2:nx+2,-2:ny+2,-2:nz+2)

    # To store ρ,u,v,w,p,e
    pq= zeros(6,nx+5,ny+5,nz+5)
    pq= OffsetArray(pq,1:6,-2:nx+2,-2:ny+2,-2:nz+2)

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

#Create a function for the ns3d run
function ns3d(nx=32,ny=32,nz=32)

    #Setup required
    Ma= 0.08
    γ = 1.4

    #Create the grid
    x,y,z = grid_init(nx,ny,nz)

    #Initialise
    q,pq = init_3d(Ma,γ,x,y,z,nx,ny,nz)

    return q
end
