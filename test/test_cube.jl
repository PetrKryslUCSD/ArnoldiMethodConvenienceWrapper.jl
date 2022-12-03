module unit_cube_tet_examples

using Test
using FinEtools
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear
using ArnoldiMethodConvenienceWrapper
using LinearAlgebra
using DataDrop

E = 1*phun("PA");
nu = 0.499;
rho = 1*phun("KG/M^3");
a = 1*phun("M"); b = a; h =  a;
n1 = 16;# How many element edges per side?
na =  n1; nb =  n1; nh  = n1;
                  # how many eigenvalues
OmegaShift = (0.01*2*pi)^2;


reffs = Dict(8 => [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.51509e-01, 2.56099e-01, 3.37083e-01, 3.39137e-01, 3.42637e-01, 3.48691e-01, 3.49059e-01, 3.50360e-01, 3.77198e-01, 3.96918e-01, 3.97757e-01, 4.27836e-01, 4.29423e-01, 4.30705e-01], 16 => [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,  2.59393e-01, 2.60642e-01, 3.51863e-01, 3.52430e-01, 3.53351e-01, 3.57428e-01, 3.57552e-01, 3.57918e-01, 3.99573e-01, 4.05492e-01, 4.05636e-01, 4.52967e-01, 4.53473e-01, 4.54334e-01], 32 => [0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.61727e-01, 2.62041e-01, 3.55982e-01, 3.56129e-01, 3.56361e-01, 3.59778e-01, 3.59813e-01, 3.59909e-01, 4.06075e-01, 4.07645e-01, 4.07676e-01, 4.59447e-01, 4.59559e-01, 4.60127e-01])

function unit_cube_esnice_ssit(N, neigvs = 20)
    tol = 0.001
    na,nb,nh = N, N, N
    # println("""
    # Vibration modes of unit cube  of almost incompressible material.

    # Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
    # tetrahedral. International Journal for Numerical Methods in
    # Engineering 67: 841-867.
    # """)

    MR = DeforModelRed3D
    fens,fes  = T4block(a,b,h, na,nb,nh)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)

    @info "size(M) = $(size(M))"
    # DataDrop.store_matrix("unit_cube_$N.h5", "/K", K)
    # DataDrop.retrieve_matrix("unit_cube_$N.h5", "/K")
    # DataDrop.store_matrix("unit_cube_$N.h5", "/M", M)


    @info "N=$(N), neigvs=$(neigvs), eigs"
    @time d,v,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    # println("Eigenvalues: $fs [Hz]")
    @test norm(fs - reffs[N]) / norm(reffs) < tol
    # reffs = fs

    # @info "N=$(N), neigvs=$(neigvs), eigs"
    # @time d, v, info = geneigsolve((K+OmegaShift*M, M), neigvs, :SR; krylovdim = 2*neigvs, issymmetric = true, verbosity = 1)
    # @show info
    # d = d .- OmegaShift;
    # fs = real(sqrt.(complex(d)))/(2*pi)
    # # println("Eigenvalues: $fs [Hz]")
    # @test norm(fs - reffs) / norm(reffs) < tol

    # @info "N=$(N), neigvs=$(neigvs), ssit"
    # @time d,v,nconv = ssit(K+OmegaShift*M, M; nev=neigvs, verbose=true)
    # d = d .- OmegaShift;
    # fs = real(sqrt.(complex(d)))/(2*pi)
    # # println("Eigenvalues: $fs [Hz]")
    # @test norm(fs - reffs) / norm(reffs) < tol

    # mode = 17
    # scattersysvec!(u, v[:,mode])
    # File =  "unit_cube_esnice.vtk"
    # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
    # @async run(`"paraview.exe" $File`)
    true
end # unit_cube_esnice

for N in (8, 16, 32)
    unit_cube_esnice_ssit(N, 20)
end


for N in (32, )
    unit_cube_esnice_ssit(N, 100)
    # unit_cube_esnice_ssit(N, 500)
end

end # module 
nothing