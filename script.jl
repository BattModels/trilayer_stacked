using Pkg; Pkg.activate("/Users/mbabar/Desktop/PhD/Analysis/Gerischer/")

loc = "/Users/mbabar/Desktop/PhD/Analysis/Gerischer/ElectrochemicalKinetics.jl/src/"
push!(LOAD_PATH, loc)

# Load libs
using MAT
using CSV
using DelimitedFiles
using Glob
using Interpolations
using QuadGK
using Plots
using ElectrochemicalKinetics
import ElectrochemicalKinetics.compute_k_cq

# Load functions
function integrand(
    mhcd::MarcusHushChidseyDOS,
    ox::Bool;
    kT::Real = 0.026,
    η::Real = 0.0,
    V_q::Real = 0.0,
)
    function marcus_term(E)
        local exp_arg
        if ox
            exp_arg = -(( E .+ mhcd.λ) .^ 2) ./ (4 * mhcd.λ * kT)
        else
            exp_arg = -(( E .- mhcd.λ) .^ 2) ./ (4 * mhcd.λ * kT)
        end
        exp.(exp_arg)
    end
    fd(E) = ox ? 1 .- fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    E -> mhcd.A .* ((mhcd.dos.interp_func.(E .+ V_q)).^ 1) .* marcus_term(E .+ η) .* fd(E)
end

function compute_k_cq(
    η,
    model::MarcusHushChidseyDOS,
    ox::Bool;
    Eo = -0.07, # E_f,red (solvent) - E_f,vac (bilayer)
    C_dl = 10.0,
    Vq_min = -0.5,
    Vq_max = 0.5,
    kT = 0.026,
    E_min = model.dos.E_min,
    E_max = model.dos.E_max,
)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    #Vappl_data, Vdl_data = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    #v_interp = LinearInterpolation(Vappl_data, Vdl_data)
    V_t = (Eo + η)
    V_dl = V_dl_interp(V_t)
    
    #V_dl = v_interp(V_t)

    V_q = V_t - V_dl
    if V_q < 0
        E_max = E_max .- 0.05
        E_min = E_min .- V_q .+ 0.05
    elseif V_q > 0
        E_min = E_min .+ 0.05
        E_max = E_max  .- V_q .- 0.05
    end
    #print(E_min, E_max)
    k_rate = quadgk(integrand(model, ox; kT = kT, η = η, V_q = V_q), E_min, E_max)[1]
    return k_rate, V_q
end

#Prefactor
function prefactor(
    mhcd::MarcusHushChidseyDOS;
    kT::Real = 0.026,
    V_q = 0.0,
)
    #fd(E) = ox ? 1 .- fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    Ft(E) = (1/4kT).*((sech.(E./2kT)).^2) # Thermal broadening function
    E -> ((mhcd.dos.interp_func.(E .+ V_q)).^ 1) .* Ft(E)
end

# Modify string
function chop_str(str::String)
         while str[length(str)] == '0'
               str = chop(str)
         end
         if str[length(str)] == '.'
            str = chop(str)
         end
         return str
end

# Get rate k
#dos = [transpose(Elist) transpose(dos_tot)]
## Define model
Vq_max = 0.6
Vq_min = -0.6
C_dl = 10.0
λ = 0.82
Eo = -0.07
A = 1.0
η = 0.0
kT = 0.026 #eV
ef = 0.0 # Fermi level eV

#η_list = LinRange(-0.3, 0.3, 13)
#kox_data = zeros(Float64, length(flist));
#kred_data = zeros(Float64, length(flist));
#q_list = zeros(Float64, length(flist));
#mhcd = MarcusHushChidseyDOS(A, λ, dos)
#k_ox = compute_k_cq(0.4, mhcd, true; Vq_min=Vq_min, Vq_max=Vq_max)
#k_red = compute_k_cq(0.4, mhcd, false; Vq_min=Vq_min, Vq_max=Vq_max)

# ABA
file = matopen("ABA_dos.mat")
Elist = transpose(read(file, "E_list"));
dos_tot = read(file, "dos");
dos_data = [Elist dos_tot]
mhcd = MarcusHushChidseyDOS(A, λ, dos_data, Ef=ef)
kox,V_q = compute_k_cq(η, mhcd, true; Eo=Eo, Vq_min=Vq_min, Vq_max=Vq_max);
kred,V_q = compute_k_cq(η, mhcd, false; Eo=Eo, Vq_min=Vq_min, Vq_max=Vq_max);
print("ABA ","k_ox ",kox," k_red ",kred,"\n")

# ABC
file = matopen("ABC_dos.mat")
Elist = transpose(read(file, "E_list"));
dos_tot = read(file, "dos");
dos_data = [Elist dos_tot]
mhcd = MarcusHushChidseyDOS(A, λ, dos_data, Ef=ef)
kox,V_q = compute_k_cq(η, mhcd, true; Eo=Eo, Vq_min=Vq_min, Vq_max=Vq_max);
kred,V_q = compute_k_cq(η, mhcd, false; Eo=Eo, Vq_min=Vq_min, Vq_max=Vq_max);
print("ABC ","k_ox ",kox," k_red ",kred,"\n")

#q_list = q_list[sortperm(q_list)]
#kox_data = kox_data[sortperm(q_list)]
#kred_data = kred_data[sortperm(q_list)]

# Write
#file = matopen("krate_tmblg.mat","w")
#write(file,"theta_list",q_list)
#write(file, "kox_data",kox_data)
#write(file, "kred_data",kred_data)

# Plots
#plot(η_list, abs.(kox_data- kred_data), label="AAA")
#xlabel!("η (eV)")
#ylabel!("kox - kred")

# Compute prefactor
#factor[i] = quadgk(prefactor(mhcd; kT = kT, V_q = V_q), -0.45, 0.45)[1]
#print(i,",",V_q,",",factor[i],"\n")


#print(k)


##
#@time begin
#     mhcd = MarcusHushChidseyDOS(20.0, 0.82, dos)
#     k = compute_k(0.4, mhcd; calc_cq=true, Vq_min=-0.45, Vq_max=0.45)
#end


