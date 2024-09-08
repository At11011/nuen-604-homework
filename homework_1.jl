using DifferentialEquations
using Plots
using PhysicalConstants.CODATA2018
using Unitful
using DataFrames

function p6_numeric()
    # Define the parameters

            # Define the initial conditions
            N0_A = 1.0  # initial amount of A
            N0_B = 0.0     # initial amount of B
            N0_C = 0.0     # initial amount of C
    
            λ = log(2)/1.0
            
            # Define the system of ODEs
            function decay_chain!(du, u, p, t)
                NA, NB, NC = u
                λ = p
                du[1] = -λ * NA        # dNA/dt
                du[2] = λ * NA - λ * NB  # dNB/dt
                du[3] = λ * NB        # dNC/dt
            end

            # Set up the ODE problem
            u0 = [N0_A, N0_B, N0_C]
            tspan = (0.0, 10.0)
            p = (λ)
            prob = ODEProblem(decay_chain!, u0, tspan, p)

            # Solve the ODE
            sol = solve(prob)

            # Plot the results
            plot!(
            sol,
            label = ["N_A" "N_B" "N_C"],
            xlabel = "Time",
            ylabel = "Amount",
            title = "Radioactive Decay Chain",
            )
        end

        function KE_rel(error, m)
            p = error / 50 - 3
            q = 2 - error / 50
            k = 1
            gamma =
            (2 * sqrt(-p / 3) * cos((1 / 3) * acos((3q / (2p)) * sqrt(-3 / p)) - 2π * k / 3))^-1

            v = sqrt(SpeedOfLightInVacuum^2 * (1 - (1 / gamma)^2))
            KE = (gamma - 1) * m * SpeedOfLightInVacuum^2
            return KE, v
        end

        function p6()
            N_A0 = 1
            T_12 = 1
            λ = log(2) / T_12
            N_A(t) = N_A0 * exp(-t * λ)
            #  N_B(t) = λ * N_A0 * (1 / (λ - 1) * exp(-t) + 1/(1 - λ) * exp(-t* λ))
            # N_B(t) = -N_A0 * (exp(λ*(1 - 2t)) - exp(λ*(1 - t)))
            N_B(t) = λ * N_A0 * t * exp(-λ*t)
            t = (0:0.01:10) ./ λ
            N_Ad = N_A.(t)
            N_Bd = N_B.(t)

            plot(t, N_Ad, xlim = (0, 10),label = "N_A(t)", xlabel = "time (1/λ)", ylabel = "Activity", title = "Activity vs time", style = :dash)
            plot!(t, N_Bd, label = "N_B(t)", style = :dash)
            savefig("~/Code/LaTeX/Fall 2024/NUEN 604/NUEN 604 Homework 1/assets/problem_6_plot.png")

        end 

        function p7_eval(error, m)
            m = ElectronMass
            c = SpeedOfLightInVacuum
            KErel, v = KE_rel(error, m)

            KEnewt = (1 / 2) * m * v^2
            KErel = ((1 / sqrt(1 - v^2 / c^2)) - 1) * m * c^2

            println("KErel: $(uconvert(u"eV", KErel))")
            println("KEnewt: $(uconvert(u"eV", KEnewt))")
            println("Velocity: $v")
            println("Error: $(100 * (KErel - KEnewt) / KErel)")

        end

        function p9()
            material = ["O", "Si", "Ca", "Al", "K", "Na", "Fe", "H", "Mg", "S"]
            Mᵢ = [15.9994, 8.0855, 40.078, 26.981538, 39.0983, 22.98977, 55.845, 1.00794, 24.305, 32.066]u"g/mol"
            wᵢ = [0.4983, 0.3158, 0.0826, 0.0456, 0.0192, 0.0171, 0.0122, 0.0056, 0.0024, 0.0012] 
            σₐⁱ = [0.00028, 0.168, 0.43, 0.23, 2.1, 0.53, 2.56, 0.333, 0.066, 0.52]u"b"
            σₛⁱ = [4.2, 1.7, 3.0, 1.4, 1.5, 4.0, 11.0, 38.0, 3.6, 1.1]u"b"
            df = DataFrame((material = material, Mᵢ = Mᵢ, wᵢ = wᵢ, σₐⁱ = σₐⁱ, σₛⁱ = σₛⁱ))
            ρ = 2.35u"g/cm^3"

            Σₜᵐⁱˣ = sum(@. df.wᵢ * ρ * AvogadroConstant / df.Mᵢ * (df.σₐⁱ + df.σₛⁱ))
            Σₐᵐⁱˣ = sum(@. df.wᵢ * ρ * AvogadroConstant / df.Mᵢ * df.σₐⁱ)

            println("Total cross-section: $(uconvert(u"cm^-1", Σₜᵐⁱˣ))")
            println("Absorption cross-section: $(uconvert(u"cm^-1", Σₐᵐⁱˣ))")
            println("Absorption probability: $(Σₐᵐⁱˣ/Σₜᵐⁱˣ)")
            return Σₜᵐⁱˣ
        end


        function p10()
            compounds = ["B", "B₄C", "BN", "ZrB₂", "TiB₂", "HfB₂", "Hf"]
            Mᵢ = [10.811, 55.24, 24.81, 112.85, 69.5, 200.1, 178.5]u"g/mol"
            ρ = [2330, 2510, 2250, 6090, 4520, 11200, 13090]u"kg/m^3"
            w = [  
                    [10.811u"g/mol"]./10.811u"g/mol", # B
                    [4 * 10.811u"g/mol", 12.011u"g/mol"]./(4 * 10.811u"g/mol" + 12.011u"g/mol"), # B₄C
                    [10.811u"g/mol", 14.007u"g/mol"]./(10.811u"g/mol" + 14.007u"g/mol"), # BN
                    [91.224u"g/mol", 2 * 10.811u"g/mol"]./(2 * 10.811u"g/mol" + 91.224u"g/mol"), # ZrB₂
                    [47.867u"g/mol", 2* 10.811u"g/mol" ]./(47.867u"g/mol" + 2* 10.811u"g/mol"), # TiB₂
                    [178.486u"g/mol", 2 * 10.811u"g/mol"]./(2 * 10.811u"g/mol" + 178.486u"g/mol"), # HfB₂
                    [178.486u"g/mol"]./178.486u"g/mol" # Hf
                    ] 
            σₐ = [790e-24, 105e-24]u"cm^2"
            Σₐ = zeros(7)u"cm^-1"
            Σₐ[1] = sum(@. w[1] * ρ[1] * AvogadroConstant * σₐ[1] / Mᵢ[1]) # B
            Σₐ[2] = sum(@. w[2] * ρ[2] * AvogadroConstant * [σₐ[1], 0u"cm^2"] / Mᵢ[2]) # B₄C
            Σₐ[3] = sum(@. w[3] * ρ[3] * AvogadroConstant * [σₐ[1], 0u"cm^2"] / Mᵢ[3]) # BN
            Σₐ[4] = sum(@. w[4] * ρ[4] * AvogadroConstant * [0u"cm^2", σₐ[1]] / Mᵢ[4]) # ZrB₂
            Σₐ[5] = sum(@. w[5] * ρ[5] * AvogadroConstant * [0u"cm^2", σₐ[1]] / Mᵢ[5]) # TiB₂
            Σₐ[6] = sum(@. w[6] * ρ[6] * AvogadroConstant * [σₐ[2], σₐ[1]] / Mᵢ[6]) # HfB₂
            Σₐ[7] = sum(@. w[7] * ρ[7] * AvogadroConstant * σₐ[2] / Mᵢ[7]) # Hf

            println("Compound\tCross-section\tMean distance")
            for (compound, crosssection) in zip(compounds, Σₐ)
                println("$compound\t\t$(round(typeof(1.0u"cm^-1"), crosssection, digits=3))\t$(round(typeof(1.0u"cm"), 1/crosssection, digits=3))")
            end

            attn = 0.2
            println("Distance (cm): $(-log(attn) / Σₐ[1])")
        end
