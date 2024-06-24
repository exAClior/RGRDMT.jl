using RGRDMT, Test
using CairoMakie
using DelimitedFiles

function my_plot(eng_filenames::Vector{String}, n_filenames::Vector{String})
    f = Figure()
    ax1 = Axis(f[1, 1],
        xlabel="n",
        xscale=log2,
        xticks=[2, 4, 6, 8, 10, 20, 30, 60],
        ylabel="ΔE_relax",
        yscale=log10,
        yticks=[1e-4, 1e-3, 1e-2, 1e-1],
        title="S = 1/2 Heisenberg",
    )
    xlims!(2, 60)
    ylims!(1e-4, 1e-1)
    # markers = (:cross, :diamond, :rect, :circle)
    plot_objs = []
    for (ii, (efname, nfname)) in enumerate(zip(eng_filenames, n_filenames))
        ΔE = readdlm(efname, ',')
        ΔE = reshape(ΔE, reduce(*, size(ΔE)))
        n = readdlm(nfname, ',')
        n = reshape(n, reduce(*, size(n)))
        pts = [Point2f(nn, δe) for (δe, nn) in zip(ΔE, n)]
        plt_obj = scatterlines!(ax1, pts)
        # plt_obj = scatterlines!(ax1, pts, marker=markers[ii])
        push!(plot_objs, plt_obj)
    end

    labelnames = vcat(["Reduced ρ"], ["Isometry, D=$D" for D in 2:2+length(n_filenames)-2])
    Legend(f[1, 2], plot_objs, labelnames)
    return f
end

# eng_filenames = vcat(["data/etfi.csv"], ["data/etfi$D.csv" for D in 2:3])

# n_filenames = vcat(["data/ntfi.csv"], ["data/ntfi$D.csv" for D in 2:3])


eng_filenames = vcat(["data/exxx.csv"], ["data/exxx$D.csv" for D in 2:3])

n_filenames = vcat(["data/nxxx.csv"], ["data/nxxx$D.csv" for D in 2:3])



cur_plt = my_plot(eng_filenames, n_filenames)

save("test_tfi.png", cur_plt, px_per_unit=2)
