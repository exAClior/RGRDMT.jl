using RGRDMT, Test
using CairoMakie

function my_plot(eng_filenames::Vector{String}, n_filenames::Vector{String})
    f = Figure()
    ax1 = Axis(f[1, 1],
        xlabel="n",
        xscale=log2,
        xticks=[0, 2, 4, 6, 8, 10, 20],
        ylabel="ΔE",
        yscale=log10,
        yticks=[1e-4, 1e-3, 1e-2, 1e-1, 1],
        title="Lower Bound",
    )
    markers = (:cross, :diamond, :rect, :circle)
    plot_objs = []
    for (ii, (efname, nfname)) in enumerate(zip(eng_filenames, n_filenames))
        ΔE = readdlm(efname, ',')
        ΔE = reshape(ΔE, reduce(*, size(ΔE)))
        n = readdlm(nfname, ',')
        n = reshape(n, reduce(*, size(n)))
        pts = [Point2f(nn, δe) for (δe, nn) in zip(ΔE, n)]
        plt_obj = scatter!(ax1, pts, marker=markers[ii])
        push!(plot_objs, plt_obj)
    end

    labelnames = ["Reduced ρ", "Isometry, D=2"]
    Legend(f[1, 2], plot_objs, labelnames)
    return f
end

# eng_filenames = ["data/etfi.csv", "data/etfi2.csv"]
# n_filenames = ["data/ntfi.csv", "data/ntfi2.csv"]

eng_filenames = ["data/exxx.csv", "data/exxx2.csv"]
n_filenames = ["data/nxxx.csv", "data/nxxx2.csv"]

cur_plt = my_plot(eng_filenames, n_filenames)
