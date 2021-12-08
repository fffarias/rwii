module utils

using SeisAcoustic, PyPlot, PyCall

function MultiSeisPlotTX(d::Vector{Matrix{T}}; style="color",
                      cmap="PuOr", pclip=98, vmin="NULL", vmax="NULL",
                      aspect="auto", interpolation="Hanning",
                      wiggle_fill_color="k", wiggle_line_color="k",
                      wiggle_trace_increment=1, xcur=1.2, scal="NULL",
                      title::Array{Ts,1}, titlesize=16, xlabel=" ", xunits=" ",
                      ylabel=" ", yunits=" ", labelsize=14, ox=0, dx=1,
                      oy=0, dy=1, xticks="NULL", yticks="NULL",
                      xticklabels="NULL", yticklabels="NULL", ticksize=11,
                      fignum="NULL", wbox=6, hbox=6, dpi=100, name="NULL") where {Ts<:String, T<:Real}

    if (vmin=="NULL" || vmax=="NULL")
        if (pclip<=100)
            a = -quantile(abs.(d[:]), (pclip/100))
        else
        a = -quantile(abs.(d[:]), 1)*pclip/100
        end
        b = -a
    else
        a = vmin
        b = vmax
    end
    plt.ion()
    if (fignum == "NULL")
        fig = plt.figure(figsize=(3 * wbox, 1 * hbox), dpi=dpi, facecolor="w",
                         edgecolor="k")
    else
        fig = plt.figure(num=fignum, figsize=(3 * wbox, 1 * hbox), dpi=dpi,
                         facecolor="w", edgecolor="k")
    end

    if (style != "wiggles")
        Tot = 1 * 3 #nrows * ncols
        Position = collect(1:Tot)
        fig = plt.figure(1,figsize=(3 * wbox, 1 * hbox), dpi=dpi, facecolor="w",
                edgecolor="k")
        for k=1:Tot
            # add every single subplot to the figure with a for loop
            ax = fig.add_subplot(1,3,Position[k])
            ax.set_title(title[k], fontsize=titlesize)                            
            ax.imshow(d[k], cmap=cmap, vmin=a, vmax=b,
                    extent=[ox - dx/2,ox + (size(d,2)-1)*dx + dx/2,
                            oy + (size(d,1)-1)*dy,oy],
                    aspect=aspect, interpolation=interpolation)

        end                            
        fig.tight_layout(pad=3.0)
    end
    if (style != "color")
        style=="wiggles" ? margin = dx : margin = dx/2
        y = oy .+ dy*collect(0:1:size(d, 1)-1)
        x = ox .+ dx*collect(0:1:size(d, 2)-1)
        delta = wiggle_trace_increment*dx
        hmin = minimum(x)
        hmax = maximum(x)
        dmax = maximum(abs.(d[:]))
        alpha = xcur*delta
        scal=="NULL" ? alpha = alpha/maximum(abs.(d[:])) : alpha=alpha*scal
        for k = 1:wiggle_trace_increment:size(d, 2)
            x_vert = Float64[]
            y_vert = Float64[]
            sc = x[k] * ones(size(d, 1))
            s  = d[:,k]*alpha + sc
            imm = plt.plot( s, y, wiggle_line_color)
            if (style != "overlay")
                plt.fill_betweenx(y, sc, s, where=s.>sc, facecolor=wiggle_line_color)
            end
        end
        plt.axis([ox - margin, ox + (size(d, 2)-1)*dx + margin,
        oy + (size(d, 1)-1)*dy, oy])
    end


    plt.xlabel(join([xlabel " " xunits]), fontsize=labelsize)
    plt.ylabel(join([ylabel " " yunits]), fontsize=labelsize)
    xticks == "NULL" ? nothing : plt.xticks(xticks)
    yticks == "NULL" ? nothing : plt.yticks(yticks)
    ax = plt.gca()
    xticklabels == "NULL" ? nothing : ax.set_xticklabels(xticklabels)
    yticklabels == "NULL" ? nothing : ax.set_yticklabels(yticklabels)
    plt.setp(ax.get_xticklabels(), fontsize=ticksize)
    plt.setp(ax.get_yticklabels(), fontsize=ticksize)
    if (name == "NULL")
        plt.show()
    else
        plt.savefig(name, dpi=dpi)
        plt.close()
    end
    return nothing # imm
end   

function readFile!(fname::AbstractString,
    data::Array{<:AbstractFloat,1})
    fp = open(fname, "r")
    read!(fp, data)
    close(fp)
end

end # module