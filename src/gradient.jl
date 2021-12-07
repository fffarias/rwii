module gradient

using SeisAcoustic

function conventional_gradient(res::Recordings, path_bnd::Ts, path_lwfd::Ts,
                               src::Source, params::TdParams) where {Ts <: String}

    # model length
    N = params.nz * params.nx

    # allocate memory for computing adjoint wavefield
    spt1 = Snapshot(params)
    spt2 = Snapshot(params)
    tmp    = zeros(params.data_format, params.Nz * params.Nx)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # allocate memory for reconstructing source-side wavefield backward
    wfd1 = Wavefield(params)
    wfd2 = read_one_wavefield(path_lwfd, 1)
    wfd_z1 = zeros(params.data_format, params.nz)
    wfd_z2 = zeros(params.data_format, params.nz)
    wfd_x1 = zeros(params.data_format, params.nx)
    wfd_x2 = zeros(params.data_format, params.nx)

    # initialize the boundary value as zero
    bnd = WavefieldBound(params)
    fid_bnd = open(path_bnd, "r")

    # initialize gradient to be zeros
    g = zeros(params.data_format, N)

    # backward time stepping
    for it = params.nt : -1 : 2

        # compute adjoint wavefield
        one_step_adjoint!(spt1, spt2, params, tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        inject_rec2spt!(spt1, res, it)

        # reconstruct source-side wavefield
        subtract_source!(wfd2, src, it)
        read_one_boundary!(bnd, fid_bnd, it-1, params)
        one_step_backward!(wfd1, wfd2, bnd, params, wfd_z1, wfd_z2, wfd_x1, wfd_x2)

        # compute the gradient
        apply_image_condition!(g, spt1, wfd1, wfd2, params)

        # prepare for the next step backward reconstruction
        copy_snapshot!(spt2, spt1)
        copy_wavefield!(wfd2, wfd1)

    end

    # close the boundary value file
    close(fid_bnd)

    # return the gradient of velocity model
    return g
end

function apply_image_condition!(g::Vector{Tv}, spt::Snapshot, wfd1::Wavefield,
                                wfd2::Wavefield, params::TdParams) where {Tv <: AbstractFloat}

    # number of elements
    N = params.nz * params.nx

    # update gradient by current adjoint snapshot
    for i = 1 : N
        p    = wfd2.p[i] 
        j    = params.spt2wfd[i]
        g[i] = g[i] + p * (spt.px[j] + spt.pz[j])
    end

    return nothing
end

end # module