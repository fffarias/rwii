module backward

using SeisAcoustic

function reduced_wavefield_imaging(res::Recordings, path_bnd::Ts, path_lwfd::Ts,
    src::Source, alpha::Tv, params::TdParams) where {Ts<:String, Tv<:AbstractFloat}

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
    wfd2 = read_one_wavefield(path_lwfd, 1) # read last wavefield
    wfd_z1 = zeros(params.data_format, params.nz)
    wfd_z2 = zeros(params.data_format, params.nz)
    wfd_x1 = zeros(params.data_format, params.nx)
    wfd_x2 = zeros(params.data_format, params.nx)

    # strength of source-side wavefield
    strength = zeros(params.data_format, N)

    # initialize the boundary value as zero
    bnd = WavefieldBound(params)
    fid_bnd = open(path_bnd, "r")

    # alpha 
    res.p[:] = alpha*res.p[:]

    # backward time stepping
    for it = params.nt : -1 : 2

        # compute adjoint wavefield
        one_step_adjoint!(spt1, spt2, params, tmp, tmp_z1, tmp_z2, tmp_x1, tmp_x2)
        inject_rec2spt!(spt1, res, it)

        # reconstruct source-side wavefield
        # remove dt * src[it]
        subtract_source!(wfd2, src, it)
        read_one_boundary!(bnd, fid_bnd, it-1, params)        
        one_step_backward!(wfd1, wfd2, bnd, params, wfd_z1, wfd_z2, wfd_x1, wfd_x2)

        # compute the strength of source-side wavefield
        reduce_source_strength!(strength, wfd1, wfd2, spt1, spt2, params)

        # prepare for the next step backward reconstruction
        copy_snapshot!(spt2, spt1)
        copy_wavefield!(wfd2, wfd1)
    end

    # close the boundary value file
    close(fid_bnd)

    return reshape(strength, params.nz, params.nx)
end

function reduce_source_strength!(strength::Vector{Tv}, wfd1::Wavefield,
                                 wfd2::Wavefield, spt1::Snapshot{Tv}, 
                                 spt2::Snapshot{Tv}, params::TdParams) where {Tv<:AbstractFloat}

    # total number of elements
    N = params.nz * params.nx

    for i = 1 : N

        j = params.spt2wfd[i]

        w2 = wfd2.p[i]         
        s2 = spt2.px[j]  + spt2.pz[j]                 
        tmp = w2 + s2 
        
        strength[i] += tmp * tmp         
    end
    
    return nothing
end


end # module

