module forward

using SeisAcoustic

function single_source_forward(src::Ts, params::TdParams;
                               rz="NULL", rx="NULL", location_flag="index", path_shot="NULL",
                               path_spt="NULL", path_wfd="NULL" , path_pre="NULL", interval=1,
                               path_bnd="NULL", path_lwfd="NULL") where {Ts<:Union{Source, Vector{Source}}}

    # allocate space for recordings
    if rz != "NULL" && rx != "NULL" && length(rz) == length(rx)
        rec = Recordings(rz, rx, params; location_flag=location_flag)
    end

    # initialize some variables
    spt1   = Snapshot(params)
    spt2   = Snapshot(params)
    tmp_z1 = zeros(params.data_format, params.Nz)
    tmp_z2 = zeros(params.data_format, params.Nz)
    tmp_x1 = zeros(params.data_format, params.Nx)
    tmp_x2 = zeros(params.data_format, params.Nx)

    # add source to the first snapshot
    add_source!(spt1, src, 1)

    # sample wavefield to get recordings
    rz != "NULL" && sample_spt2rec!(rec, spt1, 1)

    # save snapshot
    if path_spt != "NULL"
        hdr_spt   = snapshot_header(params, interval)
        fid_spt   = write_RSheader(path_spt, hdr_spt)
        append_one_snapshot(fid_spt, spt1)
    end

    # save wavefield
    if path_wfd != "NULL"
        hdr_wfd   = wavefield_header(params, interval)
        fid_wfd   = write_RSheader(path_wfd, hdr_wfd)
        append_one_wavefield(fid_wfd, spt1, params)
    end

    # save pressure
    if path_pre != "NULL"
        hdr_pre   = pressure_header(params, interval)
        fid_pre   = write_RSheader(path_pre, hdr_pre)
        append_one_pressure(fid_pre, spt1, params)
    end

    # boundary of wavefield
    if path_bnd != "NULL"
        hdr_bnd = boundary_header(params)
        fid_bnd = write_RSheader(path_bnd, hdr_bnd)
        append_one_boundary(fid_bnd, spt1, params)
    end

    # strength of source-side wavefield
    N = params.nz * params.nx
    strength = zeros(params.data_format, N)


    # applying time stepping nt-1 times
    for it = 2 : params.nt

        # one step forward
        one_step_forward!(spt2, spt1, params, tmp_z1, tmp_z2, tmp_x1, tmp_x2)

        # compute the strength of source-side wavefield
        forw_source_strength!(strength, spt1, spt2, params)

        # inject source term
        add_source!(spt2, src, it)

        # sample wavefield to get recordings
        rz != "NULL" && sample_spt2rec!(rec, spt2, it)

        # write current wavefield to disk
        if mod(it-1, interval) == 0
            path_spt != "NULL" && append_one_snapshot(fid_spt, spt2)
            path_wfd != "NULL" && append_one_wavefield(fid_wfd, spt2, params)
            path_pre != "NULL" && append_one_pressure(fid_pre, spt2, params)
        end

        # save boundary of wavefield
        path_bnd != "NULL" && append_one_boundary(fid_bnd, spt2, params)

        # prepare for next time stepping
        copy_snapshot!(spt1, spt2)
    end

    
    # save the last wavefield
    if path_lwfd != "NULL"
        wfd = sample_spt2wfd(spt1, params)
        write_wavefield(path_lwfd, wfd, params)
    end

    # close file
    path_spt != "NULL" && close(fid_spt)
    path_wfd != "NULL" && close(fid_wfd)
    path_pre != "NULL" && close(fid_pre)
    path_bnd != "NULL" && close(fid_bnd)

    if rz != "NULL"
        if path_shot == "NULL"
            return rec, reshape(strength, params.nz, params.nx)
        else
            write_recordings(path_shot, rec)
            return path_shot
        end
    else
        return nothing
    end
end

# define a auxillary function to improve efficiency
function forw_source_strength!(strength, spt1, spt2, params)

    # total number of elements
    N = params.nz * params.nx

    for i = 1 : N
        j = params.spt2wfd[i]
        p2= spt2.px[j] + spt2.pz[j]   

        strength[i] += p2 * p2
    end
    
    return nothing
end

end # module
