# download.jl
#
# Fetch a chain file to a local path. The download goes to a temporary file in
# the destination directory and is atomically renamed into place, so an
# interrupted transfer never leaves a half-written file where a valid one is
# expected. `Downloads.download` already raises on non-success HTTP status.

"""
    download_chain(url, path) -> path

Download `url` to `path` atomically, returning `path`.
"""
function download_chain(url::AbstractString, path::AbstractString)
    dir = dirname(abspath(path))
    mkpath(dir)
    tmp = tempname(dir)
    Downloads.download(url, tmp)
    mv(tmp, path; force = true)
    return path
end
