__precompile__()
module Utils

export py_sparse_matrix_to_julia
export load_qutip_object
export get_qutip_matrix

using StaticArrays
using SparseArrays
using PyCall

const qt = PyNULL()
const sp = PyNULL()

function __init__()
    copy!(qt, pyimport_conda("qutip", "qutip"))
    copy!(sp, pyimport_conda("scipy.sparse", "scipy"))
end

function py_sparse_matrix_to_julia(A; to_static=true)
    A = sp.csc_matrix(A)
    m, n = A.shape
    colptr = A.indptr .+ 1
    rowval = A.indices .+ 1
    nzval = A.data
    S = SparseMatrixCSC(m, n, colptr, rowval, nzval)
    if to_static
        return SMatrix{m, n}(S)
    else
        return S
    end
end

function load_qutip_object(path::String)
    Qobj = qt.fileio.qload(path)
    return Qobj
end

function get_qutip_matrix(path::String; to_static=true)
    Qobj = qt.fileio.qload(path)
    py_matrix = Qobj.data
    jl_matrix = py_sparse_matrix_to_julia(py_matrix;
                                          to_static=to_static)
    return jl_matrix
end

end
