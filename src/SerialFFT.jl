module SerialFFT

using LinearAlgebra
using FFTW

export fft_serinit, wvn_serinit, wvn_serinit!
export fft_for, fft_bac, fft_for!,fft_bac!
export getphase, phaseshift, phaseshift!
export getamp
export spec1d, spec2d, getcorr, specomni1d
export dealias,dealias!,filtering,filtering!, padding,padding!
export pdfx, pdfy, pdfxx, pdfyy

@deprecate fftserinit  fft_serinit
@deprecate wvnserinit  wvn_serinit
@deprecate wvnserinit! wvn_serinit!
@deprecate fftfor     fft_for
@deprecate fftfor!    fft_for!
@deprecate fftbac     fft_bac
@deprecate fftbac!    fft_bac!
@deprecate fftfor1d   fft_for
@deprecate fftbac1d   fft_bac
@deprecate fftfor1d!  fft_for!
@deprecate fftbac1d!  fft_bac!
@deprecate fftfor2d   fft_for
@deprecate fftbac2d   fft_bac
@deprecate fftfor2d!  fft_for!
@deprecate fftbac2d!  fft_bac!
@deprecate fftforx2d  fft_for
@deprecate fftbacx2d  fft_bac
@deprecate fftforx2d! fft_for!
@deprecate fftbacx2d! fft_bac!
@deprecate fftfory2d  fft_for
@deprecate fftbacy2d  fft_bac
@deprecate fftfory2d! fft_for!
@deprecate fftbacy2d! fft_bac!

"""
    wvn_serinit(nx::Integer,ny::Integer)
Calculate the 2D wavenumber corresponding to fft functions in the package SerialFFT and writes to two 1D arrays.
"""
function wvn_serinit(nx::Integer,ny::Integer)
        wvnx=zeros(nx)
        wvny=zeros(ny)
        for i=1:div(nx,2)
                    wvnx[i] = (i-1)
                end
        wvnx[div(nx,2)+1] = div(nx,2)
        for i=div(nx,2)+2:nx
                    wvnx[i] = (nx+1-i)
                end

        for j=1:div(ny,2)
                    wvny[j] = (j-1)
                end
        wvny[div(ny,2)+1] = div(ny,2)
        for j=div(ny,2)+2:ny
                    wvny[j] = (ny+1-j)
                end

        return wvnx, wvny
    end


"""
    wvn_serinit(arraysize::Tuple{<:Integer,<:Integer})
Use a tuple as input.
"""
function wvn_serinit(arraysize::Tuple{<:Integer,<:Integer})
        nx=arraysize[1]
        ny=arraysize[2]
        wvnx=zeros(nx)
        wvny=zeros(ny)
        for i=1:div(nx,2)
                    wvnx[i] = (i-1)
                end
        wvnx[div(nx,2)+1] = div(nx,2)
        for i=div(nx,2)+2:nx
                    wvnx[i] = (nx+1-i)
                end

        for j=1:div(ny,2)
                    wvny[j] = (j-1)
                end
        wvny[div(ny,2)+1] = div(ny,2)
        for j=div(ny,2)+2:ny
                    wvny[j] = (ny+1-j)
                end

        return wvnx, wvny
    end

"""
    wvn_serinit(nx::Integer)
Calculate the 1D wavenumber corresponding to fft functions in the package SerialFFT.
"""
function wvn_serinit(nx::Integer)
        wvnx=zeros(nx)
        for i=1:div(nx,2)
                    wvnx[i] = (i-1)
                end
        wvnx[div(nx,2)+1] = div(nx,2)
        for i=div(nx,2)+2:nx
                    wvnx[i] = (nx+1-i)
                end
        return wvnx
    end

"""
    wvn_serinit!(wvnx::AbstractArray,wvny::AbstractArray)
Calculate the 2D wavenumber corresponding to fft functions in the package SerialFFT and write the results to 1D arrays `wvnx` and `wvny`.
"""
function wvn_serinit!(wvnx::AbstractArray,wvny::AbstractArray)
        nx=size(wvnx,1)
        ny=size(wvny,1)
        for i=1:div(nx,2)
                    wvnx[i] = (i-1)
                end
        wvnx[div(nx,2)+1] = div(nx,2)
        for i=div(nx,2)+2:nx
                    wvnx[i] = (nx+1-i)
                end

        for j=1:div(ny,2)
                    wvny[j] = (j-1)
                end
        wvny[div(ny,2)+1] = div(ny,2)
        for j=div(ny,2)+2:ny
                    wvny[j] = (ny+1-j)
                end

        return nothing
    end


"""
    wvn_serinit!(wvnx::AbstractArray{<:Real})
Calculate the 1D wavenumber corresponding to fft functions in the package SerialFFT and write the results to 1D array `wvnx`.
"""
function wvn_serinit!(wvnx::AbstractArray{<:Real})
        nx=size(wvnx,1)
        for i=1:div(nx,2)
                    wvnx[i] = (i-1)
                end
        wvnx[div(nx,2)+1] = div(nx,2)
        for i=div(nx,2)+2:nx
                    wvnx[i] = (nx+1-i)
                end
        return nothing
    end


"""
    getamp(a::AbstractArray,k::Integer)::Real

Calculate the amplitude `A` of a signal `a(x)` at the wavenumber `k` assuming `a=Acos(kx+θ)`.
"""
function getamp(a::AbstractArray,k::Integer)::Real

        n=size(a,1)
        fa=similar(a)
        plfor = FFTW.plan_r2r(a,FFTW.R2HC,1)
        mul!(fa,plfor,a)
        fa./=n
        return 2*sqrt(fa[n-k+1]^2+fa[k+1]^2)
    end


"""
    getphase(a::AbstractArray,k::Integer)::Real

Calculate the phase angle `θ` (rad) of a signal `a(x)` at the wavenumber `k` assuming `a=Acos(kx+θ)`.
"""
function getphase(a::AbstractArray,k::Integer)::Real

        n=size(a,1)
        fa=similar(a)
        plfor = FFTW.plan_r2r(a,FFTW.R2HC,1)
        mul!(fa,plfor,a)
        return atan(fa[n-k+1],fa[k+1])
    end

"""
    phaseshift(a::AbstractMatrix,θ::Real)::AbstractMatrix

Calculate the phase-shifted signal `fa(x)=a(x-θ)`, assuming `0<x<2π`, where x is the first dimension of `a`.
"""
function phaseshift(a::AbstractArray,θ::Real)::AbstractArray

        fa=similar(a)
        phaseshift!(fa,a,θ)
        return fa
    end


"""
    phaseshift!(fa::AbstractArray{<:Real,3},a::AbstractArray{<:Real,3},θ::Real)

Calculate the 3D phase-shifted signal `fa(x,y,z)=a(x-θ,y,z)` and write results to `fa`, assuming `0<x<2π`.
"""
function phaseshift!(fa::AbstractArray{<:Real,3},a::AbstractArray{<:Real,3},θ::Real)

        fa2d=similar(a[:,:,1])
        for k=1:size(a,3)
                    phaseshift!(fa2d,a[:,:,k],θ)
                    fa[:,:,k] = fa2d
                end
        return nothing
    end


# """
#     phaseshift(a::AbstractMatrix,theta::Real)::AbstractMatrix

# Calculate the 2D phase-shifted signal `fa(x,y)=a(x-θ,y)`, assuming `0<x<2π`.
# """
# function phaseshift(a::AbstractMatrix,theta::Real)::AbstractMatrix

#     n =size(a,1)
#     fa=similar(a)
#     bufa=similar(a)

#     plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
#     plbacx = FFTW.plan_r2r(similar(a),FFTW.HC2R,1)
#     mul!(fa,plforx,a)
#     fa./=n

#     bufa[1,:] = fa[1,:]
#     for i=2:div(n,2)
#         bufa[i,:]=fa[i,:]*cos((i-1)*theta)+fa[n+2-i,:]*sin((i-1)*theta)
#     end
#     bufa[div(n,2)+1,:] = fa[div(n,2)+1,:]*cos((n-div(n,2))*theta)
#     for i=div(n,2)+2:n
#         bufa[i,:]=fa[i,:]*cos((n+1-i)*theta)-fa[n+2-i,:]*sin((n+1-i)*theta)
#     end

#     mul!(fa,plbacx,bufa)
#     return fa
# end


"""
    phaseshift!(fa::AbstractMatrix,a::AbstractMatrix,θ::Real)

Calculate the 2D phase-shifted signal `fa(x,y)=a(x-θ,y)` and write results to `fa`, assuming `0<x<2π`.
"""
function phaseshift!(fa::AbstractMatrix,a::AbstractMatrix,θ::Real)

        fa1d=similar(a[:,1])
        for k=1:size(a,2)
                    phaseshift!(fa1d,a[:,k],θ)
                    fa[:,k] = fa1d
                end
        return nothing

    end


# """
#     phaseshift(a::AbstractMatrix,theta::Real)::AbstractMatrix

# Calculate the 1D phase-shifted signal `fa(x)=a(x-θ)`, assuming `0<x<2π`.
# """
# function phaseshift(a::AbstractArray,theta::Real)::AbstractArray

#     n =size(a,1)
#     fa=similar(a)
#     bufa=similar(a)

#     plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
#     plbacx = FFTW.plan_r2r(similar(a),FFTW.HC2R,1)
#     mul!(fa,plforx,a)
#     fa./=n

#     bufa[1]=fa[1]
#     for i=2:div(n,2)
#         bufa[i]=fa[i]*cos((i-1)*theta)+fa[n+2-i]*sin((i-1)*theta)
#     end
#     bufa[div(n,2)+1]=fa[div(n,2)+1]*cos((n-div(n,2))*theta)
#     for i=div(n,2)+2:n
#         bufa[i]=fa[i]*cos((n+1-i)*theta)-fa[n+2-i]*sin((n+1-i)*theta)
#     end

#     mul!(fa,plbacx,bufa)
#     return fa
# end


"""
    phaseshift!(fa::AbstractArray,a::AbstractArray,θ::Real)

Calculate the 1D phase-shifted signal `fa(x)=a(x-θ)` and write results to `fa`, assuming `0<x<2π`.
"""
function phaseshift!(fa::AbstractArray,a::AbstractArray,θ::Real)

        n =size(a,1)
        bufa=similar(a)

        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        plbacx = FFTW.plan_r2r(similar(a),FFTW.HC2R,1)
        mul!(fa,plforx,a)
        fa./=n

        bufa[1]=fa[1]
        for i=2:div(n,2)
                    bufa[i]=fa[i]*cos((i-1)*θ)+fa[n+2-i]*sin((i-1)*θ)
                end
        bufa[div(n,2)+1]=fa[div(n,2)+1]*cos((n-div(n,2))*θ)
        for i=div(n,2)+2:n
                    bufa[i]=fa[i]*cos((n+1-i)*θ)-fa[n+2-i]*sin((n+1-i)*θ)
                end

        mul!(fa,plbacx,bufa)
        return nothing
    end


function fft_serinit(a::AbstractArray{<:Real,2})
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        plbacx = FFTW.plan_r2r(similar(a),FFTW.HC2R,1)
        plfory = FFTW.plan_r2r(similar(a),FFTW.R2HC,2)
        plbacy = FFTW.plan_r2r(similar(a),FFTW.HC2R,2)
        bufa=similar(a)
        return plforx, plbacx, plfory, plbacy, bufa
    end

function fft_serinit(a::AbstractArray{<:Real,1})
        plfor = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        plbac = FFTW.plan_r2r(similar(a),FFTW.HC2R,1)
        return plfor, plbac
    end

"""
    fft_for(a::AbstractArray{<:Real,1})::AbstractArray

Calculate the Fourier transform of a 1D array.
"""
function fft_for(a::AbstractArray{<:Real,1})::AbstractArray
        fa=similar(a)
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        mul!(fa,plforx,a)
        fa./=size(a,1)
        return fa
    end


"""
    fft_for!(a::AbstractArray{<:Real,1})

Calculate the Fourier transform of a 1D array and overwrite the input.
"""
function fft_for!(a::AbstractArray{<:Real,1})
        fa=similar(a)
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        mul!(fa,plforx,a)
        fa./=size(a,1)
        a[:]=fa
        return nothing
    end

"""
    fft_for!(fa::AbstractArray{<:Real,1},a::AbstractArray{<:Real,1})

Calculate the Fourier transform of a 1D array and write to `fa`.
"""
function fft_for!(fa::AbstractArray{<:Real,1},a::AbstractArray{<:Real,1})
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        mul!(fa,plforx,a)
        fa./=size(a,1)
        return nothing
    end

"""
    fft_for(a::AbstractArray{<:Real,2})::AbstractArray

Calculate the 2D Fourier transform of a 2D array.
"""
function fft_for(a::AbstractArray{<:Real,2})::AbstractArray
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        plfory = FFTW.plan_r2r(similar(a),FFTW.R2HC,2)
        bufa=similar(a)
        fa=similar(a)
        mul!(bufa,plforx,a)
        mul!(fa,plfory,bufa)
        fa./=size(a,1)*size(a,2)
        return fa
    end

"""
    fft_for(a::AbstractArray{<:Real,2},plforx,plfory,bufa)::AbstractArray

Calculate the 2D Fourier transform of a 2D array using external plan and buffer array.
"""
function fft_for(a::AbstractArray{<:Real,2},plforx,plfory,bufa)::AbstractArray
        fa=similar(a)
        mul!(bufa,plforx,a)
        mul!(fa,plfory,bufa)
        fa./=size(a,1)*size(a,2)
        return fa
    end

"""
    fft_for!(a::AbstractArray{<:Real,2})

Calculate the 2D Fourier transform of a 2D array and overwrite the input.
"""
function fft_for!(a::AbstractArray{<:Real,2})
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        plfory = FFTW.plan_r2r(similar(a),FFTW.R2HC,2)
        bufa=similar(a)
        mul!(bufa,plforx,a)
        mul!(a,plfory,bufa)
        a./=size(a,1)*size(a,2)
        return nothing
    end

"""
    fft_for!(a::AbstractArray{<:Real,2},plforx,plfory,bufa)

Calculate the 2D Fourier transform of a 2D array using external plan and buffer array, and overwrite the input.
"""
function fft_for!(a::AbstractArray{<:Real,2},plforx,plfory,bufa)
        mul!(bufa,plforx,a)
        mul!(a,plfory,bufa)
        a./=size(a,1)*size(a,2)
        return nothing
    end

"""
    fft_for!(fa::AbstractArray{<:Real,2},a::AbstractArray{<:Real,2})

Calculate the 2D Fourier transform of a 2D array and write to `fa`.
"""
function fft_for!(fa::AbstractArray{<:Real,2},a::AbstractArray{<:Real,2})
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        plfory = FFTW.plan_r2r(similar(a),FFTW.R2HC,2)
        bufa=similar(a)
        mul!(bufa,plforx,a)
        mul!(fa,plfory,bufa)
        fa./=size(a,1)*size(a,2)
        return nothing
    end

"""
    fft_for!(fa::AbstractArray{<:Real,2},a::AbstractArray{<:Real,2},plforx,plfory,bufa)

Calculate the 2D Fourier transform of a 2D array using external plan and buffer array, and write to `fa`.
"""
function fft_for!(fa::AbstractArray{<:Real,2},a::AbstractArray{<:Real,2},plforx,plfory,bufa)
        mul!(bufa,plforx,a)
        mul!(fa,plfory,bufa)
        fa./=size(a,1)*size(a,2)
        return nothing
    end

"""
    fft_for(a::AbstractArray{<:Real,2},axis::Integer)::AbstractArray{<:Real,2}

Calculate the 1D Fourier transform of a 2D array in the first (`axis=1`) or the second (`axis=2`) dimension.
"""
function fft_for(a::AbstractArray{<:Real,2},axis::Integer)::AbstractArray{<:Real,2}
        if axis!=1 && axis!=2
                    ArgumentError("Axis must be either 1 or 2!")
                end
        fa=similar(a)
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,axis)
        mul!(fa,plforx,a)
        fa./=size(a,axis)
        return fa
    end

"""
    fft_for(a::AbstractArray{<:Real,2},plfor)::AbstractArray{<:Real,2}

Calculate the 1D Fourier transform of a 2D array in the first dimension with external plan.
"""
function fft_for(a::AbstractArray{<:Real,2},plfor)::AbstractArray{<:Real,2}
        fa=similar(a)
        mul!(fa,plfor,a)
        fa./=size(a,plfor.region)
        return fa
    end

"""
    fft_for!(a::AbstractArray{<:Real,2},plfor)

Calculate the 1D Fourier transform of a 2D array in the first dimension with external plan, and overwrite the input.
"""
function fft_for!(a::AbstractArray{<:Real,2},plfor)
        bufa=similar(a)
        mul!(bufa,plfor,a)
        a .= bufa./size(a,plfor.region)
        return nothing
    end

"""
    fft_for!(a::AbstractArray{<:Real,2},axis::Integer)

Calculate the 1D Fourier transform of a 2D array in the first (`axis=1`) or the second (`axis=2`) dimension and overwrite the input.
"""
function fft_for!(a::AbstractArray{<:Real,2},axis::Integer)
        if axis!=1 && axis!=2
                    ArgumentError("Axis must be either 1 or 2!")
                end
        bufa=similar(a)
        plfor = FFTW.plan_r2r(similar(a),FFTW.R2HC,axis)
        mul!(bufa,plfor,a)
        a .= bufa ./ size(a,axis)
        return nothing
    end

"""
    fft_for!(fa::AbstractArray{<:Real,2},a::AbstractArray{<:Real,2},plfor)

Calculate the 1D Fourier transform of a 2D array with external plan and write to `fa`.
"""
function fft_for!(fa::AbstractArray{<:Real,2},a::AbstractArray{<:Real,2},plfor)
        mul!(fa,plfor,a)
        fa./=size(a,plfor.region)
        return nothing
    end

"""
    fft_for!(fa::AbstractArray{<:Real,2},a::AbstractArray{<:Real,2},axis::Integer)

Calculate the 1D Fourier transform of a 2D array in the first (`axis=1`) or the second (`axis=2`) dimension and write to `fa`.
"""
function fft_for!(fa::AbstractArray{<:Real,2},a::AbstractArray{<:Real,2},axis::Integer)
        if axis!=1 && axis!=2
                    ArgumentError("Axis must be either 1 or 2!")
                end
        plfor = FFTW.plan_r2r(similar(a),FFTW.R2HC,axis)
        mul!(fa,plfor,a)
        fa./=size(a,axis)
        return nothing
    end

"""
    fft_bac(fa::AbstractArray{<:Real,1})::AbstractArray

Calculate the inverse Fourier transform of a 1D array.
"""
function fft_bac(fa::AbstractArray{<:Real,1})::AbstractArray
        a=similar(fa)
        plbacx = FFTW.plan_r2r(similar(a),FFTW.HC2R,1)
        mul!(a,plbacx,fa)
        return a
    end

"""
    fft_bac!(fa::AbstractArray{<:Real,1})

Calculate the inverse Fourier transform of a 1D array and overwrite the input.
"""
function fft_bac!(fa::AbstractArray{<:Real,1})
        a=similar(fa)
        plbacx = FFTW.plan_r2r(similar(a),FFTW.HC2R,1)
        mul!(a,plbacx,fa)
        fa[:]=a
        return nothing
    end

"""
    fft_bac!(a::AbstractArray{<:Real,1},fa::AbstractArray{<:Real,1})

Calculate the inverse Fourier transform of a 1D array and write to `fa`.
"""
function fft_bac!(a::AbstractArray{<:Real,1},fa::AbstractArray{<:Real,1})
        plbacx = FFTW.plan_r2r(similar(a),FFTW.HC2R,1)
        mul!(a,plbacx,fa)
        return nothing
    end

"""
    fft_bac(fa::AbstractArray{<:Real,2})::AbstractArray

Calculate the 2D inverse Fourier transform of a 2D array.
"""
function fft_bac(fa::AbstractArray{<:Real,2})::AbstractArray{<:Real,2}
        plbacx = FFTW.plan_r2r(similar(fa),FFTW.HC2R,1)
        plbacy = FFTW.plan_r2r(similar(fa),FFTW.HC2R,2)
        bufa=similar(fa)
        a=similar(fa)
        mul!(bufa,plbacy,fa)
        mul!(a,plbacx,bufa)
        return a
    end

"""
    fft_bac(fa::AbstractArray{<:Real,2},plforx,plfory,bufa)::AbstractArray

Calculate the 2D inverse Fourier transform of a 2D array using external plan and buffer array.
"""
function fft_bac(fa::AbstractArray{<:Real,2},plbacx,plbacy,bufa)::AbstractArray{<:Real,2}
        plbacx = FFTW.plan_r2r(similar(fa),FFTW.HC2R,1)
        plbacy = FFTW.plan_r2r(similar(fa),FFTW.HC2R,2)
        bufa=similar(fa)
        a=similar(fa)
        mul!(bufa,plbacy,fa)
        mul!(a,plbacx,bufa)
        return a
    end

"""
    fft_bac!(fa::AbstractArray{<:Real,2})

Calculate the 2D inverse Fourier transform of a 2D array and overwrite the input.
"""
function fft_bac!(fa::AbstractArray{<:Real,2})
        plbacx = FFTW.plan_r2r(similar(fa),FFTW.HC2R,1)
        plbacy = FFTW.plan_r2r(similar(fa),FFTW.HC2R,2)
        bufa=similar(fa)
        mul!(bufa,plbacy,fa)
        mul!(fa,plbacx,bufa)
        return nothing
    end

"""
    fft_bac!(fa::AbstractArray{<:Real,2},plforx,plfory,bufa)

Calculate the 2D inverse Fourier transform of a 2D array using external plan and buffer array, and overwrite the input.
"""
function fft_bac!(fa::AbstractArray{<:Real,2},plbacx,plbacy,bufa)
        mul!(bufa,plbacy,fa)
        mul!(fa,plbacx,bufa)
        return nothing
    end

"""
    fft_bac!(a::AbstractArray{<:Real,2},fa::AbstractArray{<:Real,2})

Calculate the 2D Fourier transform of a 2D array and write to `fa`.
"""
function fft_bac!(a::AbstractArray{<:Real,2},fa::AbstractArray{<:Real,2})
        plbacx = FFTW.plan_r2r(similar(fa),FFTW.HC2R,1)
        plbacy = FFTW.plan_r2r(similar(fa),FFTW.HC2R,2)
        bufa=similar(fa)
        mul!(bufa,plbacy,fa)
        mul!(a,plbacx,bufa)
        return nothing
    end

"""
    fft_bac!(a::AbstractArray{<:Real,2},fa::AbstractArray{<:Real,2},plbacx,plbacy,bufa)

Calculate the 2D inverse Fourier transform of a 2D array using external plan and buffer array, and write to `a`.
"""
function fft_bac!(a::AbstractArray{<:Real,2},fa::AbstractArray{<:Real,2},plbacx,plbacy,bufa)
        mul!(bufa,plbacy,fa)
        mul!(a,plbacx,bufa)
        return nothing
    end

"""
    fft_bac(fa::AbstractArray{<:Real,2},axis::Integer)::AbstractArray

Calculate the 1D inverse Fourier transform of a 2D array in the first (`axis=1`) or the second (`axis=2`) dimension.
"""
function fft_bac(fa::AbstractArray{<:Real,2},axis::Integer)::AbstractArray
        if axis!=1 && axis!=2
                    ArgumentError("Axis must be either 1 or 2!")
                end
        a=similar(fa)
        plbacx = FFTW.plan_r2r(similar(a),FFTW.HC2R,axis)
        mul!(a,plbacx,fa)
        return a
    end

"""
    fft_bac(fa::AbstractArray{<:Real,2},plbac)::AbstractArray

Calculate the 1D Fourier transform of a 2D array in the first dimension with external plan.
"""
function fft_bac(fa::AbstractArray{<:Real,2},plbac)::AbstractArray
        a=similar(fa)
        mul!(a,plbac,fa)
        return a
    end

"""
    fft_bac!(fa::AbstractArray{<:Real,2},plbac)

Calculate the 1D inverse Fourier transform of a 2D array in the first dimension with external plan, and overwrite the input.
"""
function fft_bac!(fa::AbstractArray{<:Real,2},plbac)
        bufa=similar(fa)
        mul!(bufa,plbac,fa)
        fa .= bufa
        return nothing
    end

"""
    fft_bac!(fa::AbstractArray{<:Real,2},axis::Integer)

Calculate the 1D inverse Fourier transform of a 2D array in the first (`axis=1`) or the second (`axis=2`) dimension and overwrite the input.
"""
function fft_bac!(fa::AbstractArray{<:Real,2},axis::Integer)
        if axis!=1 && axis!=2
                    ArgumentError("Axis must be either 1 or 2!")
                end
        bufa=similar(fa)
        plbac = FFTW.plan_r2r(similar(fa),FFTW.HC2R,axis)
        mul!(bufa,plbac,fa)
        fa .= bufa
        return nothing
    end

"""
    fft_bac!(a::AbstractArray{<:Real,2},fa::AbstractArray{<:Real,2},plbac)

Calculate the 1D Fourier transform of a 2D array with external plan and write to `fa`.
"""
function fft_bac!(a::AbstractArray{<:Real,2},fa::AbstractArray{<:Real,2},plbac)
        mul!(a,plbac,fa)
        return nothing
    end

"""
    fft_bac!(a::AbstractArray{<:Real,2},fa::AbstractArray{<:Real,2},axis::Integer)

Calculate the 1D inverse Fourier transform of a 2D array in the first (`axis=1`) or the second (`axis=2`) dimension and write to `a`.
"""
function fft_bac!(a::AbstractArray{<:Real,2},fa::AbstractArray{<:Real,2},axis::Integer)
        if axis!=1 && axis!=2
                    ArgumentError("Axis must be either 1 or 2!")
                end
        plbac = FFTW.plan_r2r(similar(a),FFTW.HC2R,axis)
        mul!(a,plbac,fa)
        return nothing
    end


"""
    spec1d(a::AbstractArray{<:Real,1},dk=1::Real)

Energy density spectrum of 1D array `a` with optional wavenumber step `dk`.
"""
function spec1d(a::AbstractArray{<:Real,1},dk=1::Real)
        n=size(a,1)
        fa=similar(a)
        sk=similar(a[1:div(n,2)])
        plfor = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        mul!(fa,plfor,a)
        fa./=n

        k= [i*dk for i=1:div(n,2)]
        for i=1:div(n,2)
                    ftn = fa[i+1]^2 + fa[n+1-i]^2
                    sk[i] = 2*ftn/dk
                end
        return k,sk
    end

"""
    spec1d(a::AbstractArray{<:Real,2},axis::Integer,dk=1::Real)

When the input `a` is a 2D array, the spectra is calculated along the first (`axis=1`) or the second (`axis=2`) dimension.
"""
function spec1d(a::AbstractArray{<:Real,2},axis::Integer,dk=1::Real)
        if axis!=1 && axis!=2
                    ArgumentError("Axis must be either 1 or 2!")
                end

        n=size(a,1)
        m=size(a,2)
        nf=size(a,axis)
        fa=similar(a)
        sk=zeros(div(nf,2))
        plfor = FFTW.plan_r2r(similar(a),FFTW.R2HC,axis)
        mul!(fa,plfor,a)
        fa./=nf

        k= [i*dk for i=1:div(nf,2)]
        if axis==1
                    for j=1:m, i=1:div(nf,2)
                                    ftn = fa[i+1,j]^2 + fa[n+1-i,j]^2
                                    sk[i] = sk[i] + 2*ftn/dk/m
                                end
                end
        if axis==2
                    for j=1:n, i=1:div(nf,2)
                                    ftn = fa[j,i+1]^2 + fa[j,m+1-i]^2
                                    sk[i] = sk[i] + 2*ftn/dk/n
                                end
                end
        return k,sk
    end

"""
    spec2d(a::AbstractArray{<:Real,2},dk=(1,1)::Tuple{<:Real,<:Real})

2D energy density spectra of a 2D array `a` with `dk` being the optional wavenumber tuple (defaults to `dk=(1,1)`).
"""
function spec2d(a::AbstractArray{<:Real,2},dk=(1,1)::Tuple{<:Real,<:Real})

        n=size(a,1)
        m=size(a,2)
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        plfory = FFTW.plan_r2r(similar(a),FFTW.R2HC,2)
        bufa=similar(a)
        fa=similar(a)
        mul!(bufa,plforx,a)
        mul!(fa,plfory,bufa)
        fa./=size(a,1)*size(a,2)

        kx=zeros(div(n,2),m+1)
        ky=zeros(div(n,2),m+1)
        for i=1:div(n,2)
                    kx[i,:] .= dk[1] * i
                end
        for j=1:m+1
                    ky[:,j] .= dk[2] * (j - div(m,2) - 1)
                end

        skxy=zeros(div(n,2),m+1)

        for j=0:div(m,2)-1, i=1:div(n,2)
                    if j==0
                                    skxy[i,div(m,2)+1] = (fa[i+1,1]^2+fa[n+1-i,1]^2) / dk[1] / dk[2]
                                else
                                    a1,a2,alpha,beta=spec2phy(fa[i+1,j+1],fa[n+1-i,j+1],fa[i+1,m+1-j],fa[n+1-i,m+1-j])
                                    skxy[i,div(m,2)+1+j] = a1^2 / 2 / dk[1] / dk[2]
                                    skxy[i,div(m,2)+1-j] = a2^2 / 2 / dk[1] / dk[2]
                                end
                end

        return kx, ky, skxy
    end

"""
    specomni1d(a::AbstractArray{<:Real,2},nk::Integer,dk::Real,dkraw=(1,1)::Tuple{<:Real,<:Real})

Omnidirectional 1D energy density spectra of a 2D array `a` with `nk` being the number of spectral bins and `dk` being the wavenumber resolution.
"""
function specomni1d(a::AbstractArray{<:Real,2},nk::Integer,dk::Real,dkraw=(1,1)::Tuple{<:Real,<:Real})

        n=size(a,1)
        m=size(a,2)
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        plfory = FFTW.plan_r2r(similar(a),FFTW.R2HC,2)
        bufa=similar(a)
        fa=similar(a)
        mul!(bufa,plforx,a)
        mul!(fa,plfory,bufa)
        fa./=size(a,1)*size(a,2)

        kx=zeros(div(n,2),m+1)
        ky=zeros(div(n,2),m+1)
        for i=1:div(n,2)
                    kx[i,:] .= dkraw[1] * i
                end
        for j=1:m+1
                    ky[:,j] .= dkraw[2] * (j - div(m,2) - 1)
                end

        sk=zeros(nk)
        allk=[i*dk for i=1:nk]

        for j=0:div(m,2)-1, i=1:div(n,2)
                    wav=sqrt((kx[i,div(m,2)-j])^2+(ky[i,div(m,2)-j])^2)
                    for k=1:nk
                                    if wav>=(k-0.5)*dk && wav<(k+0.5)*dk
                                                        if j==0
                                                                                sk[k]+=4*(fa[i+1,1]^2+fa[n+1-i,1]^2)/dk
                                                                            else
                                                                                sk[k]+=4*(fa[i+1,j+1]^2+fa[n+1-i,j+1]^2+fa[i+1,m+1-j]^2+fa[n+1-i,m+1-j]^2)/dk
                                                                            end
                                                    end
                                end
                end

        return allk, sk
    end

"""
    function getcorr(a::AbstractArray{<:Real,2},axis::Integer,pexy::Real,n::Integer,dr::Real)
Calculate the 2-point correlations of a 2D array in dimension `axis`. The base wavenumber is given by
`pexy`, `n` is the length of the distance and `dr` is the distance interval. Returns `r` (distance)
and `R` correlation.
"""
function getcorr(a::AbstractArray{<:Real,2},axis::Integer,pexy::Real,n::Integer,dr::Real)

        R=zeros(n)
        allk,sk=spec1d(a,axis,pexy)
        nk=size(allk,1)
        dk=allk[2] - allk[1]

        r=[(i-1)*dr for i=1:n]
        for i=1:n
                    for k=1:nk
                                    R[i]=R[i]+sk[k]*cos(allk[k]*r[i])*dk
                                end
                end
        return r,R./R[1]
    end


"""
    pdfx(a::AbstractArray{<:Real,1},pex::Real)
Calculate the first derivative of a 1D array.
"""
function pdfx(a::AbstractArray{<:Real,1},pex::Real)
        fa=similar(a)
        nx=size(a,1)
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        plbacx = FFTW.plan_r2r(similar(a),FFTW.HC2R,1)

        bufx=similar(a)
        mul!(fa,plforx,a)
        fa./=nx
        bufx[1]=0
        bufx[div(nx,2)+1] = 0
        for i=2:div(nx,2)
                    wvn=(i-1)*pex
                    bufx[nx+2-i] = fa[i] * wvn
                    bufx[i] = -fa[nx+2-i] * wvn
                end
        mul!(fa,plbacx,bufx)
        return fa
    end

"""
    pdfx(a::AbstractArray{<:Real,2},pex::Real)::AbstractArray{<:Real,2}
Calculate the first derivative of a 2D array in the first dimension.
"""
function pdfx(a::AbstractArray{<:Real,2},pex::Real)::AbstractArray{<:Real,2}
        fa=similar(a)
        nx=size(a,1)
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        plbacx = FFTW.plan_r2r(similar(a),FFTW.HC2R,1)

        bufx=similar(a)
        mul!(fa,plforx,a)
        fa./=nx
        bufx[1,:] .=0
        bufx[div(nx,2)+1,:] .= 0
        for i=2:div(nx,2)
                    wvn=(i-1)*pex
                    bufx[nx+2-i,:] = fa[i,:] .* wvn
                    bufx[i,:] = -fa[nx+2-i,:] .* wvn
                end
        mul!(fa,plbacx,bufx)
        return fa
    end

"""
    pdfx(a::AbstractArray{<:Real,3},pex::Real)::AbstractArray{<:Real,3}
Calculate the first derivative of a 3D array in the first dimension.
"""
function pdfx(a::AbstractArray{<:Real,3},pex::Real)::AbstractArray{<:Real,3}
        fa=similar(a)
        nz=size(a,3)
        for k=1:nz
                    fa[:,:,k] = pdfx(a[:,:,k],pex)
                end
        return fa
    end

"""
    pdfy(a::AbstractArray{<:Real,2},pey::Real)::AbstractArray{<:Real,2}
Calculate the first derivative of a 2D array in the second dimension.
"""
function pdfy(a::AbstractArray{<:Real,2},pey::Real)::AbstractArray{<:Real,2}
        return Array(transpose(pdfx(Array(transpose(a)),pey)))
    end


"""
    pdfy(a::AbstractArray{<:Real,3},pey::Real)::AbstractArray{<:Real,3}
Calculate the first derivative of a 3D array in the second dimension.
"""
function pdfy(a::AbstractArray{<:Real,3},pey::Real)::AbstractArray{<:Real,3}
        fa=similar(a)
        nz=size(a,3)
        for k=1:nz
                    fa[:,:,k] = pdfy(a[:,:,k],pey)
                end
        return fa
    end

"""
    pdfxx(a::AbstractArray{<:Real,2},pex::Real)::AbstractArray{<:Real,2}
Calculate the second derivative of a 2D array in the first dimension.
"""
function pdfxx(a::AbstractArray{<:Real,2},pex::Real)::AbstractArray{<:Real,2}
        fa=similar(a)
        nx=size(a,1)
        plforx = FFTW.plan_r2r(similar(a),FFTW.R2HC,1)
        plbacx = FFTW.plan_r2r(similar(a),FFTW.HC2R,1)

        bufx=similar(a)
        mul!(fa,plforx,a)
        fa./=nx
        bufx[1,:] .=0
        bufx[div(nx,2)+1,:] .= 0
        for i=2:div(nx,2)
                    wvn=(i-1)*pex
                    bufx[nx+2-i,:] = -fa[nx+2-i,:] .* wvn^2
                    bufx[i,:] = -fa[i,:] .* wvn^2
                end
        mul!(fa,plbacx,bufx)
        return fa
    end

"""
    pdfyy(a::AbstractArray{<:Real,2},pey::Real)::AbstractArray{<:Real,2}
Calculate the second derivative of a 2D array in the second dimension.
"""
function pdfyy(a::AbstractArray{<:Real,2},pey::Real)::AbstractArray{<:Real,2}
        return Array(transpose(pdfxx(Array(transpose(a)),pey)))
    end

"""
    dealias(ah::AbstractArray)

Dealiasing of the first and second dimension of array `ah` by removing the top 2/3 wavenumbers.

"""
function dealias(ah::AbstractArray)
        al=copy(ah)
        dealias!(al)
        return al
    end

"""
    dealias!(ah::AbstractArray{<:Real,3})

Dealiasing of a 3D array `ah` by removing the top 2/3 wavenumbers and write to `ah`.

"""
function dealias!(ah::AbstractArray{<:Real,3})
        nh1=size(ah,1)
        nh2=size(ah,2)

        wvnx, wvny=wvn_serinit(nh1,nh2)
        wvnxm=div(nh1,2)
        wvnym=div(nh2,2)
        for k=1:size(ah,3)
                    al = fft_for(ah[:,:,k],1)
                    for i=1:nh1
                                    if wvnx[i]>=div(2*wvnxm,3)
                                                        al[i,:] .= 0
                                                    end
                                end
                    fft_for!(al,2)
                    for j=1:nh2
                                    if wvny[j]>=div(2*wvnym,3)
                                                        al[:,j] .= 0
                                                    end
                                end
                    ah[:,:,k]=fft_bac(al)
                end
        return nothing
    end

"""
    dealias!(ah::AbstractArray{<:Real,3},dims::Integer)

Dealiasing of a 3D array `ah` by removing the top 2/3 wavenumbers of
dimension `dims` and write to `ah`.

"""
function dealias!(ah::AbstractArray{<:Real,3},dims::Integer)

        if dims!=1 && dims!=2
                    error("Dims must be either 1 or 2!")
                end
        nh1=size(ah,1)
        nh2=size(ah,2)
        wvnx, wvny=wvn_serinit(nh1,nh2)
        wvnxm=div(nh1,2)
        wvnym=div(nh2,2)
        if dims==1
                    for k=1:size(ah,3)
                                    al = fft_for(ah[:,:,k],1)
                                    for i=1:nh1
                                                        if wvnx[i]>=div(2*wvnxm,3)
                                                                                al[i,:] .= 0
                                                                            end
                                                    end
                                    ah[:,:,k]=fft_bac(al,1)
                                end
                else
                    for k=1:size(ah,3)
                                    al = fft_for(ah[:,:,k],2)
                                    for j=1:nh2
                                                        if wvny[j]>=div(2*wvnym,3)
                                                                                al[:,j] .= 0
                                                                            end
                                                    end
                                    ah[:,:,k]=fft_bac(al,2)
                                end
                end
        return nothing
    end

"""
    dealias!(ah::AbstractArray{<:Real,2},dims::Integer)

Dealiasing of a 2D array `ah` by removing the top 2/3 wavenumbers of
dimension `dims` and write to `ah`.

"""
function dealias!(ah::AbstractArray{<:Real,2},dims::Integer)

        if dims!=1 && dims!=2
                    error("Dims must be either 1 or 2!")
                end
        nh1=size(ah,1)
        nh2=size(ah,2)
        wvnx, wvny=wvn_serinit(nh1,nh2)
        wvnxm=div(nh1,2)
        wvnym=div(nh2,2)
        if dims==1
                    al = fft_for(ah,1)
                    for i=1:nh1
                                    if wvnx[i]>=div(2*wvnxm,3)
                                                        al[i,:] .= 0
                                                    end
                                end
                    ah[:,:]=fft_bac(al,1)
                else
                    al = fft_for(ah,2)
                    for j=1:nh2
                                    if wvny[j]>=div(2*wvnym,3)
                                                        al[:,j] .= 0
                                                    end
                                end
                    ah[:,:]=fft_bac(al,2)
                end
        return nothing
    end


"""
    dealias!(ah::AbstractArray{<:Real,2})

Dealiasing of a 2D array `ah` by removing the top 2/3 wavenumbers and write to `ah`.

"""
function dealias!(ah::AbstractArray{<:Real,2})
        nh1=size(ah,1)
        nh2=size(ah,2)

        wvnx, wvny=wvn_serinit(nh1,nh2)
        wvnxm=div(nh1,2)
        wvnym=div(nh2,2)
        fft_for!(ah,1)
        for i=1:nh1
                    if wvnx[i]>=div(2*wvnxm,3)
                                    ah[i,:] .= 0
                                end
                end
        fft_for!(ah,2)
        for j=1:nh2
                    if wvny[j]>=div(2*wvnym,3)
                                    ah[:,j] .= 0
                                end
                end
        fft_bac!(ah)

        return nothing
    end


"""
    filtering(ah::AbstractArray{<:Real,2},nl::Tuple{<:Integer,<:Integer})
Apply a low-pass filter on the 2D array `ah` and return a 2D array with size `nl`.
"""
function filtering(ah::AbstractArray{<:Real,2},nl::Tuple{<:Integer,<:Integer})

        al=zeros(nl)
        nh1=size(ah,1)
        nh2=size(ah,2)
        nl1=nl[1]
        nl2=nl[2]

        tmp1=fft_for(ah,1)
        tmp2=zeros(nl1,nh2)
        tmp2[1:div(nl1,2),:]=tmp1[1:div(nl1,2),:]
        tmp2[div(nl1,2)+1,:]=tmp1[div(nh1,2)+1,:]
        tmp2[div(nl1,2)+2:nl1,:]=tmp1[nh1-div(nl1,2)+2:nh1,:]
        tmp3=fft_for(tmp2,2)
        tmp4=zeros(nl1,nl2)
        tmp4[:,1:div(nl2,2)]=tmp3[:,1:div(nl2,2)]
        tmp4[:,div(nl2,2)+1]=tmp3[:,div(nh2,2)+1]
        tmp4[:,div(nl2,2)+2:nl2]=tmp3[:,nh2-div(nl2,2)+2:nh2]
        fft_bac!(al,fft_bac(tmp4,2),1)
        return al
    end

"""
    filtering(ah::AbstractArray{<:Real,2},ratio::Real)

Apply the low-pass filter to remove components with wavenumbers higher than `ratio` of the maximum wavenumber and return a 2D array of the same size.

"""
function filtering(ah::AbstractArray{<:Real,2},ratio::Real)
        al=copy(ah)
        filtering!(al,ratio)
        return al
    end

"""
    filtering!(ah::AbstractArray{<:Real,2},ratio::Real)

Apply the low-pass filter to remove components with wavenumbers higher than `ratio` of the maximum wavenumber and write to `ah`.

"""
function filtering!(ah::AbstractArray{<:Real,2},ratio::Real)
        nh1=size(ah,1)
        nh2=size(ah,2)

        wvnx, wvny=wvn_serinit(size(ah))
        wvnxm=div(nh1,2)
        wvnym=div(nh2,2)
        fft_for!(ah,1)
        for i=1:nh1
                    if wvnx[i]/wvnxm>=ratio
                                    ah[i,:] .= 0
                                end
                end
        fft_for!(ah,2)
        for j=1:nh2
                    if wvny[j]/wvnym>=ratio
                                    ah[:,j] .= 0
                                end
                end
        fft_bac!(ah)
        return nothing
    end

"""
    filtering(ah::AbstractArray{<:Real,3},ratio::Real)

Apply the low-pass filter in the first two dimensions of `ah` to
remove components with wavenumbers higher than `ratio` of the maximum
wavenumber and return a 3D array of the same size.

"""
function filtering(ah::AbstractArray{<:Real,3},ratio::Real)
        al=similar(ah)
        nh1=size(ah,1)
        nh2=size(ah,2)

        wvnx, wvny=wvn_serinit(nh1,nh2)
        wvnxm=div(nh1,2)
        wvnym=div(nh2,2)
        for k=1:size(ah,3)
                    al[:,:,k] =fft_for(ah[:,:,k],1)
                    for i=1:nh1
                                    if wvnx[i]/wvnxm>=ratio
                                                        al[i,:,k] .= 0
                                                    end
                                end
                    al[:,:,k]=fft_for(al[:,:,k],2)
                    for j=1:nh2
                                    if wvny[j]/wvnym>=ratio
                                                        al[:,j,k] .= 0
                                                    end
                                end
                    al[:,:,k]=fft_bac(al[:,:,k])
                end
        return al
    end

"""
    filtering!(ah::AbstractArray{<:Real,3},ratio::Real)

Apply the low-pass filter in the first two dimensions of `ah` to
remove components with wavenumbers higher than `ratio` of the maximum
wavenumber and write to `ah`.

"""
function filtering!(ah::AbstractArray{<:Real,3},ratio::Real)
        nh1=size(ah,1)
        nh2=size(ah,2)

        wvnx, wvny=wvn_serinit(nh1,nh2)
        wvnxm=div(nh1,2)
        wvnym=div(nh2,2)
        for k=1:size(ah,3)
                    al = fft_for(ah[:,:,k],1)
                    for i=1:nh1
                                    if wvnx[i]/wvnxm>=ratio
                                                        al[i,:] .= 0
                                                    end
                                end
                    fft_for!(al,2)
                    for j=1:nh2
                                    if wvny[j]/wvnym>=ratio
                                                        al[:,j] .= 0
                                                    end
                                end
                    ah[:,:,k]=fft_bac(al)
                end
        return nothing
    end

"""
    filtering!(al::AbstractArray{<:Real,2},ah::AbstractArray{<:Real,2})

In-place version.
"""
function filtering!(al::AbstractArray{<:Real,2},ah::AbstractArray{<:Real,2})
        nh1=size(ah,1)
        nh2=size(ah,2)
        nl1=size(al,1)
        nl2=size(al,2)

        tmp1=fft_for(ah,1)
        tmp2=zeros(nl1,nh2)
        tmp2[1:div(nl1,2),:]=tmp1[1:div(nl1,2),:]
        tmp2[div(nl1,2)+1,:]=tmp1[div(nh1,2)+1,:]
        tmp2[div(nl1,2)+2:nl1,:]=tmp1[nh1-div(nl1,2)+2:nh1,:]
        tmp3=fft_for(tmp2,2)
        tmp4=zeros(nl1,nl2)
        tmp4[:,1:div(nl2,2)]=tmp3[:,1:div(nl2,2)]
        tmp4[:,div(nl2,2)+1]=tmp3[:,div(nh2,2)+1]
        tmp4[:,div(nl2,2)+2:nl2]=tmp3[:,nh2-div(nl2,2)+2:nh2]
        fft_bac!(al,fft_bac(tmp4,2),1)
        return nothing
    end

"""
    filtering!(al::AbstractArray{<:Real,2},ah::AbstractArray{<:Real,2},ratio::Real)

Apply the low-pass filter to remove components with wavenumbers higher than `ratio` of the maximum wavenumber and write to `al`.
"""
function filtering!(al::AbstractArray{<:Real,2},ah::AbstractArray{<:Real,2},ratio::Real)
        nh1=size(ah,1)
        nh2=size(ah,2)

        wvnx, wvny=wvn_serinit(size(ah))
        wvnxm=div(nh1,2)
        wvnym=div(nh2,2)
        al .=fft_for(ah,1)
        for i=1:nh1
                    if wvnx[i]/wvnxm>=ratio
                                    al[i,:] .= 0
                                end
                end
        fft_for!(al,2)
        for j=1:nh2
                    if wvny[j]/wvnym>=ratio
                                    al[:,j] .= 0
                                end
                end
        fft_bac!(al)
        return nothing
    end

"""
    filtering(ah::AbstractArray{<:Real,1},nl::Integer)

Apply a low-pass filter on the 1D array `ah` and return a 1D array with size `nl`.

"""
function filtering(ah::AbstractArray{<:Real,1},nl::Integer)

        al=zeros(nl)
        nh=size(ah,1)

        tmp1=fft_for(ah)
        tmp2=zeros(nl)
        tmp2[1:div(nl,2)]=tmp1[1:div(nl,2)]
        tmp2[div(nl,2)+1]=tmp1[div(nh,2)+1]
        tmp2[div(nl,2)+2:nl]=tmp1[nh-div(nl,2)+2:nh]
        fft_bac!(al,tmp2)
        return al
    end

"""
    filtering!(al::AbstractArray{<:Real,1},ah::AbstractArray{<:Real,1})

Apply a low-pass filter on the 1D array `ah` and write to a 1D array with size `nl`.

"""
function filtering!(al::AbstractArray{<:Real,1},ah::AbstractArray{<:Real,1})

        nh=size(ah,1)
        nl=size(al,1)

        tmp1=fft_for(ah)
        tmp2=zeros(nl)
        tmp2[1:div(nl,2)]=tmp1[1:div(nl,2)]
        tmp2[div(nl,2)+1]=tmp1[div(nh,2)+1]
        tmp2[div(nl,2)+2:nl]=tmp1[nh-div(nl,2)+2:nh]
        fft_bac!(al,tmp2)
        return al
    end

"""
    filtering(ah::AbstractArray{<:Real,1},ratio::AbstractFloat)

Apply a low-pass filter on the 1D array `ah` and keep the percentage assigned by `ratio`.

"""
function filtering(ah::AbstractArray{<:Real,1},ratio::AbstractFloat)

        n=size(ah,1)
        al=fft_for(ah)
        wvnx=wvn_serinit(n)
        for i=1:n
                    if wvnx[i] >= ratio*(n/2)
                                    al[i] = 0
                                end
                end
        fft_bac!(al)
        return al
    end

"""
    filtering!(ah::AbstractArray{<:Real,1},ratio::AbstractFloat)

Apply a low-pass filter on the 1D array `ah` and keep the percentage assigned by `ratio`.

"""
function filtering!(ah::AbstractArray{<:Real,1},ratio::AbstractFloat)

        n=size(ah,1)
        fft_for!(ah)
        wvnx=wvn_serinit(n)
        for i=1:n
                    if wvnx >= ratio*(n/2)
                                    ah[i] = 0
                                end
                end
        fft_bac!(ah)
        return nothing
    end

"""
    filtering!(al::AbstractArray{<:Real,1},ah::AbstractArray{<:Real,1},ratio::AbstractFloat)

Apply a low-pass filter on the 1D array `ah` and keep the percentage assigned by `ratio`.

"""
function filtering!(al::AbstractArray{<:Real,1},ah::AbstractArray{<:Real,1},ratio::AbstractFloat)

        n=size(ah,1)
        fft_for!(al,ah)
        wvnx=wvn_serinit(n)
        for i=1:n
                    if wvnx[i] >= ratio*(n/2)
                                    al[i] = 0
                                end
                end
        fft_bac!(al)
        return nothing
    end


"""
    padding(al::AbstractArray{<:Real,1},nh::Integer})

Return a 1D array of size `nh` with padded zeros at the high wavenumbers of a 1D array `al`.
"""
function padding(al::AbstractArray{<:Real,1},nh::Integer)

        ah=zeros(nh)
        nl=size(al,1)

        tmp1=fft_for(al)
        ah[1:div(nl,2)]=tmp1[1:div(nl,2)]
        ah[nh-div(nl,2)+2:nh]=tmp1[nl-div(nl,2)+2:nl]
        fft_bac!(ah)
        return ah
    end


"""
    padding(al::AbstractArray{<:Real,2},nh::Tuple{<:Integer,<:Integer})

Return a 2D array of size `nh` with padded zeros at the high wavenumbers of a 2D array `al`.
"""
function padding(al::AbstractArray{<:Real,2},nh::Tuple{<:Integer,<:Integer})

        ah=zeros(nh)
        nh1=nh[1]
        nh2=nh[2]
        nl1=size(al,1)
        nl2=size(al,2)

        tmp1=fft_for(al,1)
        tmp2=zeros(nh1,nl2)
        tmp2[1:div(nl1,2),:]=tmp1[1:div(nl1,2),:]
        tmp2[nh1-div(nl1,2)+2:nh1,:]=tmp1[nl1-div(nl1,2)+2:nl1,:]
        tmp3=fft_for(tmp2,2)
        tmp4=zeros(nh1,nh2)
        tmp4[:,1:div(nl2,2)]=tmp3[:,1:div(nl2,2)]
        tmp4[:,nh2-div(nl2,2)+2:nh2]=tmp3[:,nl2-div(nl2,2)+2:nl2]
        fft_bac!(ah,fft_bac(tmp4,2),1)
        return ah
    end


"""
    padding!(ah::AbstractArray{<:Real,2},al::AbstractArray{<:Real,2})

In-place version of `padding`.
"""
function padding!(ah::AbstractArray{<:Real,2},al::AbstractArray{<:Real,2})
        nh1=size(ah,1)
        nh2=size(ah,2)
        nl1=size(al,1)
        nl2=size(al,2)

        tmp1=fft_for(al,1)
        tmp2=zeros(nh1,nl2)
        tmp2[1:div(nl1,2),:]=tmp1[1:div(nl1,2),:]
        tmp2[nh1-div(nl1,2)+2:nh1,:]=tmp1[nl1-div(nl1,2)+2:nl1,:]
        tmp3=fft_for(tmp2,2)
        tmp4=zeros(nh1,nh2)
        tmp4[:,1:div(nl2,2)]=tmp3[:,1:div(nl2,2)]
        tmp4[:,nh2-div(nl2,2)+2:nh2]=tmp3[:,nl2-div(nl2,2)+2:nl2]
        fft_bac!(ah,fft_bac(tmp4,2),1)
        return nothing
    end

function spec2phy(a,b,c,d)

        #     by xuanting hao,
        #     04/23/2014   first edition

        #     cc/4, -sc/4      a   b
        #     -cs/4,ss/4       c   d

        #      4a  =  cc = a1*cos(alpha) + a2*cos(beta)
        #     -4b  =  sc = -a1*sin(alpha) - a2*sin(beta)
        #     -4c  =  cs = -a1*sin(alpha) + a2*sin(beta)
        #      4d  =  ss = -a1*cos(alpha) + a2*cos(beta)
        #     a1^2 + a2^2 = 8 * (a^2+b^2+c^2+d^2)

        a1 = 2*sqrt((b + c)^2 + (a - d)^2)
        a2 = 2*sqrt((b - c)^2 + (a + d)^2)
        alpha = abs(acos((a - d) / a1)) * sign(b + c)
        beta  = abs(acos((a + d) / a2)) * sign(b - c)

        return a1,a2,alpha,beta
    end


end



