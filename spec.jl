using SpecialFunctions
using Test
using YAML

function readwaveparams(fname="wave.yml"::String)
    data = YAML.load_file(fname)

    key="physical"
    u10 = data[key]["u10"]
    γ = data[key]["gamma"]
    fetch = data[key]["fetch"]
    nd = data[key]["nd"]

    key="numerical"
    nx = data[key]["nx"]
    ny = data[key]["ny"]
    pex = data[key]["pex"]
    pey = data[key]["pey"]
    
    @test u10 > 0
    @test γ >= 1 && γ <=3.3
    @test fetch >= 0
    @test nd>0
    @test nx > 1
    @test ny > 1
    @test pex > 0
    @test pey > 0
    return u10,γ,fetch,nd,nx,ny,pex,pey
end

function writewaveparams(u10::Union{AbstractFloat,Int},γ::Union{AbstractFloat,Int},
                         fetch::Union{AbstractFloat,Int},nd::Integer,nx::Integer,ny::Integer,
                         pex::Union{AbstractFloat,Int},pey::Union{AbstractFloat,Int},fname="wave.yml"::String)
    data=Dict("physical"=>Dict("u10"=>u10,"gamma"=>γ,"fetch"=>fetch,"nd"=>nd),
              "numerical"=>Dict("nx"=>nx,"ny"=>ny,"pex"=>pex,"pey"=>pey))
    writewaveparams(data,fname)
    return nothing
end

function writewaveparams(data::Dict,fname="wave.yml"::String)
    YAML.write_file(fname,data)
    return nothing
end

function jonsparam(a::Union{AbstractFloat,Int},b::Union{AbstractFloat,Int},mode::Integer)

    g=9.8
    if mode==1
        α=a
        fp=b
        ωp=2*π*fp
    end
    if mode==2
        α=a
        λp=b
        ωp=sqrt(g*2*π/λp)
    end
    if mode==3
        α=a
        kp=b
        ωp=sqrt(g*kp)
    end
    u10=22*g/ωp*(α/0.076)^(1/0.66)
    fetch=(22/ωp)^3*g^2/u10
    return u10,fetch
end

function jonswap2d(u10::Union{AbstractFloat,Int},fetch::Union{AbstractFloat,Int},
                   γ::Union{AbstractFloat,Int},nx::Integer,nd::Union{AbstractFloat,Int},
                   pex::Union{AbstractFloat,Int};itest=false::Bool)

    g=9.8
    α = 0.076 * (u10^2 / fetch / g)^0.22
    ωp = 22.0 * (g^2 / u10 / fetch)^(1/3)
    fr2=1/g

    dx=2*pi/pex/nx
    η=zeros(nx)
    swk=zeros(nx)
    allkx=[kx*pex for kx=1:nx]
    for kx=1:div(nx,2)-1
        wvn=kx*pex
        ϕ=rand()*2*pi
        ω=sqrt(wvn/fr2)
        if ω<ωp
            σ=0.07
        else
            σ=0.09
        end
        sw = (α/fr2^2/ω^5) * exp(-1.25*(ω/ωp)^(-4))*γ^(exp(-(ω-ωp)^2/2/σ^2/ωp^2))
        sw = sw/fr2/2/ω
        swk[kx]=sw 
        amp = sqrt(2*sw*pex)
        for i=1:nx
            x=(i-1)*dx
            η[i] = η[i] + amp*cos(wvn*x+ϕ)
        end
    end
    if itest
        return η,allkx,swk
    else
        return η
    end
end

function linear3d(ka::Real,nkx::Integer,nky::Integer,nx::Integer,ny::Integer,
                  pex::Union{AbstractFloat,Int},pey::Union{AbstractFloat,Int},θ=0::Real)

    kx=nkx*pex
    ky=nky*pey
    amp = ka / sqrt(kx^2+ky^2)

    lx=2*π/pex
    ly=2*π/pey
    x=[(i-1)*lx/nx for i=1:nx]
    y=[(j-1)*ly/ny for j=1:ny]
    η=zeros(nx,ny)

    for j=1:ny, i=1:nx
        η[i,j] = amp * cos(nkx*pex*x[i]+nky*pey*y[j]+θ)
    end
    return η
end


function jonswap3d(u10::Union{AbstractFloat,Int},fetch::Union{AbstractFloat,Int},
                   γ::Union{AbstractFloat,Int},nx::Integer,ny::Integer,
                   nd::Union{AbstractFloat,Int},pex::Union{AbstractFloat,Int},
                   pey::Union{AbstractFloat,Int};itest=false::Bool)

    g=9.8
    α = 0.076 * (u10^2 / fetch / g)^0.22
    ωp = 22.0 * (g^2 / u10 / fetch)^(1/3)
    fr2=1/g
    lx=2*π/pex
    ly=2*π/pey
    x=[(i-1)*lx/nx for i=1:nx]
    y=[(j-1)*ly/ny for j=1:ny]
    if itest
        swk=zeros(div(nx,2),ny+1)
        allkx=zeros(div(nx,2),ny+1)
        allky=zeros(div(nx,2),ny+1)
        for i=1:div(nx,2)
            allkx[i,:].= i*pex
        end
        for j=1:ny+1
            allky[:,j].= (j-div(ny,2)-1)*pey 
        end
    end
    η=zeros(nx,ny)
    for kx=0:div(nx,2)-1
        for ky=1:2:ny
            kyr=div(ky+1,2) - 1
            if kx + kyr > 0 
                wvn=sqrt((kx*pex)^2+(kyr*pey)^2)
                ϕ1,ϕ2=rand(2) .* 2 * pi
                ω=sqrt(wvn/fr2)
                if ω<ωp
                    σ=0.07
                else
                    σ=0.09
                end
                θ=atan(kyr,kx)
                Dθ=spreading(θ,nd)
                sw = (α/fr2^2/ω^5) * exp(-1.25*(ω/ωp)^(-4))*γ^(exp(-(ω-ωp)^2/2/σ^2/ωp^2))*Dθ
                sw = sw/fr2^2/2/ω^3                
                if itest && kx>0
                    swk[kx,kyr+div(ny,2)+1]=sw
                    if kyr<div(ny,2)
                        swk[kx,-kyr+div(ny,2)+1]=sw
                    end
                end
                amp = sqrt(2*sw*pex*pey)

                if kyr==0
                    ϕ2 = ϕ1
                end
                amp = amp / 4
                if kx==0 || kyr==0
                    amp=amp*sqrt(2)
                end

                #Physical space
                # amp=amp*sqrt(2)*2
                # if kx==0 || kyr==0
                #     for i=1:nx,j=1:ny
                #         η[i,j]=η[i,j]+amp*cos(kx*pex*x[i]+kyr*pey*y[j]+ϕ1)
                #     end               
                # else
                #     for i=1:nx,j=1:ny
                #         η[i,j]=η[i,j]+amp*cos(kx*pex*x[i]+kyr*pey*y[j]+ϕ1)+amp*cos(kx*pex*x[i]-kyr*pey*y[j]+ϕ2)
                #     end
                # end

                i1 = kx+1
                i2 = nx+1-kx
                j1 = kyr+1
                j2 = ny+1-kyr    
                η[i1,j1] = amp * cos(ϕ1) + amp * cos(ϕ2)
                if i2 <= nx 
                    η[i2,j1] = amp * sin(ϕ1) + amp * sin(ϕ2)
                end
                if j2 <= ny
                    η[i1,j2] = amp * sin(ϕ1) - amp * sin(ϕ2)
                    if i2 <= nx
                        η[i2,j2] = -amp * cos(ϕ1) + amp * cos(ϕ2)
                    end
                end 
            end
        end
    end
    fftbac!(η)
    if itest
        return η,allkx,allky,swk
    else
        return η
    end
end

function eckv2d(u10::Union{AbstractFloat,Int},Ωc::Union{AbstractFloat,Int},
                nx::Integer,pex::Union{AbstractFloat,Int};itest=false::Bool)

    α = 0.0081
    β = 1.25
    g = 9.82
    Cd10N=0.00144
    ustar=sqrt(Cd10N)*u10
    ao=0.1733
    ap=4.0
    km=370
    cm=0.23
    am=0.13*ustar/cm
    if Ωc<=1
        γ=1.7
    else
        γ=1.7+6*log10(Ωc)
    end
    σ=0.08*(1+4/Ωc^3)
    αp=0.006*Ωc^0.55
    if ustar<=cm
        αm=0.01*(1+log(ustar/cm))
    else
        αm=0.01*(1+3*log(ustar/cm))
    end
    ko=g/u10^2
    kp=ko*Ωc^2
    cp=sqrt(g/kp)
    ωp=kp*cp

#    dx=2*pi/pex/nx              
    η=zeros(nx)
    if itest
        swk=zeros(div(nx,2))
        allkx=[kx*pex for kx=1:div(nx,2)]
    end
    for kx=1:div(nx,2)-1
        k=kx*pex
        ϕ=rand()*2*pi
        c=sqrt((g/k)*(1+(k/km)^2))

        LPM=exp(-1.25*(kp/k)^2)
        Γ=exp((sqrt(k/kp)-1)^2/(-2*σ^2))
        Jp=γ^Γ
        Fp=LPM*Jp*exp(-0.3162*Ωc*(sqrt(k/kp)-1))
        Fm=LPM*Jp*exp(-0.25*(k/km-1)^2)
        Bl=0.5*αp*(cp/c)*Fp
        Bh=0.5*αm*(cm/c)*Fm
        if Bl<0
            Bl=0
        end
        if Bh<0
            Bh=0
        end
        sw=(Bl+Bh)/k^3
        if itest
            swk[kx]=sw
        end
        amp = sqrt(2*sw*pex)
        # for i=1:nx
        #     x=(i-1)*dx
        #     η[i] = η[i] + amp*cos(k*x+ϕ)
        # end

        amp=amp/2
        i1 = kx+1
        i2 = nx+1-kx
        η[i1] = amp * cos(ϕ) 
        if i2 <= nx
            η[i2] = amp * sin(ϕ)
        end
    end
    fftbac!(η)
    if itest
        return η,allkx,swk
    else
        return η
    end   
    
end

function eckv3d(u10::Union{AbstractFloat,Int},Ωc::Union{AbstractFloat,Int},
                nx::Integer,ny::Integer,nd::Union{AbstractFloat,Int},
                pex::Union{AbstractFloat,Int},pey::Union{AbstractFloat,Int};itest=false::Bool)

    α = 0.0081
    β = 1.25
    g = 9.82
    Cd10N=0.00144
    ustar=sqrt(Cd10N)*u10
    ao=0.1733
    ap=4.0
    km=370
    cm=0.23
    am=0.13*ustar/cm
    if Ωc<=1
        γ=1.7
    else
        γ=1.7+6*log10(Ωc)
    end
    σ=0.08*(1+4/Ωc^3)
    αp=0.006*Ωc^0.55
    if ustar<=cm
        αm=0.01*(1+log(ustar/cm))
    else
        αm=0.01*(1+3*log(ustar/cm))
    end
    ko=g/u10^2
    kp=ko*Ωc^2
    cp=sqrt(g/kp)
    ωp=kp*cp

    # dx=2*pi/pex/nx
    # dy=2*pi/pey/ny
    η=zeros(nx,ny)
    if itest
        swk=zeros(div(nx,2),ny+1)
        allkx=zeros(div(nx,2),ny+1)
        allky=zeros(div(nx,2),ny+1)
        for i=1:div(nx,2)
            allkx[i,:].= i*pex
        end
        for j=1:ny+1
            allky[:,j].= (j-div(ny,2)-1)*pey
        end
    end
    for kx=0:div(nx,2)-1
        for ky=1:2:ny
            kyr=div(ky+1,2) - 1
            if kx + kyr > 0
                k=sqrt((kx*pex)^2+(kyr*pey)^2)
                θ=atan(kyr,kx)
                Dθ=spreading(θ,nd)
                ϕ1,ϕ2=rand(2) .*2*pi
                c=sqrt((g/k)*(1+(k/km)^2))
                
                LPM=exp(-1.25*(kp/k)^2)
                Γ=exp((sqrt(k/kp)-1)^2/(-2*σ^2))
                Jp=γ^Γ
                Fp=LPM*Jp*exp(-0.3162*Ωc*(sqrt(k/kp)-1))
                Fm=LPM*Jp*exp(-0.25*(k/km-1)^2)
                Bl=0.5*αp*(cp/c)*Fp
                Bh=0.5*αm*(cm/c)*Fm
                if Bl<0
                    Bl=0
                end
                if Bh<0
                    Bh=0
                end

                sw=(Bl+Bh)/k^3                
                sw=(sw/k)*Dθ
                
                if itest && kx>0
                    swk[kx,kyr+div(ny,2)+1]=sw
                    if kyr<div(ny,2)
                        swk[kx,-kyr+div(ny,2)+1]=sw
                    end
                end
                amp = sqrt(2*sw*pex*pey)
                
                if kyr==0
                    ϕ2 = ϕ1
                end
                amp = amp / 4
                if kx==0 || kyr==0
                    amp=amp*sqrt(2)
                end

                i1 = kx+1
                i2 = nx+1-kx
                j1 = kyr+1
                j2 = ny+1-kyr
                η[i1,j1] = amp * cos(ϕ1) + amp * cos(ϕ2)
                if i2 <= nx
                    η[i2,j1] = amp * sin(ϕ1) + amp * sin(ϕ2)
                end
                if j2 <= ny
                    η[i1,j2] = amp * sin(ϕ1) - amp * sin(ϕ2)
                    if i2 <= nx
                        η[i2,j2] = -amp * cos(ϕ1) + amp * cos(ϕ2)
                    end
                end
            end
        end
    end
    fftbac!(η)
    if itest
        return η,allkx,allky,swk
    else
        return η
    end       
end

function donelan3d(u10::Union{AbstractFloat,Int},ucp::Union{AbstractFloat,Int},
                   nx::Integer,ny::Integer,pex::Union{AbstractFloat,Int},pey::Union{AbstractFloat,Int})

    if ucp>5 || ucp < 0.83
        error("Parameter ucp not satisfying: 0.83<U/cp<5!")
    end
    g=9.8
    fr2=1/g
    cp=u10/ucp
    kp=1/cp^2/fr2
    ωp=sqrt(kp/fr2)
    α=0.006*ucp^0.55
    if ucp<1
        γ=1.7
    else
        γ=1.7+6*log10(ucp)
    end
    σ=0.08*(1+4/ucp^3)

    η=zeros(nx,ny)
    for kx=0:div(nx,2)-1
        for ky=1:2:ny
            kyr=div(ky+1,2) - 1
            if kx + kyr > 0 
                wvn=((kx*pex)^2+(kyr*pey)^2)^0.5
                θ=atan(kyr,kx)
                ϕ1,ϕ2=rand(2) .* 2 * pi
                ω=(wvn/fr2)^0.5
                if ω/ωp>0.56 && ω/ωp<0.95
                    β=2.61*(ω/ωp)^1.3
                else 
                    if ω/ωp>0.95 && ω/ωp<1.6
                        β=2.28*(ω/ωp)^(-1.3)
                    else
                        β=1.24
                    end
                end
                Dθ=0.5*β*(sech(β*θ))^2
                sw = α/fr2^2/ω^5 * (ω/ωp)* exp(-1*(ωp/ω)^4)*γ^(exp(-(ω-ωp)^2/2/σ^2/ωp^2))*Dθ
                sw = sw/fr2^2/2/ω^3
                amp = sqrt(2*sw*pex*pey)

                i1 = kx+1
                i2 = nx+1-kx
                j1 = kyr+1
                j2 = ny+1-kyr    
                η[i1,j1] = amp * cos(ϕ1) + amp * cos(ϕ2)
                if i2 <= nx 
                    η[i2,j1] = amp * sin(ϕ1) + amp * sin(ϕ2)
                end
                if j2 <= ny
                    η[i1,j2] = amp * sin(ϕ1) - amp * sin(ϕ2)
                    if i2 <= nx
                        η[i2,j2] = -amp * cos(ϕ1) + amp * cos(ϕ2)
                    end
                end 
            end
        end
    end
    fftbac!(η)
    return η
end

function spreading(θ::Union{AbstractFloat,Int},nd::Union{AbstractFloat,Int})
    return (1/sqrt(pi))*gamma(nd/2+1.0)/gamma(nd/2+0.5)*cos(θ)^nd
end

