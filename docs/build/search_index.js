var documenterSearchIndex = {"docs":
[{"location":"module/refraction.html#Refraction","page":"Refraction","title":"Refraction","text":"","category":"section"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"Photons refraction between air and water ","category":"page"},{"location":"module/refraction.html#Snell's-Law","page":"Refraction","title":"Snell's Law","text":"","category":"section"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"angle of reflection [1]","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"theta_r = cos^-1hatxicdothatn","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"When our incoming photon is directly downward: hatxi=beginbmatrix001endbmatrix and the normal vector can be defined by hatn=dfrac1sqrt1+left(eta_xright)^2+left(eta_yright)^2left(-eta_xhati-eta_yhatj+hatkright) When, eta_x is the partial derivative the partial derivative of eta with respect to x: dfracpartialetapartial x and eta_y is the partial derivative the partial derivative of eta with respect to y: dfracpartialetapartial y","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"Therefore, the angle of reflection, in this module, can be described by ","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"theta_r = cos^-1left(frac1sqrt1+left(eta_xright)^2+left(eta_yright)^2right)","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"angle of transmission","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"theta_t = sin^-1left(frac1n_wsintheta_rright)","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"Substitute the theta_r that we found above. ","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"theta_t = sin^-1left(frac1n_wsqrtfracleft(eta_xright)^2+left(eta_yright)^21+left(eta_xright)^2+left(eta_yright)^2right)","category":"page"},{"location":"module/refraction.html#Fresnel-Reflectance","page":"Refraction","title":"Fresnel Reflectance","text":"","category":"section"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"In our package, we calculate the energy proportion of the light ray that being transmitted to the water, but first we identify the amplitude transmission coefficient or the ratio between electric field amplitude of the transmitted light ray to the intirial light ray. [2]","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"t_perp=left(fracE_tE_0right)_perp=frac2sin(theta_t)cos(theta_r)sintheta_t+theta_r","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"t_parallel=left(fracE_tE_0right)_parallel=frac2sin(theta_t)cos(theta_r)sintheta_t+theta_rcostheta_r-theta_t","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"When t_perp is corresponding to the amplitude transmission coefficient of the light ray in which the electric field, that constitute the electro magnetic wave, perpendicular to the plane-of-incident , and t_parallel is corresponding to the amplitude transmission coefficient of the light ray in which the electric field, that constitute the electro magnetic wave, parallels to the plane-of-incident.","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"The intensity of the light is proportional to the square of the amplitude, I propto E^2.  All conceivable azimuths of waves that are polarized combine to form natural or unpolarized light. Each wave can be broken down into its constituent parts. Each component will have an equal amount due to symmetry. Then, half of the amplitude transmission coefficient yields the transmission coefficient of a surface in natural light. [3]","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"t = fracI_tI_0 = frac12leftleftfrac2sin(theta_t)cos(theta_r)sin(theta_r+theta_t)right^2+leftfrac2sin(theta_t)cos(theta_r)sin(theta_r+theta_t)cos(theta_r-theta_t)right^2right","category":"page"},{"location":"module/refraction.html#result","page":"Refraction","title":"result","text":"","category":"section"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"We, then, transform the reflection angle into the azimuthal angle and polar angle in the spherical coordination. ","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"temx = -fraceta_xsqrt(eta_x)^2+(eta_y)^2","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"temy = -fraceta_ysqrt(eta_x)^2+(eta_y)^2","category":"page"},{"location":"module/refraction.html#Reference","page":"Refraction","title":"Reference","text":"","category":"section"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"[1]: Mobley, C. (1994). Across the Surface. Light and Water: Radiative Transfer in Natural Waters (pp. 155-157). Academic Press. ","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"[2]: Hecht, E. (2001). The Propagation of Light. Optics (pp. 113-115). Addison-Wesley. ","category":"page"},{"location":"module/refraction.html","page":"Refraction","title":"Refraction","text":"[3]: Sears, F.W. (1949). Polarization Optics (pp. 173-174). Addison-Wesley","category":"page"},{"location":"QuickStart/Center.html#Quick-Start","page":"Photons at the center","title":"Quick Start","text":"","category":"section"},{"location":"QuickStart/Center.html","page":"Photons at the center","title":"Photons at the center","text":"import Pkg \nPkg.activate(\"https://github.com/haoboatlab/Light.git\")\nimport Light","category":"page"},{"location":"QuickStart/Center.html#Initial-Condition","page":"Photons at the center","title":"Initial Condition","text":"","category":"section"},{"location":"QuickStart/Center.html","page":"Photons at the center","title":"Photons at the center","text":"For the input initial condition, the function readparam will automatically read the input from light.yml file.  Below, is the structure contained in the light.yml file.","category":"page"},{"location":"QuickStart/Center.html","page":"Photons at the center","title":"Photons at the center","text":"using YAML\n\nYAML.write_file(\"light.yml\",\"\nirradiance:\n  nz: 200\n  dz: 1\n  nxe: 512\n  nye: 512\n  num: 31\n  ztop: 10\nphoton:\n  nphoton: 100000\n  kr: 10\n  nxp: 512\n  kbc: 0\n  b: 0.006\n  nyp: 512\n  a: 0.007\nwave:\n  pey: 0.07853981633974483\n  nxeta: 512\n  nyeta: 512\n  pex: 0.07853981633974483\n\")","category":"page"},{"location":"QuickStart/Center.html","page":"Photons at the center","title":"Photons at the center","text":"We, then, store all the input variable in the struct p, by using the function below. ","category":"page"},{"location":"QuickStart/Center.html","page":"Photons at the center","title":"Photons at the center","text":"p = readparams()","category":"page"},{"location":"QuickStart/Center.html","page":"Photons at the center","title":"Photons at the center","text":"allind=1:p.nphoton\n\"η is the height(z axis) corresponding to each grid point on the surface, 2d array of grid number in x and y direction\"\nη=zeros(p.nxs,p.nys)\n\"ηx is the x coordination corresponding to each grid point on the surface, 2d array of grid number in x and y direction\"\nηx=zeros(p.nxs,p.nys)\n\"ηy is the y coordination corresponding to each grid point on the surface, 2d array of grid number in x and y direction\"\nηy=zeros(p.nxs,p.nys)\nϕps,θps=phasePetzold()","category":"page"},{"location":"QuickStart/Center.html","page":"Photons at the center","title":"Photons at the center","text":"ed=zeros(p.nx,p.ny,p.nz)\nesol=zeros(p.num,p.nz)\narea=zeros(4)\ninteri=zeros(Int64,4)\ninterj=zeros(Int64,4)\nxpb=zeros(p.nxp,p.nyp)\nypb=zeros(p.nxp,p.nyp)\nzpb=zeros(p.nxp,p.nyp)\nθ=zeros(p.nxp,p.nyp)\nϕ=zeros(p.nxp,p.nyp)\nfres=zeros(p.nxp,p.nyp)\nix=div(p.nxη,2)+1\niy=div(p.nyη,2)+1","category":"page"},{"location":"QuickStart/Center.html#Run-the-Monte-Carlo-Simulation","page":"Photons at the center","title":"Run the Monte Carlo Simulation","text":"","category":"section"},{"location":"module/Parameters.html#Parameters","page":"Parameters","title":"Parameters","text":"","category":"section"},{"location":"module/Parameters.html#Coordinate-system-and-notation","page":"Parameters","title":"Coordinate system and notation","text":"","category":"section"},{"location":"module/Parameters.html","page":"Parameters","title":"Parameters","text":"LightMC.jl is formulated in the spherical system hatzeta = (thetaphi), where polar angle theta is measured from the direction of hatz and the azimuthal ange phi is measured positive counter clockwise from hatx, when looking toward the origin along hatz. Let hatxi denoted a unit vector pointing in the desired direction, when hatxi=left(xi_xxi_yxi_zright), and becasue hatxi is of unit length, its component satisfy hatxi_x^2+hatxi_y^2+hatxi_z^2=1. [1] Therefore, ","category":"page"},{"location":"module/Parameters.html","page":"Parameters","title":"Parameters","text":"hatxi = beginbmatrix sin(theta)cos(phi) sin(theta)sin(phi) cos(theta) endbmatrix","category":"page"},{"location":"module/Parameters.html","page":"Parameters","title":"Parameters","text":"To simplify the term above, we simplify hatxi by using the cosine parameter.","category":"page"},{"location":"module/Parameters.html","page":"Parameters","title":"Parameters","text":"hatxi = beginbmatrixmu_x mu_y mu_z endbmatrix = beginbmatrixsin(theta)cos(phi) sin(theta)sin(phi) cos(theta) endbmatrix ","category":"page"},{"location":"module/Parameters.html#Local-Coordinate-system","page":"Parameters","title":"Local Coordinate system","text":"","category":"section"},{"location":"module/Parameters.html","page":"Parameters","title":"Parameters","text":"When we calculate for the scattering direction, our result is in the local coordination system (hat(theta)hat(phi)hat(r)), when radial unit vector hat(r) is the same initial direction of photons before scattering hatxi, the azimuthal unit vector hat(phi) is defined by the cross product of the ocean coordinate system hatz and the incident vector's direction hatphi=frachatztimeshatrhatztimeshatr, and polar unit vector is given by hattheta=hatphitimeshatr.  Therefore, the unit vector of the scattered direction of photons hatxi_s can be described in the local coordination system (hat(theta)hat(phi)hat(r)) as,","category":"page"},{"location":"module/Parameters.html","page":"Parameters","title":"Parameters","text":"hatxi_(s) = beginbmatrix sin(theta_s)cos(phi_s) sin(theta_s)sin(phi_s) cos(theta_s) endbmatrix","category":"page"},{"location":"module/Parameters.html","page":"Parameters","title":"Parameters","text":"when, theta_s and phi_s is polar angle and azimuthal angle in relative to the local coordinate system (hat(theta)hat(phi)hat(r)). ","category":"page"},{"location":"module/Parameters.html","page":"Parameters","title":"Parameters","text":"To change from the local coordinate system to the cartesian coordination in the global system, we multiply hatxi_(s) by the basis of our local coordinate system.","category":"page"},{"location":"module/Parameters.html","page":"Parameters","title":"Parameters","text":"beginbmatrix mu_x mu_y mu_z endbmatrix = beginbmatrixhattheta  hatphi  hatr endbmatrixbeginbmatrixhatxi_(s)endbmatrix","category":"page"},{"location":"module/Parameters.html","page":"Parameters","title":"Parameters","text":"And, after we do the cross product and substitute hatxi_(s). ","category":"page"},{"location":"module/Parameters.html","page":"Parameters","title":"Parameters","text":"beginbmatrix mu_x mu_y mu_z endbmatrix = beginbmatrixfracmu_xmu_zsqrt1-mu_z^2-fracmu_ysqrt1-mu_z^2mu_xfracmu_ymu_zsqrt1-mu_z^2fracmu_xsqrt1-mu_z^2mu_y-sqrt1-mu_z^20mu_z endbmatrixbeginbmatrix sin(theta_s)cos(phi_s) sin(theta_s)sin(phi_s) cos(theta_s) endbmatrix","category":"page"},{"location":"module/Parameters.html#Reference","page":"Parameters","title":"Reference","text":"","category":"section"},{"location":"module/Parameters.html","page":"Parameters","title":"Parameters","text":"[1]: Mobley, C. (2021). Light and Radiometry. In A Ocean Optics Web Book. https://www.oceanopticsbook.info","category":"page"},{"location":"reference.html#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference.html","page":"Reference","title":"Reference","text":"Mobley, C. (2021). Light and Radiometry. In A Ocean Optics Web Book. https://www.oceanopticsbook.info","category":"page"},{"location":"module/MonteCarlo.html#Fundamental-Principle-of-Monte-Carlo-Simulation","page":"Monte Carlo Simulation","title":"Fundamental Principle of Monte Carlo Simulation","text":"","category":"section"},{"location":"module/MonteCarlo.html","page":"Monte Carlo Simulation","title":"Monte Carlo Simulation","text":"p_Re (Re) equiv begincases\n1  quad 0 leq Re leq 1 \n0  quad Re  0 cup Re  1\nendcases","category":"page"},{"location":"module/MonteCarlo.html#The-Optical-Path-Length","page":"Monte Carlo Simulation","title":"The Optical Path Length","text":"","category":"section"},{"location":"module/MonteCarlo.html","page":"Monte Carlo Simulation","title":"Monte Carlo Simulation","text":"r=-frac1clnRe","category":"page"},{"location":"module/MonteCarlo.html#Sampling-Scattering-direction","page":"Monte Carlo Simulation","title":"Sampling Scattering direction","text":"","category":"section"},{"location":"module/MonteCarlo.html","page":"Monte Carlo Simulation","title":"Monte Carlo Simulation","text":"azimuthal angle in the plane of the scattering event relative to the direction of photons before scattering varphi is uniformly distribute between 0 and 2pi Therefore,","category":"page"},{"location":"module/MonteCarlo.html","page":"Monte Carlo Simulation","title":"Monte Carlo Simulation","text":"varphi = 2piRe","category":"page"},{"location":"module/MonteCarlo.html","page":"Monte Carlo Simulation","title":"Monte Carlo Simulation","text":"To find the angle between new trajectory and the direction of photons before scattering, we use the Petzold..","category":"page"},{"location":"module/MonteCarlo.html","page":"Monte Carlo Simulation","title":"Monte Carlo Simulation","text":"phasePetzold()","category":"page"},{"location":"module/MonteCarlo.html#LightMC.phasePetzold-Tuple{}","page":"Monte Carlo Simulation","title":"LightMC.phasePetzold","text":"phasePetzold()\n\nreturn 2 arrays: ϕps and θps. When ϕps is the cumulation distribution of scattering angle and θps is the angle between new trajectory and the direction of the photon before scattering corresponding to each ϕps For example, the probability to find the scattering photon from the direction of the photon before scattering to θps[1] is ϕps[1] The data come from Kirk,1981, Monte Carlo Procedure for Simulating the Penetration of Light into Natural Waters\n\n\n\n\n\n","category":"method"},{"location":"index.html#LightMC.jl-Documentation","page":"Home","title":"LightMC.jl Documentation","text":"","category":"section"},{"location":"index.html#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"LightMC.jl is a collection of functions that can be used to calculate the downaward irradiance field and run the Monte Carlo Simulation","category":"page"},{"location":"index.html#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"using Pkg \nPkg.activate(\"https://github.com/haoboatlab/Light.git\") ","category":"page"}]
}
