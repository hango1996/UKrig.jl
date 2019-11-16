module UKrig

using SpecialFunctions: besselk, gamma
import DynamicPolynomials: @polyvar, monomials
using FixedPolynomials
using LinearAlgebra
using Statistics
using Random

export Mnu, Gnu, generate_Gnu_krigY, generate_Mnu_krigY

tν𝒦t(t,ν) = t^ν * besselk(ν, t)

function Mnu(t, ν)::Float64
	pt, pν, p0, p1 = promote(t, ν, Float64(0), Float64(1)) # note besselk always returns a Float64 apparently
	return (pt==p0) ? p1 : tν𝒦t(√(2pν)*pt,pν) * 2^(1-pν) / gamma(pν)
end

# the const on the principle irregular term
scν(ν)       = - (2ν)^ν * gamma(ν + 1//2) * gamma(1-ν) / gamma(2ν+1) / sqrt(π)
scν(ν::Int)  = - 2 * (-2ν)^ν * gamma(ν + 1//2) / gamma(ν) / gamma(2ν+3) / sqrt(π)

function Gnu(t::T, ν::Int) where T<:Real
	if t==0
		return T(0)
	end
	return scν(ν) * t^(2ν) * log(t)
end

# Gnu(t::T, ν) where T<:Real = scν(ν) * t^(2ν)
function Gnu(t::T, ν) where T<:Real
	return scν(ν) * t^(2ν)
end


"""
`generate_Gnu_krigY(;fdata, xdata, ν, σg, σε) -> (x::Array->krigY(x), x::Array->fp(x'), b, c)`
"""
function generate_Gnu_krigY(Ydata,Xdata1,Xdata2,Xdata3, ν)
	m   =  7
	n₁  = length(Xdata1)
	#kriging matrix F
	@polyvar x y z
	p=1+x+y+x^2+y^2+x*y+z
	F=[t(x=>vx,y=>vy,z=>vz) for (vx, vy,vz) in zip(Xdata1,Xdata2,Xdata3), t in p]
	#distance matrix and Gnu
 	dist_mat= sqrt.(norm.( (Xdata1.-Xdata1').^2  + (Xdata2.-Xdata2').^2) )
	Gnu_mat=Gnu.(dist_mat, ν) |>Symmetric
	#generate the nullspace of F
	Mᵀ = nullspace(F')
	M = transpose(Mᵀ)
	#pre-compute matrix
	#This one takes so much time
 	Mat₁=Symmetric(M*Gnu_mat*Mᵀ)
	Mat₂=Symmetric(M*Mᵀ)

	Ydata_M=M*Ydata

    function KrigY(σg,σε,X1,X2,X3)
		G₁₁ = (σg^2) .* Gnu_mat
		Ξ   = [
			G₁₁ .+ σε^2*I(n₁)    	 	F
			F'       zeros(m,m)
		]
		cb = Ξ \ vcat(Ydata, zeros(m))
		c  = cb[1:length(Ydata)]
		b  = cb[length(Ydata)+1:end]
		K= (σg^2) .* Gnu.( sqrt.(norm.( (X1.-Xdata1').^2  + (X2.-Xdata2').^2) )      , ν)
		FF = [t(x=>vx,y=>vy,z=>vz) for (vx,vy,vz) in zip(X1,X1,X3), t in p]
		return K*c.+FF*b
	end

	function loglike(σg,σε)

		Σ₁=σg^2*Mat₁+σε^2*Mat₂
		ch=cholesky(Σ₁)
 		return -1/2*sum( (ch.L\(Ydata_M)).^2)- sum(log.(diag(ch.L)))

	end
	return KrigY, loglike
end

end
