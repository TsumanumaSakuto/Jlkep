module Jlkep


using LinearAlgebra
using SPICE

#各種定数
const AU = 1.49597870700e8
const DAY2SEC = 86400
const SEC2DAY = 1.1574074074074073e-5
const DAY2YEAR = 31556952
const DEG2RAD = 0.017453292519943295
const RAD2DEG = 57.29577951308232
const MU_SUN = 1.327124400419393e11
const MU_EARTH = 3.98600435436e5

#軌道6要素形式
struct Par
	a::Float64
	e::Float64
	i::Float64
	Ω::Float64
	ω::Float64
	E::Float64
end

#位置と速度の形式
struct Ic
	r::Vector{Float64}
	v::Vector{Float64}
end

#単位ベクトル
î = [1, 0, 0]
ĵ = [0, 1, 0]
k̂ = [0, 0, 1]

#位置と速度の形式から軌道6要素に変換
function ic2par(r, v, μ)

	r_norm = norm(r)
	ϵ = (v ⋅ v) / 2 - μ / r_norm
	a = -μ / 2ϵ

	h = r × v
	P = (-μ / r_norm) .* r - h × v
	e = norm(P) / μ

	E = 0
	if e == 0
		E = 0
	elseif e < 1
		ecosE = 1 - r_norm / a
		esinE = r ⋅ v / √(μ * a)
		E = atan(esinE, ecosE)
	else
		E = asinh(r ⋅ v / (e * √(-μ * a)))
	end

	Ŵ = h / norm(h)
	N̂ = k̂ × h / norm(k̂ × h)

	i = acos(Ŵ ⋅ k̂)
	Ω = i == 0 ? 0 : atan(Ŵ ⋅ î, -Ŵ ⋅ ĵ)
	ω = e == 0 ? 0 : 2π - acos(N̂ ⋅ (P / norm(P)))
	return Par(a, e, i, Ω, ω, E)
end

function ic2par(IC::Ic, μ)
	ic2par(IC.r, IC.v, μ)
end


#軌道6要素から位置と速度の形式に変換
function par2ic(par::Par, μ)
	r = [0, 0, 0]
	v = [0, 0, 0]

	P̂ = [
		cos(par.ω) * cos(par.Ω) - sin(par.ω) * sin(par.Ω) * cos(par.i),
		cos(par.ω) * sin(par.Ω) + sin(par.ω) * cos(par.Ω) * cos(par.i),
		sin(par.ω) * sin(par.i),
	]
	Q̂ = [
		-sin(par.ω) * cos(par.Ω) - cos(par.ω) * sin(par.Ω) * cos(par.i),
		-sin(par.ω) * sin(par.Ω) + cos(par.ω) * cos(par.Ω) * cos(par.i),
		cos(par.ω) * sin(par.i),
	]
	if par.e < 1
		p = par.a * (1 - par.e^2)
		r_norm = par.a * (1 - par.e * cos(par.E))
		r = par.a * (cos(par.E) - par.e) * P̂ + √(par.a * p) * sin(par.E) * Q̂
		v = -(√(μ * par.a) / r_norm) * sin(par.E) * P̂ + (√(μ * p) / r_norm) * cos(par.E) * Q̂
	else
		p = -par.a * (1 - par.e^2)
		r_norm = -par.a * (par.e * cosh(par.E) - 1)
		r = -par.a * (par.e - cosh(par.E)) * P̂ + √(-par.a * p) * sinh(par.E) * Q̂
		v = -(√(μ * -par.a) / r_norm) * sinh(par.E) * P̂ + (√(μ * p) / r_norm) * cosh(par.E) * Q̂
	end
	return r, v
end

#ケプラー問題を解く
function solve_kepler_equation(a, e, tof, μ, err = 10^-15)

	if e < 1
		#楕円軌道の場合
		M = √(μ / a^3) * tof
		Eₙ = M

		f(E) = E - e * sin(E) - M
		df_dE(E) = 1 - e * cos(E)

		while true
			Eₙ₊₁ = Eₙ - f(Eₙ) / df_dE(Eₙ)
			if abs(Eₙ₊₁ - Eₙ) < err
				break
			end
			Eₙ = Eₙ₊₁
		end
		return Eₙ
	else
		#双曲線軌道の場合
		M = √(μ / (-a)^3) * tof
		Hₙ = M

		F(H) = e * sinh(H) - H - M
		df_dH(H) = e * cosh(H) - 1

		while true
			Hₙ₊₁ = Hₙ - F(Hₙ) / df_dH(Hₙ)
			if abs(Hₙ₊₁ - Hₙ) < err
				break
			end
			Hₙ = Hₙ₊₁
		end
		return Hₙ
	end
end

#軌道伝播を行う関数．初期位置と初期速度と時刻から現在の速度と位置を求める

#一度軌道6要素に変換してから元に戻して伝播を行う方法．
#ただiが0のときや離心率が0の時など特殊な場合うまくいかない可能性がある
"""
function propagate_lagrangian(r₀, v₀, tof, μ)
	par = ic2par(r₀, v₀, μ)
	par.E = solve_kepler_equation(par.a, par.e, tof, μ)
	return par2ic(par, μ)
end
"""

#ラグランジュの係数を求めてから軌道伝播を行う方法．
#ただeが1に近い時などうまくいかない可能性がある
#動作未確認
function propagate_lagrangian(r₀, v₀, tof, μ, err = 10^-15)

	a = 1 / (1 / norm(r₀) - v₀ ⋅ v₀ / μ)

	if a > 0
		M = √(μ / a^3) * tof
		ΔEₙ = M

		F(ΔE) = ΔE - (1 - norm(r₀) / a) * sin(ΔE) + r₀ ⋅ v₀ * (1 - cos(ΔE)) / √(μ * a) - M
		dF_dΔE(ΔE) = 1 - (1 - norm(r₀) / a) * cos(ΔE) + r₀ ⋅ v₀ * sin(ΔE) / √(μ * a)

		while true
			ΔEₙ₊₁ = ΔEₙ - F(ΔEₙ) / dF_dΔE(ΔEₙ)
			if abs(ΔEₙ₊₁ - ΔEₙ) < err
				break
			end
			ΔEₙ = ΔEₙ₊₁
		end
		ΔE = ΔEₙ
		f = 1 - a * (1 - cos(ΔE)) / norm(r₀)
		g = tof - √(a^3 / μ) * (ΔE - sin(ΔE))

		r = f * r₀ + g * v₀

		ḟ = -√(μa) * sin(ΔE) / (norm(r) * norm(r₀))
		ġ = 1 - a * (1 - cos(ΔE)) / norm(r)

		v = ḟ * r₀ + ġ * v₀
		return r, v
	else
		M = √(μ / (-a)^3) * tof
		ΔHₙ = M

		F(ΔH) = (1 + norm(r₀) / (-a)) * sinh(ΔH) + r₀ ⋅ v₀ * (cosh(ΔH) - 1) / √(μ * (-a)) - ΔH - M
		dF_dΔE(ΔE) = (1 + norm(r₀) / (-a)) * cosh(ΔH) + r₀ ⋅ v₀ * sinh(ΔH) / √(μ * (-a)) - 1

		while true
			ΔHₙ₊₁ = ΔHₙ - F(ΔHₙ) / dF_dΔE(ΔHₙ)
			if abs(ΔHₙ₊₁ - ΔHₙ) < err
				break
			end
			ΔHₙ = ΔHₙ₊₁
		end
		ΔH = ΔHₙ
		f = 1 + (-a) * (1 - cosh(ΔH)) / norm(r₀)
		g = tof - √((-a)^3 / μ) * (ΔH - sinh(ΔH))

		r = f * r₀ + g * v₀

		ḟ = -√(μ * (-a)) * sinh(ΔH) / (norm(r) * norm(r₀))
		ġ = 1 + (-a) * (1 - cosh(ΔH)) / norm(r)

		v = ḟ * r₀ + ġ * v₀
		return r, v
	end
end


end
