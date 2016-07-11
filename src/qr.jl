import Base.qr

function qr(a::Matrix{AlgebraicNumber})
	n = size(a,1)
	proj(ei,ai) = (dot(ei,ai)/dot(ei,ei))*ei
	u = Array(AlgebraicNumber{BigInt,BigFloat},n,n)
	q = Array(AlgebraicNumber{BigInt,BigFloat},n,n)
	for i = 1 : n
		u0 = a[:,i]
		for j = 1 : (i-1)
			u0 = u0 - proj(u[:,j],a[:,i])
		end
		q[:,i] = u0 / sqrt(dot(u0,u0))
		u[:,i] = u0
	end
	return q
end