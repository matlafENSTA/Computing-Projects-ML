### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 3fd3bb54-095c-47aa-9a94-9532b7b6e09b
using LinearAlgebra, Plots, PlutoUI, Test

# ╔═╡ 606540b0-9846-11eb-2a09-8112c5081854
md"""
# ANN 203
"""

# ╔═╡ 28280ba3-789f-40ec-b731-cbc43334b839
md"""!!! danger "Submission instructions"
	**Due date:** 5pm on April 4, 2024.
	
	**Instructions:** Submit your homework on [ecampus](https://ecampus.paris-saclay.fr). 
"""

# ╔═╡ 6ea0197b-5af3-4116-b214-a27f28508c33
md"""## Introduction

This homework is intended to give you a *hands-on* experience on coding some of the fundamental numerical linear algebra algorithms that we have covered in the lectures so far. You are encouraged to consult the course material available [here](https://perso.ensta-paris.fr/~mbonnet/aln.pdf). Feel free to reference the formulas in the course material if needed. 

You will have *one week* to submit the homework on *ecampus*. Unless a valid excuse is presented, **late homeworks will not be accepted.** 

Let's get started. 
"""

# ╔═╡ 8c26975d-ec0c-423c-9644-daf031c0aabb
md""" ## Part 1: QR factorization

Throughout this section we will assume that $A \in \mathbb{C}^{m\times n}$ with $m \geq n$, and $\text{rank}(A) = n$ (i.e. $A$ has full rank). This is a simplyfying assumption which is not necessarily needed in many of the algorithms (such as QR).  

The main goal of this section is to write a program that can solve the following least squares problem:
```math
	\text{least squares} \rightarrow \underset{x \in \mathbb{C}^n}{\mathrm{argmin}}|| Ax - b ||_2
```
Since we assume that $A$ has full rank, this problem has a unique solution. We will proceed to solve this problem by performing a QR factorization of $A$ through Householder reflectors, as described next. """

# ╔═╡ c707bd38-32f2-4adb-bdaf-8c8325d40ab9
md"""### Householder reflectors
Recall the definition of the Householder reflector:
```math
F(u) = I - 2 \frac{uu^*}{u^* u},
```
where $I$ is the identity matrix and $u$ is a complex vector. The notation $(\cdot)^*$ is used for the adjoint. The first task to code precisely this:
"""

# ╔═╡ 989434b2-0f36-4998-b213-5ba6d2ce7166
md"""!!! tip 
	You can use `I` for representing an identity matrix in *Julia*. It is similar to *Matlab*'s `eye`, but is it is *lazy*, meaning that not matrix is actually formed. Check also `diagm` for a more general method. If you are running this notebook interactively through *Pluto.jl*, dont forget to check the live docs for more information!
"""

# ╔═╡ 06ddb55b-e7a7-4ee8-91dd-171d5665df5e
function householder_matrix(v)
	return I
end

# ╔═╡ 5a6887a8-e651-4e4b-96c6-d1d7058cd26d
md"""
Now, at each step of the *QR* factorization, we need to compute a Householder vectors $u$ given a vector $x$ so that 
```math 
F(u)x = \textrm{sign}(x_1)||x||e_1
```
where $x_1$ is the first entry of $x$ and $e_1$ is the cartesian basis vector. This is your second task:
"""

# ╔═╡ 1186d2ae-f51e-4237-a14f-e4275b5b304c
function householder_vector(x)
    u = copy(x)
    # complete the code
    return u
end

# ╔═╡ 21796460-d43b-4554-adbc-013f58cf92c3
# write some tests here to make sure your function is correct in this cell

# ╔═╡ ef865aba-8f80-4f38-9560-feecdb9013a7
md"""### *Naive* QR algorithm
We can now code a version of the Householder QR factorization similar to [Algorithm 3.1](https://perso.ensta-paris.fr/~mbonnet/aln.pdf) on your notes. One difference will be that instead of simply reducing the matrix $A$ to an upper triangular form using orthogonal triangularization, we will actually form the matrix $Q$ in the process. In general this is inneficient, but we will do that anyways to learn the steps. 

Recall that the main idea is to gradually convert $A$ into an upper triangular matrix $R$ as follows:
```math
\underbrace{F(u_n)F(u_{n-1}) \ldots F(u_2)F(u_1)}_{Q^*} A = R.
```
Remember: the matrices $F(u_k)$ act only on the $k$ through $n$ columns of the matrix to its right. On your lectures notes that was achieved by some *zero padding* of the vectors $u_k$. Here we will only compute the non-zero part of these matrices. 
"""

# ╔═╡ 65523ce4-ebb8-4fb2-9869-61890e0b7e2d
function householder_qr(A)
	R   = copy(A)
    m,n = size(A)
	kmax = m>n ? n : m-1
    Qadj  = diagm(0=>ones(eltype(A),m)) # m×m identity matrix
    for k in 1:kmax
		# modify Qadj and R columns by column
    end
    Q = adjoint(Qadj)
    return Q,R
end

# ╔═╡ 665919fd-592f-494a-9d17-adde3434a7e0
md"""
If you got this far you have succeeded in implementing a $QR$ factorization in *Julia*! 
"""

# ╔═╡ 7c40b5d3-8afd-4648-91ef-ec6235ec376e
md"""### A compact QR algorithm
Forming the matrix $Q$ requires some extra computation (and extra memory), and having $Q$ explicitly formed is often not necessary. In what follows we will dive a little into a more efficient *compact* representation of the *QR* factorization. The idea is to use the upper triangular part of the input matrix $A$ to store $R$, and the lower triangular part (diagonal excluded) to store the Householder vectors $u_k$. We will normalize $u_k$ so that its first value is always one, meaning only $m-k$ values are needed to define it instead of $m-k+1$ values. This *trick* means $u_k$ will fit in the lower triangular part (diagonal excluded) of the $k$-th column of $A$, and therefore we should be able to perform our QR factorization using *only the memory* allocated in $A$. This will require some significant rewrite. In particular

* The `householder_vector` function should normalize the vector so that the first entry is `1`.
* No Householder matrices should ever be stored in memory. This means we need to write code to *apply* $H(u)$ to a matrix given the vector $u$ that defines it.  
* Some care has to be taken regarding the order of the operations so that we do not overwrite memory that is still needed for other computations. 

The function `householder_vector_normalized!(x)` shown below modifies `x[2:end]` to store the `2:end` entries of the Householder vector. Note that the first entry is not modified since it is (implicitly) one. This memory space will be reserved to the diagonal part of $R$.
"""

# ╔═╡ e24e5718-1a1a-47bc-8f25-3268ce9a5826
function householder_vector_normalized!(x)
    # norm of x
    d2 = x' * x
    d  = sqrt(d2)
    # normalized reflection vector
    scale = 1/(x[1] + sign(x[1])*d)
	m     = length(x)
	for k in 2:m
		x[k] *= scale
	end
    return x
end

# ╔═╡ 13c5b2d6-7f30-4afd-8a08-3dfe4f2db184
md"""
Since only the action of Householder matrices $H(u)$ on a matrix/vector is needed, it makes sense to avoid creating a dense matrix to store it. The `apply_hh_matrix` function below uses the special structure of $H(u)$ to efficiently computes $H(u)A$ where $u \in \mathbb{C}^m$ and $A \in \mathbb{C}^{m\times n}$: 
"""

# ╔═╡ 355bd86a-a278-4ade-96a7-2a18e0c8f7db
function apply_hh_matrix!(A,u)
	m,n = size(A)
	@assert length(u) == m
	beta = 2/(u'*u)
	for j in 1:n
		col   = @views A[:,j]
		coeff = u'*col*beta
		for i in 1:m
			col[i] -= coeff*u[i]
		end
	end
	return A
end

# ╔═╡ e24068de-07a2-41df-a686-01b3e21ca6b0
md"""!!! tip
	In `julia`, *slicing* an array (e.g. `x[5:10]`) creates a new array with the 	sliced entries. This means memory is allocated to store the new object. The `@views` macro used in the `apply_hh_matrix!` function above avoids such allocations by creating a `view` of the array which uses the same memory. Look for `view` in the *Live docs* for more information.
"""

# ╔═╡ 4bbfd29d-832b-4e17-8f59-25463466fa37
md"It is your turn now to explain the code above:"

# ╔═╡ b7dc76de-d6df-414e-bc76-fa9c5e2af727
function apply_hh_matrix_normalized!(A,u)
	# a good starting point is the code for apply_hh_matrix!
	return A
end

# ╔═╡ ae805230-574e-4818-bdbc-1dda19fa33b0
function householder_qr_compact!(A)
    m,n = size(A)
    for k=1:n
		# iteratively modify A in place
    end
    return A
end

# ╔═╡ a0af5567-02c6-4c3a-a1f2-649c622d1319
let
	m,n = 50,50
	A = rand(ComplexF64,m,n)
	t1 = @elapsed householder_qr(A)
	t2 = @elapsed householder_qr_compact!(A)
	t3 = @elapsed qr!(A)
	m1 = @allocated householder_qr(A)
	m2 = @allocated householder_qr_compact!(A)
	m3 = @allocated qr!(A)
	md"""You case you got this far, the table below will give you a rough idea of execution time of the *naive* and *compact* implementations of your algorithm for a matrix of $A \in \mathbb{C}^{m\times n}$ with m=$m and n = $n. We will also compare to the default `qr` algorithm provided by the `LinearAlgebra` package. The dynamically allocated memory is also shown. These numbers only make sense **if your algorithm is working**. 
	
Version | Time (s) | Memory (bytes)
:------------ | :-------------: | :----------:
naive   | $t1 | $m1
compact | $t2 | $m2
`LinearAlgebra`  | $t3 | $m3
	"""
end

# ╔═╡ a1cfd851-c006-4506-be18-7571927fb5ef
md"""
You are now (almost) ready to solve the original least squares problem!
"""

# ╔═╡ f89404e4-9287-4014-bd2d-93f6833a1985
md"""## Part 2: Least squares
"""

# ╔═╡ 8efa2b75-f85a-4be9-85cd-6cccd193bded
md"""Now that everthing is in place, solve the least squares problem:
"""

# ╔═╡ 2eb53b12-d20d-411c-9da1-96a0b6ecac39
# function provided to solve Rx = b 
function solve_uppertriangular(R,b)
    m,n = size(R)
    @assert m==n
    @assert n==length(b)
    x = copy(b)
    for i in n:-1:1
        for j in i+1:n
            x[i] -= R[i,j]*x[j]
        end
        x[i] /= R[i,i]
    end
    return x
end

# ╔═╡ 7c2e9206-8403-4eab-a39e-6625945dd0c1
function leastsquares_qr(A,b)
	# compute the solution x and return it
end

# ╔═╡ e4350897-af1a-4081-93aa-2eb24b3b5a0a
function leastsquares_qr_compact!(A,x)
	return nothing
end

# ╔═╡ c50531af-4176-40a9-9812-1fae66ffb9f0
md"""
## Part 3: Polynomial approximation

In this third and last part we will focus on *polynomial interpolation* and its connection to *numerical linear algebra*. More precisely, we will investigate the following problem:


!!! danger "Polynomial interpolation"
	Given a continuous function $f : [-1,1] \to \mathbb{R}$, and a set of $n+1$ interpolation points 
	```math
		-1 \leq x_0 < \ldots < x_n \leq 1,
	```
	find $p_n \in \mathbb{P}_n = \text{span}\{x^i\}_{i=0}^n$ such that $p_n(x_i) = f(x_i)$ for $0 \leq i \leq n$. 
"""

# ╔═╡ a52638e9-3115-4951-917a-a7b325ea94c1
function monomial_interpolation(x,y)
	n = length(x)
	A = zeros(n,n) # FIXME: compute the vandermond matrix
	c = ones(n) # FIXME: find the monomial coefficients
	p = (x) -> evalpoly(x,c) # this creates a polynomial. See docs for evalpoly
	return p, cond(A), norm(c,Inf)
end

# ╔═╡ d7e82187-d508-4447-a151-314f85749517
md"""
The choice of interpolation points ``\mathcal{I}_n = \{ x_i \}_{i=0}^n`` plays an important role in the quality of the approximation (assuming exact arithmetic), while the choice of basis for ``\mathbb{P}_k`` has affects the stability of the interpolation method (i.e. how well can you actually find the unique polynomial you seek in the presence of rounding errors). 

Defining ``\Pi : \mathbb{R}^{n+1} \to \mathbb{P}_n`` as the map from the ``n+1`` function values to the resulting (unique) polynomial in ``\mathbb{P}_n``, we may define its operator norm as
```math
\Lambda(\mathcal{I}_n) = \sup_{y \in \mathbb{R}^{n+1}} \frac{|| \Pi(y)||}{||y||}
```
Often the $L^\infty$ norm is used here, but that is not very important for our discussion. ``\Lambda(\mathcal{I}_n)`` is commonly refered to as its **Lebesgue constant** for ``\mathcal{I}_n``. It is a simple exercise to show that
```math
|| f - p_n || \leq (1 + \Lambda(\mathcal{I_n})) \inf_{q_n \in \mathbb{P}_n}|| f - q_n ||, 
```
and therefore a "small" Lebesgue constant means we are close to the best possible choice. 

As it turns out, some points are significantly better than others for polynomial interpolation. Perhaps counterintuitively, equispaced points are very bad in that their Lebesgue constant diverges exponentially as ``n \to \infty``. Among the *good* choice of points, Chebyshev nodes are particularly useful and easy to implement, as we shall see next.
"""

# ╔═╡ d619f2e2-31ad-4a56-9225-c99a0b635331
function chebyshev_points(n) 
	x = zeros(n) # FIXME
	return x
end

# ╔═╡ 4eef11e8-b650-423d-8288-f11f84c7bfc9
md"""!!! note
	A common feature with all *good* set of interpolation points is that they accumulate values near the endpoints, as you can see below:
"""

# ╔═╡ 2a350f62-caca-4139-96b1-06186b043755
let
	n = 20
	scatter([-cos((k-1)/(n-1) * π) for k in 1:n], zeros(n); label="chebyshev")
	scatter!(collect(range(-1,1,n)), zeros(n); label="uniform",m=:cross)
	ylims!(-0.5,0.5)
end

# ╔═╡ bd2e2f6b-0b66-4907-abc3-fac25f50d89f
md"""!!! note
	In case your code for `monomial_interpolation` does not work, the behavior that would be observed for a correctly implemented `monomial_interpolation` is the following: the interpolation error initially decrease, but will plateau way above the expected machine precision. This is the *puzzling* behaviour you must try to explain.  
"""

# ╔═╡ cc4ee8ce-68fe-43b6-92ae-e59b05f87eb9
let
	f = (x) -> atan(4x) # an analytic function
	xtest = -1:0.01:1
	println("Estimated interpolation error:")
	for n in [2^i for i in 1:10]
		x = chebyshev_points(n)
		y = [f(x) for x in x]
		p, κ, cmax = monomial_interpolation(x,y)
		er = norm(f.(xtest)-p.(xtest), Inf)
		println("n = $n") 
		println("|-- error      = $er")
		println("|-- condition  = $κ")
		println("|-- coefs_norm = $cmax")
	end
end

# ╔═╡ d19e8ee2-f546-4667-b38a-6af8006158dc
md"""
The failure that we (hopefully) observed above is *not* attributed to the fact that the error $|| f - p_n ||_{\infty}$ fails to converge as $n \to \infty$, but rather to the fact that we have failed to compute $p_n$ in inexact arithmetic. 

In what follows we use a **Lagrange basis** instead, with a stable barycentric formula for the polynomial evaluation. In particular, the `lagrange_interpolation` serves as an alternative to the `monomial_interpolation` method you implented above: 
"""

# ╔═╡ 109a5754-511c-4df1-b79b-a145e35504f8
function barycentric_lagrange_weights(x)
	n = length(x) - 1
	w = similar(x)
	w[1] = 1.0
	for j in 1:n
		for k in 0:(j - 1)
			w[k + 1] = (x[k + 1] - x[j + 1]) * w[k + 1]
		end
		w[j + 1] = prod(0:(j - 1)) do k
			return x[j + 1] - x[k + 1]
		end
	end
	for j in 0:n
		w[j + 1] = 1 / w[j + 1]
	end
	return w
end

# ╔═╡ 671f56ef-689b-4300-8ca0-5fc9e32d2146
function lagrange_basis(nodes)
	n = length(nodes)
	w = barycentric_lagrange_weights(nodes)
	l = x -> prod(xi->x-xi,nodes)
    (x) -> map(1:n) do i
		if x == nodes[i]
			1.0
		else
			l(x)*w[i]/(x-nodes[i])
		end
	end
end

# ╔═╡ 14427835-cace-4a5f-94c9-bac9463cc14c
function lagrange_interpolation(x,y)
	b = lagrange_basis(x)
	p = (x) -> sum(dot(b(x),y))
	return p
end

# ╔═╡ 8fde2a0a-9dcd-4b81-b7ae-4249587ebaac
md"""
Here is a figure that shows how the interpolation error changes as a function of $n$ when using two ways of computing the interpolating polynomial $p_n$ given a set of Chebyshev nodes. ``\kappa(V)`` denotes the condition number of the Vandermond matrix.

![](https://github.com/maltezfaria/ANN203/blob/main/convergence.png?raw=true)
"""

# ╔═╡ 1bf61994-f5d6-4ec6-b088-4f60a44953b3
md"""!!! note
	If your code for `monomial_interpolation` works, you can comment the cell below to generate the figure above yourself.
"""

# ╔═╡ 1425d9ea-d1ad-4b86-ab9d-ecea3a73a3df
# let
# 	f = (x) -> atan(4x)
# 	#f = (x) -> x + x^2
# 	xtest = -1:0.01:1
# 	mm = 4:1:200
# 	ee_mon = []
# 	ee_lag = []
# 	cc = [] 
# 	for m in mm
# 		x = chebyshev_points(m)
# 		y = [f(x) for x in x]
# 		p, κ, cmax = try 
# 			monomial_interpolation(x,y)
# 		catch
# 			@warn "monomial iterpolation not working" maxlogs=1
# 			0*x
# 		end
# 		l = lagrange_interpolation(x,y)
# 		er_mon = norm(f.(xtest)-p.(xtest),Inf)
# 		er_lag = norm(f.(xtest) - l.(xtest), Inf)
# 		push!(ee_mon,er_mon)	
# 		push!(ee_lag,er_lag)
# 		push!(cc,κ)
# 	end
# 	plot(mm,ee_mon, yscale=:log10, label="monomial", lw=4, legend=:topleft, ylabel="maxerror", ylims=(eps(),1), xlabel="interpolation order")
# 	plot!(mm,ee_lag, label="lagrange", lw=4)
# 	#plot!(mm,eps()*ones(length(cc)), label="lagrange", lw=4)
# 	plot!(twinx(),cc, lc=:black, lw=4, yscale=:log10, label="κ(V)", ylabel="condition")
# 	plot!(twinx(),1/eps()*ones(length(cc)), lc=:red, lw=4, yscale=:log10, label="", ls=:dash)
# end

# ╔═╡ 254bfbc4-117b-4518-a188-c19a7c123035
md"""### Utility functions
"""

# ╔═╡ 57160c99-cccc-4aca-99fd-df9356da76a1
PlutoUI.TableOfContents(;depth=2)

# ╔═╡ 363968d9-23c9-45e7-b47d-5d7e8ad5648f
question(text,number="") = Markdown.MD(Markdown.Admonition("tip", "Question $number", [text]))

# ╔═╡ 2cded3aa-8661-4477-bcc8-ac375e01622c
question(md"Complete the `householder_matrix` function below to compute the Householder reflector.",1)

# ╔═╡ 4720b6c3-d9bf-4a46-9d38-2cee67e35142
question(md"""Complete the `householder_vector` function below to implement the  	  *Householder vector* $u$ given an input vector $x$. Check that the property 
	```math 
		F(u)x = \textrm{sign}(x_1)e_1||x||_2
	```
	is satisfied, where $e_1 = [1,0,\ldots,0]$
	""",2
)

# ╔═╡ 9cfa6c75-cb75-470d-b037-2f01c2f22ae9
question(md"Using your functions `householder_vector` and `householder_matrix`, complete the code below to implement a modified version of Algorithm 3.1 of the lecture notes where the matrices Q and R are explicitly constructed. You should return the *full QR*, so that $Q \in \mathbb{C}^{m\times m}$ and $R \in \mathbb{C}^{m \times n}$",3)

# ╔═╡ ee32ef8b-3b5f-4bd8-bafc-588166899782
question(md"Describe in words and formulas the steps in the code for `apply_hh_matrix!` above. Feel free to include snippets of the code if you find that helps the explanation.",4)

# ╔═╡ 0f634e30-4f59-48d4-9843-8bbb1e244db5
question(md"""Modify `apply_hh_matrix!` so that it implicitly assumes `u[1]==1`. We will call the new function `apply_hh_matrix_normalized`. Note that this function should **not** use the value of `u[1]`.
""",5)

# ╔═╡ 5ed2c935-22d6-4e3c-ab2a-a52dd8ee462b
question(md"""
Complete the code below to perform a compact QR factorization in place.
""",6)

# ╔═╡ 4e5aefd4-b965-40d6-ae04-cda4f36945d2
question(md"Can you make your $QR$ code even faster?",
"(bonus)")

# ╔═╡ dbf408ba-bf6e-440c-8855-3f7b7f254d6e
question(md"""Explain how you may reformulate the least squares problem 
```math
\begin{align}
\underset{x \in \mathbb{R}^n}{\mathrm{argmin}}||Ax - b||_2
\end{align}
```
using the QR factorization of A. Then complete the function `least_squares_qr` below to calculate $x$ given $A$ and $b$. You can use the function `solve_uppertriangular(R,b)`, written below, as well as your implementation of `householder_qr`, solve the least squares problem.
	""",7)

# ╔═╡ 6217aa0d-46e1-4dd6-a1e5-5eb2696757e4
question(md"""Explain how to solve the least squares problem using the `householder_qr_compact!` factorization, where you do not explicitly form `Q`. What method(s) do you need to implement? 
""",8)

# ╔═╡ 13aff5d5-e813-40de-adc1-21fa8d07c812
question(md"""Implement `leastsquares_qr_compact!(A,x)` which uses the in-place `QR` factorization of Part 1 to solve the least-squares problem without ever forming the matrix `Q`. 
""","(bonus)")

# ╔═╡ 2143cef3-f43c-422d-b738-1aadbd72660a
question(md"""
Prove that in exact arithmetic the **interpolation problem** above has a unique solution $p_n \in \mathbb{P}_n$ by considering a *monomial basis* for $\mathbb{P}_k$ and writting down the resulting linear system for the coefficients. 
""", 9)

# ╔═╡ 4292d332-ac6e-4b40-96e6-dbf086c1118d
question(md"""
Based on your previous answer, complete the function below to generate the interpolating polynomial on a monomial basis given a set of nodes `x` and a set of values `y`. Because we will later need it, this function will also return the condition number of the interpolation/Vandermond matrix, and the maximum norm of the coefficients.
""",10)

# ╔═╡ 954fbbeb-8ca5-4e3d-96de-1b4a61c645fe
question(md"""
Complete the function below to return the Chebyshev nodes of the second kind, given by the projection onto the $x$ axis of nodes equispaced on a half-circle of radius 1. In ascending order, the set of nodes are given by the following formula:
```math
	\mathcal{I}_{n-1} = \left\{-\cos\left(\frac{k}{n}\pi \right) \right\}_{k=0}^{n-1}
``` 
""")

# ╔═╡ c3472420-ad6f-4146-89bf-3adcdd57200a
question(md"""
Interporlating a smooth function on Chebyshev nodes should converge exponentially fast as $n \to \infty$. Explain the behaviour of the code below given what you have learned about linear algebra. Does the error decay to machine precision? Is it even true that ``p(x_i) = f(x_i)`` for the interpolation nodes? You can create *code cells* to explore e.g. how the residue behaves for your interpolation problem.
""",11)

# ╔═╡ 10116502-8a73-45ca-ab72-8e5b6ce6d115
question(md"""
Given what you know about linear algebra, assuming a *backward stable* method is used to solve the linear system ``Vc = b``, where ``V`` is the Vandermond matrix and ``c`` the sought polynomial coefficients, what can you say about the residue ``||V\tilde{c} - b||``, where $\tilde{c}$ denotes the output of the algorithm (i.e. the approximate solution)?
""", 12)

# ╔═╡ b061d7f4-1d38-4f93-abce-8be49363204c
question(md"""
Given what you have learned so far, what general conclusions can you draw from stability and convergence of polynomial interpolation on Chebyshev nodes? Is there a difference between theory and practice?
""", 13)

# ╔═╡ c65c438b-d97f-4857-9897-20180cc53d25
answer(text=md"Write your answer here!
") = Markdown.MD(Markdown.Admonition("warning", "Answer", [text]))

# ╔═╡ 21f916fe-21b9-43c1-b9d1-a8e541d24e43
answer(md""" Your answer goes here.
""")

# ╔═╡ f31d5edd-722e-4c55-9ef3-190b4996fd29
answer(md"""
Write your explanation here (and don't forget to complete the code below).
""")

# ╔═╡ 165cfb6b-63f6-4857-849b-559ba9221768
answer(md"""
Write your answer here!
""")

# ╔═╡ c6bca558-46de-40c1-b54d-7c61812a62c7
answer()

# ╔═╡ 0348cadb-80c0-425a-b2b3-4266fee723ec
answer()

# ╔═╡ f0c70407-13be-449c-ac0c-fce644caaf2f
answer()

# ╔═╡ e4b5d949-c047-48b6-8f51-bf17bd152f5f
answer()

# ╔═╡ 78e2208d-f9f8-44f7-806f-4590db85fe1a
hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]))

# ╔═╡ 5fc13f1f-5818-4ebd-8531-fa5900c47bc9
hint(md"Think about the operation $Q^*b$ and how to implement that when only the normalized Householder vectors have been stored.")

# ╔═╡ a7aa6fe4-3ee2-4dff-a792-4289d4fca83c
keep_working(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("danger", "Keep working on it!", [text]))

# ╔═╡ e25e883e-e0f9-4848-adfd-24f155ecf902
correct(text=md"Great! You got the right answer! Let's move on to the next section.") = Markdown.MD(Markdown.Admonition("note", "Got it!", [text]))

# ╔═╡ 21f882ba-fa33-4695-9324-b058ba64d258
note(text) = Markdown.MD(Markdown.Admonition("note", "Note", [text]))

# ╔═╡ 8c0961b4-4982-4355-9f7a-6aa37e342f3c
note(md"The cell below will add the dependencies that you need for this notebook. This may take some time the first time you run it since some packages will be downloade and precompiled.")

# ╔═╡ 92791a5f-02f8-461e-9bb1-4e1077f24be9
note(md"""In *Julia*, functions whose name end with `!` typically mutate their arguments. This is merely a convention, but sticking to it is considered good practice. """)

# ╔═╡ f5506679-584a-4f57-93c5-5b71f43b0896
note(md"For more serious benchmarks, you should the `BenchmarkTools` package.")

# ╔═╡ 74097a18-b4fe-4ec8-b47c-b279fce28598
note(md"""
In case your `householder_qr` function does not work, you may use `LinearAlgebra.qr` for perfoming the $QR$ factorization. Check its documentation for more details on how to use it.
""")

# ╔═╡ 167a446c-24b7-414c-9bdf-bae7c90df2ab
warning(text) = Markdown.MD(Markdown.Admonition("warning", "Warning", [text]))

# ╔═╡ 7e94c5d7-ba40-4850-9602-18fe63381ee8
function check_answer(f)
	try
		t = f()
		correct()
	catch E
		keep_working()
	end
end

# ╔═╡ bb7b2a43-4e93-4bfb-a41f-65ff6376babd
check_answer() do
	@testset "Householder matrix" begin
		v = rand(10)
		F = householder_matrix(v)
		@test F*v ≈ -v # reflection
		@test F*F' ≈ F' * F ≈ I # unitary
	end
end

# ╔═╡ 15410cfb-db47-48e0-a7e0-49a061aaaec5
check_answer() do 
	A    = rand(5,5) .- 0.5 # a random matrix for testing
	a1   = A[:,1]
	u1   = householder_vector(a1)
	F    = householder_matrix(u1)
	A1   = F*A
	@testset "Householder vector" begin
		@test A1[1,1] ≈ -sign(a1[1])*norm(a1)
		@test norm(A1[2:end,1],Inf) < 1e-10
	end
end


# ╔═╡ 9fc30c4d-e5c2-4359-9bec-e93ca5876763
check_answer() do 
	@testset "Householder qr" begin
		for (m,n) in [(10,5),(10,10)]
			A   = rand(ComplexF64,m,n)
			Q,R = householder_qr(A)
			@test size(Q) == (m,m)
			@test size(R) == (m,n)
			# test ortogonality
			@test Q'*Q ≈ I
			# test factorization
			@test Q*R  ≈ A
			# test that R is lower triangular
			@test norm(tril(R,-1),Inf) < 1e-15
		end
	end
end

# ╔═╡ 297f8790-3973-462e-b14a-bff39f399e8e
check_answer() do 
	A = rand(10,5)
	x = rand(10)
	exact = (I - 2*x*x'/(x'*x))*A
	apply_hh_matrix!(A,x)
	@testset "Appling householder matrix" begin
		@test exact ≈ A
		mem = @allocated apply_hh_matrix!(A,x)
		@test mem == 0
	end
end

# ╔═╡ e15fdd6e-ed3e-49bf-a3a7-7530cda0d479
check_answer() do 
	A    = rand(10,5)
	x    = rand(10)
	x[1] = 1
	ex = (I - 2*x*x'/(x'*x))*A
	x[1] = rand()
	apply_hh_matrix_normalized!(A,x)
	@testset "Appling householder matrix normalized" begin
		@test ex ≈ A
		mem = @allocated apply_hh_matrix_normalized!(A,x)
		@test mem == 0
	end
end

# ╔═╡ 2e865221-74b1-4023-962c-77203686cf03
check_answer() do 
	A = rand(10,5)
	Q,R = householder_qr(A)
	householder_qr_compact!(A)
	@testset "Compact QR" begin
		@test triu(A) ≈ R
		mem = @allocated householder_qr_compact!(A)
		@test mem == 0
	end
end

# ╔═╡ c727e963-6afb-41f5-9de8-38026ab126e0
check_answer() do
	@testset "Least squares QR" begin
		m,n = 10,5
		A = rand(ComplexF64,10,5)
		b = rand(10)
		x = leastsquares_qr(A,b)
		xe = A\b
		@test norm(x-xe,Inf) < 1e-10
	end
end

# ╔═╡ 6b571796-a8c5-40fe-838e-fbbadd2469ad
check_answer() do
	@testset "Monomial interpolation" begin
		n = 10
		x = [-cos((2k-1)/2n*π) for k in 1:n]
		f = x -> x^10 - x^9 + x^3
		y = f.(x)
		p, _ = monomial_interpolation(x,y)
		@test norm(p.(x)-f.(x),Inf) < 1e-10
	end
end

# ╔═╡ 7dde9e2b-462f-4f52-be7d-a1e50300efa1
check_answer() do
	@testset "Chebyshev nodes" begin
		@test chebyshev_points(5)[1] ≈ - 1
		@test chebyshev_points(5)[end] ≈ 1
		@test chebyshev_points(5) ≈ [-1, -sqrt(2)/2, 0, sqrt(2)/2, 1]
	end
end

# ╔═╡ 8b8a87c1-871d-4fab-bf1a-6b3d8bc9a0fd
function show_tests(f)
	with_terminal() do
		try
			f()
		catch
		end
	end
end

# ╔═╡ 2f9ccd3c-14a1-486e-8b6c-707deda8b4bc
show_tests() do 
	x     = rand(ComplexF64,5) .- 0.5 # a random matrix for testing
	v     = copy(x)
	v     = householder_vector_normalized!(v)
	v[1]  = 1 # explicitly store 1 here for testing
	H     = I - 2*v*v'/ (v'*v)
	p     = H*x
	@testset "Householder vector normalized" begin
			@test norm(p[2:end,1],Inf) < 1e-10
			@test p[1] ≈ -sign(x[1])*norm(x)
			mem = @allocated householder_vector_normalized!(v)
			@test mem == 0
	end
end

# ╔═╡ 9ad83c33-f843-4d63-b0f3-1a0401068c80
show_tests() do 
	try
		@testset "Upper triangular solver" begin
			A = rand(10,10)
			b = rand(10)
			R = triu(A)
			xe = R\b
			x = solve_uppertriangular(R,b)
			@test norm(x - xe,Inf) < 1e-10
		end
	catch
		@warn "Exception encoutered"
	end
end

# ╔═╡ 0ef60b26-4769-45d8-b1a3-e2a87bb0d8a2
show_tests() do
	@testset "Lagrange interpolation" begin
		x = rand(5)
		p = lagrange_basis(x)
		@test p(x[1]) ≈ [1, 0, 0, 0, 0]
		@test p(x[2]) ≈ [0, 1, 0, 0, 0]
		@test p(x[5]) ≈ [0, 0, 0, 0, 1]
	end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
Plots = "~1.40.2"
PlutoUI = "~0.7.58"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "31ecda1a2ba2ff286361f1d3509d9eab0cbf362a"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0f748c81756f2e5e6854298f11ad8b2dfae6911a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a4c43f59baa34011e303e76f5c8c91bf58415aaf"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "6cbbd4d241d7e6579ab354737f4dd95ca43946e1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "0f4b5d62a88d8f59003e43c25a8a90de9eb76317"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.18"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "3437ade7073682993e092ca570ad68a2aba26983"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.3"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a96d5c713e6aa28c242b0d25c1347e258d6541ab"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.3+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "359a1ba2e320790ddbe4ee8b4d54a305c0ea2aff"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.0+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "8e59b47b9dc525b70550ca082ce85bcd7f5477cd"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.5"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3336abae9a713d2210bb57ab484b1e065edd7d23"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cad560042a7cc108f5a4c24ea1431a9221f22c1b"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.2"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dae976433497a2f841baadea93d27e68f1a12a97"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.39.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0a04a1318df1bf510beb2562cf90fb0c386f58c4"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.39.3+1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "af81a32750ebc831ee28bdaaba6e1067decef51e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60e3045590bd104a16fefb12836c00c0ef8c7f8c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "3c403c6590dd93b36752634115e20137e79ab4df"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.2"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "14389d51751169994b2e1317d5c72f7dc4f21045"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.6"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "07e470dabc5a6a4254ffebc29a1b3fc01464e105"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "31c421e5516a6248dfb22c194519e37effbf1f30"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─606540b0-9846-11eb-2a09-8112c5081854
# ╟─8c0961b4-4982-4355-9f7a-6aa37e342f3c
# ╠═3fd3bb54-095c-47aa-9a94-9532b7b6e09b
# ╟─28280ba3-789f-40ec-b731-cbc43334b839
# ╟─6ea0197b-5af3-4116-b214-a27f28508c33
# ╟─8c26975d-ec0c-423c-9644-daf031c0aabb
# ╟─c707bd38-32f2-4adb-bdaf-8c8325d40ab9
# ╟─2cded3aa-8661-4477-bcc8-ac375e01622c
# ╟─989434b2-0f36-4998-b213-5ba6d2ce7166
# ╠═06ddb55b-e7a7-4ee8-91dd-171d5665df5e
# ╟─bb7b2a43-4e93-4bfb-a41f-65ff6376babd
# ╟─5a6887a8-e651-4e4b-96c6-d1d7058cd26d
# ╟─4720b6c3-d9bf-4a46-9d38-2cee67e35142
# ╠═1186d2ae-f51e-4237-a14f-e4275b5b304c
# ╠═21796460-d43b-4554-adbc-013f58cf92c3
# ╟─15410cfb-db47-48e0-a7e0-49a061aaaec5
# ╟─ef865aba-8f80-4f38-9560-feecdb9013a7
# ╟─9cfa6c75-cb75-470d-b037-2f01c2f22ae9
# ╠═65523ce4-ebb8-4fb2-9869-61890e0b7e2d
# ╟─9fc30c4d-e5c2-4359-9bec-e93ca5876763
# ╟─665919fd-592f-494a-9d17-adde3434a7e0
# ╟─7c40b5d3-8afd-4648-91ef-ec6235ec376e
# ╟─92791a5f-02f8-461e-9bb1-4e1077f24be9
# ╠═e24e5718-1a1a-47bc-8f25-3268ce9a5826
# ╟─2f9ccd3c-14a1-486e-8b6c-707deda8b4bc
# ╟─13c5b2d6-7f30-4afd-8a08-3dfe4f2db184
# ╠═355bd86a-a278-4ade-96a7-2a18e0c8f7db
# ╟─297f8790-3973-462e-b14a-bff39f399e8e
# ╟─e24068de-07a2-41df-a686-01b3e21ca6b0
# ╟─4bbfd29d-832b-4e17-8f59-25463466fa37
# ╟─ee32ef8b-3b5f-4bd8-bafc-588166899782
# ╟─21f916fe-21b9-43c1-b9d1-a8e541d24e43
# ╟─0f634e30-4f59-48d4-9843-8bbb1e244db5
# ╠═b7dc76de-d6df-414e-bc76-fa9c5e2af727
# ╟─e15fdd6e-ed3e-49bf-a3a7-7530cda0d479
# ╟─5ed2c935-22d6-4e3c-ab2a-a52dd8ee462b
# ╠═ae805230-574e-4818-bdbc-1dda19fa33b0
# ╟─2e865221-74b1-4023-962c-77203686cf03
# ╟─a0af5567-02c6-4c3a-a1f2-649c622d1319
# ╟─f5506679-584a-4f57-93c5-5b71f43b0896
# ╟─a1cfd851-c006-4506-be18-7571927fb5ef
# ╟─4e5aefd4-b965-40d6-ae04-cda4f36945d2
# ╟─f89404e4-9287-4014-bd2d-93f6833a1985
# ╟─8efa2b75-f85a-4be9-85cd-6cccd193bded
# ╟─dbf408ba-bf6e-440c-8855-3f7b7f254d6e
# ╟─74097a18-b4fe-4ec8-b47c-b279fce28598
# ╟─2eb53b12-d20d-411c-9da1-96a0b6ecac39
# ╟─9ad83c33-f843-4d63-b0f3-1a0401068c80
# ╟─f31d5edd-722e-4c55-9ef3-190b4996fd29
# ╠═7c2e9206-8403-4eab-a39e-6625945dd0c1
# ╟─c727e963-6afb-41f5-9de8-38026ab126e0
# ╟─6217aa0d-46e1-4dd6-a1e5-5eb2696757e4
# ╟─5fc13f1f-5818-4ebd-8531-fa5900c47bc9
# ╟─165cfb6b-63f6-4857-849b-559ba9221768
# ╟─13aff5d5-e813-40de-adc1-21fa8d07c812
# ╟─e4350897-af1a-4081-93aa-2eb24b3b5a0a
# ╟─c50531af-4176-40a9-9812-1fae66ffb9f0
# ╟─2143cef3-f43c-422d-b738-1aadbd72660a
# ╟─c6bca558-46de-40c1-b54d-7c61812a62c7
# ╟─4292d332-ac6e-4b40-96e6-dbf086c1118d
# ╠═a52638e9-3115-4951-917a-a7b325ea94c1
# ╟─6b571796-a8c5-40fe-838e-fbbadd2469ad
# ╟─d7e82187-d508-4447-a151-314f85749517
# ╠═954fbbeb-8ca5-4e3d-96de-1b4a61c645fe
# ╠═d619f2e2-31ad-4a56-9225-c99a0b635331
# ╟─7dde9e2b-462f-4f52-be7d-a1e50300efa1
# ╟─4eef11e8-b650-423d-8288-f11f84c7bfc9
# ╠═2a350f62-caca-4139-96b1-06186b043755
# ╟─c3472420-ad6f-4146-89bf-3adcdd57200a
# ╟─bd2e2f6b-0b66-4907-abc3-fac25f50d89f
# ╟─cc4ee8ce-68fe-43b6-92ae-e59b05f87eb9
# ╟─0348cadb-80c0-425a-b2b3-4266fee723ec
# ╠═10116502-8a73-45ca-ab72-8e5b6ce6d115
# ╠═f0c70407-13be-449c-ac0c-fce644caaf2f
# ╟─d19e8ee2-f546-4667-b38a-6af8006158dc
# ╠═109a5754-511c-4df1-b79b-a145e35504f8
# ╠═671f56ef-689b-4300-8ca0-5fc9e32d2146
# ╠═14427835-cace-4a5f-94c9-bac9463cc14c
# ╟─0ef60b26-4769-45d8-b1a3-e2a87bb0d8a2
# ╠═8fde2a0a-9dcd-4b81-b7ae-4249587ebaac
# ╟─1bf61994-f5d6-4ec6-b088-4f60a44953b3
# ╠═1425d9ea-d1ad-4b86-ab9d-ecea3a73a3df
# ╟─b061d7f4-1d38-4f93-abce-8be49363204c
# ╠═e4b5d949-c047-48b6-8f51-bf17bd152f5f
# ╟─254bfbc4-117b-4518-a188-c19a7c123035
# ╟─57160c99-cccc-4aca-99fd-df9356da76a1
# ╟─363968d9-23c9-45e7-b47d-5d7e8ad5648f
# ╟─c65c438b-d97f-4857-9897-20180cc53d25
# ╟─78e2208d-f9f8-44f7-806f-4590db85fe1a
# ╟─a7aa6fe4-3ee2-4dff-a792-4289d4fca83c
# ╟─e25e883e-e0f9-4848-adfd-24f155ecf902
# ╟─21f882ba-fa33-4695-9324-b058ba64d258
# ╟─167a446c-24b7-414c-9bdf-bae7c90df2ab
# ╟─7e94c5d7-ba40-4850-9602-18fe63381ee8
# ╟─8b8a87c1-871d-4fab-bf1a-6b3d8bc9a0fd
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
