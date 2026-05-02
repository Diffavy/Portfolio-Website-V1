### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 554b9b10-6f27-11ee-2678-f10719029989
using Plots ; using PlutoUI ;using Statistics ; plotly() ; gr()

# ╔═╡ 56c20ce8-e208-46db-80ae-ecbe7df7ad9e
md"# *Investigation of electrons tunnelling through potential barriers*"

# ╔═╡ cead5407-18d3-4ca6-b98e-e852308e1937
md"### F132357 - 12/11/2023" 

# ╔═╡ 92fac372-efd4-4dc9-80f0-fded906ec865
md"Throughout this notebook, I will go through how electrons tunnel through different potential barriers and how the wavefunctions, transmission and reflection probabilities describe quantum tunnelling." 

# ╔═╡ 9bf9ca86-e6ca-43b8-bac4-bea979cf728b
md"To describe the motion of the electron and how it behaves in different regions, I use the Time Independent Schrodinger Equation (TISE) however I will look at this for a constant potential **U:**"

# ╔═╡ d50eff2a-55f2-4681-95c6-31d076f1c709
md"""
$$-\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2} + U\psi(x) = E\psi(x) \tag{5}$$"""

# ╔═╡ a841a6ea-ff03-474c-a942-4c3cb74257f5
md"The general solution gives the wavefunction showing how the electron behaves in different regions, with different potentials (**U**) in each region:"

# ╔═╡ 2a2068d1-6419-44ef-89e4-2dc192991b9e
md"""
$$\psi(x) = A exp(ikx) + B exp(-ikx) \tag{5}$$"""

# ╔═╡ e25552ec-811f-48b9-a1f6-9615ae05e554
md"A and B are coefficients that depend on the initial conditions. And the wavevector, k:"

# ╔═╡ 0ba047a4-b6f9-4ca3-9597-c3e2f7418bfb
md"""
$$k =\frac{\sqrt{2m(E - U)}}{\hbar} \tag{5}$$"""

# ╔═╡ 4aeaeaea-10ae-4d29-aba6-3bb47c5ebdf0
md"### Wavefunctions for two regions"

# ╔═╡ ff4f1c98-5b3d-4823-b900-d9ba5dde07ef
md"For this case there are two regions, with two separate potentials that are constant, $U_1$ and $U_2$, for their respective regions $1$ and $2$. Due to the separate potentials, as shown in the relations of the wavevector and potential, there are separate wavevectors and in turn separate wavefunctions. The functions operate in two different domains (Regions) $R_1$ and $R_2$. Take the boundary value between regions $1$ and $2$ to be $x=a$."

# ╔═╡ 896d5022-0e2d-4610-89f2-fda508b6d41c
md"$x<a = R_1$
$\psi_1(x) = A_1 exp(ik_1x) + B_1 exp(-ik_1x)$ 
$$\tag{5}$$
$x>a = R_2$
$\psi_1(x) = A_2 exp(ik_2x) + B_2 exp(-ik_2x)$"

# ╔═╡ 9f4f6842-89a3-4cc9-955a-60b484e1a8c3
md"### The transfer matrix method"

# ╔═╡ ac810460-6acf-4d2c-bd2e-6d97cebefed7
md"The coefficients $A_2$ and $B_2$ are values put in when initialising the wavefunction. The terms with the coefficients $A$ and $B$ define the components of the wavefunction going left to right and right to left respectively. In the case of region $2$ as there are no more subsequent boundaries there will be no reflected waves which means there will be no components going right to left, hence: $B_2 = 0$. To define the coefficients $A_1$ and $B_1$, the transfer matrix method is helpful.

By equating the two wavefunctions, taking their derivatives you can relate the coefficients $A_1$ and $B_1$ to the existing $A_2$ and $B_2$. By equating them, you ensure that the wavefunction 1 & 2 are equal at the boundary of the two regions. Taking the derivatives ensures the transition between the regions is smooth:"

# ╔═╡ c85bad5c-0491-429f-b02e-998c59da6b23
md"""
$$T_{1,2} =\frac{1}{2k_1} \begin{pmatrix}
(k_1 + k_2)exp(-ia(k_1 - k_2)) & (k_1 - k_2)exp(-ia(k_1 + k_2))\\
(k_1 - k_2)exp(ia(k_1 + k_2)) & (k_1 + k_2)exp(ia(k_1 - k_2))
\end{pmatrix}$$

$$\begin{pmatrix}
A_1\\
B_1
\end{pmatrix} = T_{1,2} \begin{pmatrix}
A_2\\
B_2
\end{pmatrix} \tag{5}$$"""


# ╔═╡ 332f42be-7674-4717-b2c3-e522ac0970f3
md"The transfer matrix method is efficient, but it works less beyond one dimensional systems, therefore it cannot accurately be related to real world situations. It works well for situations such as a thin layered film, but it assumes no thickness. Moreover, it's commonly used for the (TISE) and less so for the Time Dependent Schrodinger equation (TDSE)."

# ╔═╡ fe63d788-3a0a-4748-b577-00e67e031576
md"I have defined the different packages I use to create the plots and the different functions I use from the Statistics package. The plotly and gr backends are for analysing and outputting the graphs."

# ╔═╡ 8af54e36-f964-426b-b790-2ad1a0f1ac94
md"I then defined constants that I used throughout, $i$, $\hbar$ and the mass of an electron as variables used throughout."

# ╔═╡ 82789e60-7dc4-4452-b328-57c9e7cdab22
i = 1im ; hbar = (6.63*10^(-34))/(2*pi) ; m = 9.11*10^(-31)

# ╔═╡ c264e0ac-f50c-4fce-acf2-421bf9d8f01e
md"The following function is to ensure that arrays are fully defined (i.e. not including an infinity or NaN elements) as a form of testing the function when building it to find errors. It will not be used throughout the examples as they have already been tested."

# ╔═╡ 26e91094-9a06-4a9b-a58e-d9510c07af8c
function isArrayDefined(arr)
    return !any(isnan.(arr)) && !any(isinf.(arr))
end

# ╔═╡ 9195cd37-9de5-4797-96ba-0e2be3002dad
md"This function is for an arbitrary transfer matrix between two arbitrary regions ($R_i$ and $R_{i+1}$), where the wavevectors $k_i$ and $k_{i+1}$ are the wavevectors for the respective regions. And the $x_i$ is the boundary values between the regions:

$T_{i,i+1} =\frac{1}{2k_i} \begin{pmatrix}
(k_i + k_{i+1})exp(-ia(k_i - k_{i+1})) & (k_i - k_{i+1})exp(-ia(k_i + k_{i+1}))\\
(k_i - k_{i+1})exp(ia(k_i + k_{i+1})) & (k_i + k_{i+1})exp(ia(k_i - k_{i+1}))
\end{pmatrix}$
$$\tag{5}$$"

# ╔═╡ 68cc8a11-f3d6-4ac7-baf0-b55cec5e3ee7
function transferMatrix(kiplus1,ki,xi)
	(1/(2*ki))*[(ki+kiplus1)*exp(-i*xi*(ki-kiplus1)) (ki - kiplus1)*exp(-i*xi*(ki + kiplus1)) ; (ki - kiplus1)*exp(i*xi*(ki + kiplus1)) (ki + kiplus1)*exp(i*xi*(ki - kiplus1))]
end

# ╔═╡ 0941fa2c-3249-4bc2-9885-899effc0e3fe
md"This function is to convert between (electron volts) eV and  (joules) J units:"

# ╔═╡ c207319c-dee6-4a87-b968-7ac8c86f047d
function evtoJ(ev::Float64)
	return ev*1.6*10^(-19)
end

# ╔═╡ 49ae3345-2d58-4424-877f-6f5f1de44546
md"Finally, this function is to define the wavevector for an arbitrary region. With these defined the wavefunctions for the different scenarios can be constructed."

# ╔═╡ dc480e3f-828a-4398-a944-68064c708443
function k(m::Float64,E::Float64,U::Float64)
	sqrt(Complex(2*m*(E - U)))/(hbar)
end

# ╔═╡ 81284779-be86-46d4-9e68-817225790be6
md"### 1. Potential step (Two regions)"

# ╔═╡ 309e640a-0fef-4209-893f-3c6068d0bfd2
md"The function starts by defining the wavevectors $k_1$ and $k_2$ using the wavevector function previously defined with the initialised values of $U_1$ and $U_2$. I then define '*deltax*' which is the width of each region. The $dx$ value is an arbitrary value that controls how many equidistance points there are between the boundary values that enclose the two regions defined by $R_1$ and $R_2$."

# ╔═╡ bb9d49af-1772-4c85-ac90-29e19f533c86
md"Using the defined wavevectors, and the previously defined transfer matrix function, the transfer matrix acts on the vector of the $A_2$ and $B_2$ coefficients and produces the coefficients of $A_1$ and $B_1$. Now using the coefficients, wavevectors and the domains of the Regions $1$ and $2$ to obtain an array of $\psi_1$ and $\psi_2$ to plot against their respective domains $R_1$ and $R_2$."

# ╔═╡ 839ad551-f112-4899-a0dd-ba63fdfe5788
md"The transmission and reflection coefficients determine the probability that the electron goes through the potential boundary between the regions.

Transmission coefficient :

$Tr = |\frac{1}{t_{11}}|^2 \frac{k_N}{k_1}$

Reflection coefficient :

$$R  = |\frac{t_{21}}{t_{11}}|^2 \tag{5}$$
"

# ╔═╡ da9904d7-3ab6-4cb5-b66a-f2bc94c0585b
md"Where $t_{ij}$ are the elements of the transfer matrix and $k_N$ is the wavevector of the final region $N$"

# ╔═╡ 0a135be0-38ce-45d7-80d5-e10ec09f81ba
function waveFunction(m::Float64,U1::Float64,U2::Float64,E::Float64,A2::Float64,B2::Float64,a::Float64,RorI::String)
	k1 = k(m,E,U1) #wavevector in region 1
	k2 = k(m,E,U2) #wavevector in region 2

	deltax = 5*10^(-9)
	dx = 0.01*10^(-9)
	#Region 1 x<a
	R1 = (a*10^(-9)-deltax):dx:a*10^(-9)
	#Region 2 x>a
	R2 = a*10^(-9):dx:(a*10^(-9) + deltax)

	T12 = transferMatrix(k2, k1, a)

	#Transmission and reflection coefficients
	Tr = (1/T12[1][1])^(2) * (k2/k1)
	R = (T12[2][1]/T12[1][1])^(2)
	
	println("Transmission =") ; println(real(Tr))
	println("Reflection =") ; println(real(R))
	
	A2B2 = [A2 ; B2]
	A1B1 = T12 * A2B2
	A1 = A1B1[1]
	B1 = A1B1[2]
	
	#Region 1 when x<a
	psi1 = A1.*exp.(i*k1.*R1) .+ B1.*exp.(-i*k1.*R1)
	#Region 2 when x>a
	psi2 = A2.*exp.(i*k2.*R2) .+ B2.*exp.(-i*k2.*R2)

	if RorI == "Real"
		#real components
		plot(R1 ,real(psi1) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi1" )
		plot!(R2 ,real(psi2) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi2")
		plot!(title="Wave function (real components) in regions 1 & 2")

	elseif RorI == "Imag"
		#imaginary components
		plot(R1 ,imag(psi1) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi1" )
		plot!(R2 ,imag(psi2) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi2")
		plot!(title="Wave function (imaginary components) in regions 1 & 2")

	elseif RorI != "Real" || "Imag"
		md"##### Input either the string 'Real' or 'Imag' to plot either the real or imaginary components of the wavefunction"
	end
end

# ╔═╡ d9bb8263-5b0d-485a-9e6a-bdb02f3c8443
md"Syntax to call the function: waveFunction($m,U_1,U_2,E,A_2,B_2,a,RorI$)"

# ╔═╡ e2e9f666-f840-47d8-9a73-04157395d822
md"The wavefunction defined above takes in the following arguments the mass of the given particle, the constant potentials in the regions $1$ & $2$, the energy of the particle, the coefficients of region 2 $A_2$ and $B_2$, the boundary value between the regions 'a' and the string 'RorI' which links to an if statement that plots either the real or imaginary components of the wavefunction."

# ╔═╡ 3a4dba24-705b-4329-a8b4-6a31e2d5eb87
md"To initialise the function I plotted $(1.1)$ the behaviour of an electron incident on a potential barrier where the energy of the electron $(E =0.5eV)$ is less than the potential $(U_2 =2.0eV)$ of the barrier. As expected the electron is incident and propagates in Region $1$ and it begins to decay once it hits the potential barrier and enters Region $2$."

# ╔═╡ 136e351f-f373-43a9-a799-3bb4e820d49e
waveFunction(m,0.0,evtoJ(2.0),evtoJ(0.5),1.0,0.0,0.0,"Real")

# ╔═╡ 5728e695-b021-4ba8-b2b0-d1c1e4592de8
md"""
$$\tag{Figure 1.1}$$"""

# ╔═╡ d9977075-05ea-4625-9002-7f2a767544c7
md"When increasing the energy of the incident electron (e.g. to $E = 1.8eV$ [as shown below]) yet not greater than the potential energy of the potential step, the wavelength of the incident electron decreases showing more oscillations in Region $1$. However, interestingly the point of decay stays constant despite the changing of the incident electron's energy. This implies that the two following expressions do not depend on the energy of the incident particle.

$\psi_1(a) = \psi_2(a)$
$$\frac{d\psi_1}{dx}(a) = \frac{d\psi_2}{dx}(a)$$"

# ╔═╡ 8ac53ff1-1b91-44a8-ace2-81134f67743b
waveFunction(m,0.0,evtoJ(2.0),evtoJ(1.8),1.0,0.0,0.0,"Real")

# ╔═╡ 6c102db4-74e7-42fc-b310-03ec11d53642
md"""
$$\tag{Figure 1.2}$$"""

# ╔═╡ 11b9d64a-f4f1-4571-ad14-2a8959a81ff0
md"When the energy of the incident electron is equal to the energy of the potential barrier $E = U_2$ $($shown in $1.3)$, the transmission coefficient and the reflection coefficient are $0$  and $1$ respectively. This means there's complete reflection and no transmission. This phenomenon is called resonance/particle trapping. Which as the name suggests traps the particle within a region (in this case Region $1$)."

# ╔═╡ 64618b84-5fa2-45f1-a697-57cc3fa8378c
waveFunction(m,0.0,evtoJ(2.0),evtoJ(2.0),1.0,0.0,0.0,"Real")

# ╔═╡ 0a509173-c9ff-431c-a5ff-2f50393465da
md"""
$$\tag{Figure 1.3}$$"""

# ╔═╡ 960ce184-befd-41f0-8f13-84ea7c3cd007
md"Finally, as expected when the energy of the electron is greater than that of the potential step it oscillates beyond the barrier (with less potential energy) and does not decay as shown (1.4):"

# ╔═╡ 9d59944a-5a39-4157-beb0-d2a4d6c5dbfd
waveFunction(m,0.0,evtoJ(2.0),evtoJ(2.5),1.0,0.0,0.0,"Real")

# ╔═╡ 4aae4592-19c1-41dd-aadf-2c96e68df96f
md"""
$$\tag{Figure 1.4}$$"""

# ╔═╡ 583a7e06-2a9e-481e-b21d-923f4012bf9a
md"### 2. Single Finite Barrier (Three regions)"

# ╔═╡ c8c81211-f8c5-412b-a297-14ae8f12359e
md"To model the transmission probability as a function of a variable energy (energy of the electron), I modelled the potential step in three regions. The code is for the most part the same as for the potential step, with some minor changes to upscale for three regions. The changes being two transfer matrices take the coefficients from Region $3$ to $2$ and then Region $2$ to $1$. There was a single boundary value for the potential step (Section $1$) $a$, however there are now two $'x_1$ and $x_2'$ which are the boundaries between regions $1$ & $2$ and $2$ & $3$ respectively."

# ╔═╡ 22edd46c-a68b-4e6b-a419-7940ea18efc0
md"To model this as a potential step as before $U_2$ is the same constant potential as before and $U_1$ is $0$ as done before, $U_3$ is also $0$. The main difference is the introduction of $E1$ which is the range of values upon which $T(E)$ is defined, $a$ which is the number of elements that $E1$ and $T(E)$ have. The 'ones' initializes an array of ones with $a$ elements to overwrite to with calculated values."

# ╔═╡ c5918e70-99f8-488f-996e-7153dc6c12e2
md"As the energy of the incident electron changes there are three cases that change how it relates to the transmission probability. The $\delta$ is the width of the region. From $[1]$ Page $38$:

$k = \sqrt{2mE/\hbar^2}$

when: $E<V$ 

$$\kappa = \sqrt{2m(V-E)/\hbar^2} \tag{1}$$

$T(E) = \frac{1}{1 + \frac{1}{4}(\frac{\kappa}{k} + \frac{k}{\kappa})^2\sinh^2{\kappa\delta}}$

when: $$E=V$$

$$\frac{1}{1+(\frac{k\delta}{2})^2} \tag{1}$$

when: $E>V$

$$k^\prime = \sqrt{2m(E - V)/\hbar^2} \tag{1}$$

$T(E) = \frac{1}{1 + \frac{1}{4}(\frac{k^\prime}{k} - \frac{k}{k^\prime})^2\sin^2{k^\prime\delta}}$"

# ╔═╡ e25e8620-8713-401e-83de-083fc269763c
#=
	#The energy dependence of the transmission (Tr) probability [T(E)]
	dE = evtoJ(0.01)
	E1 = dE:dE:evtoJ(20.0)
	a = Int(E1[end]/dE)
	TE = ones(a)

	for i in 1:a
		kappa = i*k(m,E1[i],U2)
		kprime = k(m,E1[i],U2)
		ka = k(m,E1[i],0.0)
		
		if E1[i] < U2
			TE[i] = 1/(1 + (0.25)*(kappa/ka + ka/kappa)^(2)*(sinh(kappa*deltax))^2)
		elseif E1[i] == U2
			TE[i] = 1/(1 + (ka*deltax/2)^2)
		elseif E1[i] > U2
			TE[i] = 1/(1 + (0.25)*(kprime/ka - ka/kprime)^(2)*(sin(kprime*deltax))^2)
		end
	end
=#

# ╔═╡ f936738a-42d9-4e9b-bcbb-86a252998c79
function waveFunction2(m::Float64,U1::Float64,U2::Float64,U3::Float64,E::Float64,A3::Float64,B3::Float64,x1::Float64,RorIorTr::String)
	deltax = 5*10^(-9)
	#Wavevectors in each region
	k1 = k(m,E,U1) ; k2 = k(m,E,U2) ; k3 = k(m,E,U3)
	#Boundary values
	x1 = x1*10^(-9) ; x2 = x1 + deltax
	
	dx = (0.01*10^(-9))
	
	#Region 1
	R1 = (x1-deltax):dx:x1
	#Region 2 
	R2 = x1:dx:x2
	#Region 3
	R3 = x2:dx:(x2+deltax)
	
	#Defining the coefficients of the wavefunctions
	T12 = transferMatrix(k2, k1, x1)
	T23 = transferMatrix(k3, k2, x2)
	
	A3B3 = [A3 ; B3]
	A2B2 = T23 * A3B3
	A1B1 = T12 * A2B2
	A1 = A1B1[1] ; B1 = A1B1[2] ; A2 = A2B2[1] ; B2 = A2B2[2]
	
	#when Region 1 (0<R1<x1)
	psi1 = A1.*exp.(i*k1.*R1) .+ B1.*exp.(-i*k1.*R1)
	#when Region 2 (x1<R2<x2)
	psi2 = A2.*exp.(i*k2.*R2) .+ B2.*exp.(-i*k2.*R2)
	#when Region 3 (x2<R3<(x2+deltax))
	psi3 = A3.*exp.(i*k3.*R3) .+ B3.*exp.(-i*k3.*R3)

	#The energy dependence of the transmission (Tr) probability [T(E)]
	dE = evtoJ(0.01)
	E1 = dE:dE:evtoJ(20.0)
	a = Int(E1[end]/dE)
	TE = ones(a)

	for i in 1:a
		kappa = i*k(m,E1[i],U2)
		kprime = k(m,E1[i],U2)
		ka = k(m,E1[i],0.0)
		
		if E1[i] < U2
			TE[i] = 1/(1 + (0.25)*(kappa/ka + ka/kappa)^(2)*(sinh(kappa*deltax))^2)
		elseif E1[i] == U2
			TE[i] = 1/(1 + (ka*deltax/2)^2)
		elseif E1[i] > U2
			TE[i] = 1/(1 + (0.25)*(kprime/ka - ka/kprime)^(2)*(sin(kprime*deltax))^2)
		end
	end
	
	if RorIorTr == "Real"
		#real components
		plot(R1 ,real(psi1) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi1" )
		plot!(R2 ,real(psi2) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi2")
		plot!(R3 ,real(psi3) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi3")
		plot!(title="Psi")

	elseif RorIorTr == "Imag"
		#imaginary components
		plot(R1 ,imag(psi1) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi1" )
		plot!(R2 ,imag(psi2) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi2")
		plot!(R3 ,imag(psi3) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi3")
		plot!(title="Psi")


	elseif RorIorTr == "Tr"
		plot(E1, TE,xlabel= "E (J)",ylabel = "Tr(E[J])", label= "T(E)")
		plot!(title="Transmission Probability")
	elseif RorIorTr != "Real" || "Imag" || "Tr"
		md"##### Input either the string 'Real' or 'Imag' to plot either the real or imaginary components of the wavefunction. And to plot the transmission as a function of energy, 'Tr'."
	end
	
end


# ╔═╡ 71c797bf-abd2-4791-a367-a52b7c083b52
md"Syntax to call the function: $waveFunction2(m,U_1,U_2,U_3,E,A_3,B_3,a,RorIorTr)$"

# ╔═╡ 1b6d2638-4a48-4eb3-a95e-5ecdab18bc73
md"With these initial conditions, the finite barrier outputs the same wavefunction as for the potential step. See, $(2.1)$ :"

# ╔═╡ 8757956d-6797-40ba-ae3a-d087fba12999
waveFunction2(9.11*10^(-31) ,0.0, evtoJ(2.0), 0.0, evtoJ(0.5), 1.0, 0.0, 0.0,"Real")

# ╔═╡ 71f4c76b-fc1e-4828-86e7-db49fbd311ee
md"""
$$\tag{Figure 2.1}$$"""

# ╔═╡ ba34511b-367f-49d6-9cee-3da45070ad75
waveFunction2(m ,0.0, evtoJ(2.0), 0.0, evtoJ(0.5), 1.0, 0.0, 0.0,"Tr")

# ╔═╡ 6cda3ee2-eebf-4330-a713-d1ad8f6e6b7f
md"""
$$\tag{Figure 2.2}$$"""

# ╔═╡ 4007531f-5c53-47bc-8400-5d2c8fd7786c
md"The graph above $(2.2)$ shows how the transmission probability as a function of energy and how it changes for a finite barrier. The shape follows the expected outcome, where the probability increases as the energy increases as it oscillates. If you were to fit a curve through the troughs of the oscillation you would obtain the lower bound of transmission probability and if you were to fit a curve to the peaks you would have the upper bounds of the transmission probability. The also lines up with the literature when compared with the same initial conditions [1] -Page 39 (Figure 10)."

# ╔═╡ 6e0fc9bf-8abb-444b-9653-af7a97e61f2b
md"### 3. Potential step : (Five regions) [Double Barrier]"

# ╔═╡ 3dca1e0b-fbe5-42e6-87e2-81d283491829
md"The following function is for a double barrier system of an electron incident travelling from left to right. Defining the distance between the two barriers $w$, the width of both barriers $b$, the start of the system's domain $x_0$ and the potential of each barrier $U_1$ and $U_2$ respectively."

# ╔═╡ 2726ed0f-d546-43a0-b9f2-3c1728054126
md"The code follows a similar structure upscaled to the necessary five regions. Where $R_1$ is the first region of an arbitrary width $a$ and $U=0.0eV$, followed by $R_2$ which contains the first barrier with width $b$, then $R_3$ which is the gap between the two barriers of width $w$ with the potential as $U$. Finally, $R_4$ is the same as $R_2$ except with a potential of $U_2$ and $R_5$ is the same as $R_1$ just on the other side (they act as reference regions to see how the particle enters the double barrier system $[$in $R_1]$ and leaves the double barrier $[$in $R_5]$)."

# ╔═╡ 39e540d3-d039-4f09-8469-252055162ae9
md"As the potentials of regions $1$, $3$ & $5$ are all zero $(U)$, they all have the same wavevector defined by $k_i$ [shown below]. To create a continuous wavefunction, showing the system in its entirety, the boundary values have to be defined relative to the start of systems domain $x_0$ which is $0.0$ for ease of understanding in my plots. Finally, $dx$ is used to define the equidistant data points in each region."

# ╔═╡ 059713dc-c14b-4261-822f-84650f595775
#=
	#Potential and wavevector in Regions 1,3 & 5 (No potential)
	U = 0.0 ; ki = k(m,E,U)
	#Wavevectors in Regions 2 & 4 (Barrier regions)
	k1 = k(m,E,U1) ; k2 = k(m,E,U2)

	#Boundary values
	x1 = a ; x2= a+b ; x3= a+b+w ; x4= a+2*b+w ; x5= 2*a+2*b+w

	#Defining each region
	dx = 0.01 * 10^(-9)
	R1 = x0:dx:x1 ; R2 = x1:dx:x2 ; R3 = x2:dx:x3 ; R4 = x3:dx:x4 ; R5 = x4:dx:x5
=#

# ╔═╡ 6d9824c2-a1d3-4168-933b-d8a8b7a5a7f9
md"As explained in Section $2$ the final region has coefficients of $1.0$ and $0.0$, which in this case are the coefficients $A_5$ and $B_5$ respectively. To obtain the coefficients in each region I define the transfer matrices once again (this time $4$ of them as there are $5$ regions). Using those coefficients and the defined regions, the wavefunctions are formed. "

# ╔═╡ 87154415-4a49-474e-bb21-a5424d4de711
md"The code section for the energy dependence of the transmission probability '$Tr(E)$'  has not changed, and the if statement '$RorIorTr$' is for the input string which lets you plot either the 'Real', 'Imag' or 'Tr' as done in the previous two systems."

# ╔═╡ 0e676326-bc3e-4f73-9d62-99d3e02169c7
function waveFunction3(m::Float64,U1::Float64,U2::Float64,E::Float64,a::Float64,b::Float64,w::Float64,A5::Float64,B5::Float64,x0::Float64,RorIorTr::String)
	#Potential and wavevector in Regions 1,3 & 5 (No potential)
	U = 0.0 ; ki = k(m,E,U)
	#Wavevectors in Regions 2 & 4 (Barrier regions)
	k1 = k(m,E,U1) ; k2 = k(m,E,U2)

	#Boundary values
	x1 = a ; x2= a+b ; x3= a+b+w ; x4= a+2*b+w ; x5= 2*a+2*b+w

	#Defining each region
	dx = 0.01 * 10^(-9)
	R1 = x0:dx:x1 ; R2 = x1:dx:x2 ; R3 = x2:dx:x3 ; R4 = x3:dx:x4 ; R5 = x4:dx:x5
	
	#Defining the coefficients of the wavefunctions
	T12 = transferMatrix(k1,ki,x1)
	T23 = transferMatrix(ki,k1,x2)
	T34 = transferMatrix(k2,ki,x3)
	T45 = transferMatrix(ki,k2,x4)
	
	A5B5 = [A5 ; B5]
	A4B4 = T45*A5B5 ; A3B3 = T34*A4B4 ; A2B2 = T23*A3B3 ; A1B1 = T12*A2B2
	A1=A1B1[1] ; B1=A1B1[2] ; A2=A2B2[1] ; B2=A2B2[2] ; A3=A3B3[1] ; B3=A3B3[2]
	A4=A4B4[1] ; B4=A4B4[2] ; A5=A5B5[1] ; B5=A5B5[2]

	#Defining the wavefunctions in each region
	psi1 = A1.*exp.(i*ki.*R1) .+ B1.*exp.(-i*ki.*R1)
	psi2 = A2.*exp.(i*k1.*R2) .+ B2.*exp.(-i*k1.*R2)
	psi3 = A3.*exp.(i*ki.*R3) .+ B3.*exp.(-i*ki.*R3)
	psi4 = A4.*exp.(i*k2.*R4) .+ B4.*exp.(-i*k2.*R4)
	psi5 = A5.*exp.(i*ki.*R5) .+ B5.*exp.(-i*ki.*R5)
	
	#The energy dependence of the transmission (Tr) probability [Tr(E)]
	dE = evtoJ(0.01)
	E1 = dE:dE:evtoJ(20.0)
	a = Int(E1[end]/dE)
	TE = ones(a)
	
	for i in 1:a
		kappa = k(m,E1[i],U2)
		kprime = k(m,E1[i],U2)
		ka = k(m,E1[i],0.0)
		
		if E1[i] < U2
			TE[i] = 1/(1 + (0.25)*(kappa/ka + ka/kappa)^(2)*(sinh(kappa*b))^2)
		elseif E1[i] == U2
			TE[i] = 1/(1 + (ka*b/2)^2)
		elseif E1[i] > U2
			TE[i] = 1/(1 + (0.25)*(kprime/ka - ka/kprime)^(2)*(sin(kprime*b))^2)
		end
	end
	
	if RorIorTr == "Real"
		#real components
		plot(R1 ,real(psi1) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi1" )
		plot!(R2 ,real(psi2) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi2")
		plot!(R3 ,real(psi3) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi3")
		plot!(R4 ,real(psi4) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi4")
		plot!(R5 ,real(psi5) ,xlabel= "x (m)",ylabel = "psi m^(-1/2)", label= "psi5")
		plot!(title="Psi")

	elseif RorIorTr == "Imag"
		#imaginary components
		plot(R1 ,imag(psi1) ,xlabel= "x",ylabel = "psi", label= "psi1" )
		plot!(R2 ,imag(psi2) ,xlabel= "x",ylabel = "psi", label= "psi2")
		plot!(R3 ,imag(psi3) ,xlabel= "x",ylabel = "psi", label= "psi3")
		plot!(R4 ,imag(psi4) ,xlabel= "x",ylabel = "psi", label= "psi4")
		plot!(R5 ,imag(psi5) ,xlabel= "x",ylabel = "psi", label= "psi5")
		plot!(title="Psi")


	elseif RorIorTr == "Tr"
		plot(E1, TE,xlabel= "E(J)",ylabel = "Tr(E[J])", label= "T(E)")
		plot!(title="Transmission Probability")
		
	elseif RorIorTr != "Real" || "Imag" || "Tr"
		md"##### Input either the string 'Real' or 'Imag' to plot either the real or imaginary components of the wavefunction. And to plot the transmission as a function of energy, 'Tr'."
	end
	
end

# ╔═╡ 269cf590-06cc-402d-8c07-77c85055fbdf
md"Syntax to call the function: $waveFunction3(m,U_1,U_2,E,a,b,w,A_5,B_5,x_0,RorIorTr)$"

# ╔═╡ 5dcd0bf5-b13f-4ea3-92c7-9c66f2fac085
md"In $(3.1)$ I initialised the wavefunction such that the potential of the first barrier $'U_1 = 2.0eV'$ is marginally less than the energy of the particle $'E = 2.5eV'$, and the potential of the second barrier is marginally greater $'U_2 = 3.0eV'$. If you first look at the energy of the particle before reaching the first barrier in blue $'\psi_1'$ and compare it after going through the first barrier in Region $3$ (green $'\psi_3'$) its energy has decreased. This is the expected outcome (as shown in $(1.4)$) when the energy of the particle is marginally greater than that of the barrier."

# ╔═╡ d393922b-a0e9-4524-b71f-0310369a0e4d
md"Now as the particle continues towards the second barrier, it begins to decay into the barrier $\psi_4$ , as was shown in $(1.1)$. However due to the smaller width of the potential barrier $(b = 1 * 10^{-9})$, as opposed to the one shown in $(1.1)$, the particle does not fully decay within $\psi_4$ leaving the second barrier with some energy. In $\psi_5$, you can see the electron propagates with some small energy, instead of levelling off."

# ╔═╡ 8ac27cfd-9c9b-4f81-b6ba-adcae7e2e13e
waveFunction3(m,evtoJ(2.0),evtoJ(3.0),evtoJ(2.5),3*10^(-9),1*10^(-9),5*10^(-9),1.0,0.0,0.0,"Real") 

# ╔═╡ a03ef348-17b8-4058-9fea-5e36c3f9a1ec
md""" 
$$\tag{Figure 3.1}$$"""

# ╔═╡ 9335ffb6-098e-44db-8ac4-554b3a3495dc
md"In $(3.2)$ I tripled the width of the potentials '$b$' to reaffirm the conclusion I came to. Here the particle decays further, as it has more to tunnel through, to the point of negligible energy. Clearly, here there are no oscillations in $\psi_5$ when compared to the same region in $(3.1)$."

# ╔═╡ 4794b945-9ac0-465c-ad65-713331acf94e
waveFunction3(m,evtoJ(2.0),evtoJ(3.0),evtoJ(2.5),3*10^(-9),3*10^(-9),5*10^(-9),1.0,0.0,0.0,"Real") 

# ╔═╡ 9be1351c-8fd5-412f-b852-6226c14b5f7d
md""" 
$$\tag{Figure 3.2}$$"""

# ╔═╡ 49dd843d-9791-427b-842b-2fb1d3642486
md"In the example below $(3.3)$ I show how the system evolves when the potentials of both barriers are equal $U_1 = U_2 = 2.5eV$, but the energy of the electron is marginally greater $E = 2.6eV$ than the potential of both barriers. The energy of the electron is dampened as it goes through each barrier."

# ╔═╡ 5674456b-a714-4b33-aa17-ec5656328e30
waveFunction3(m,evtoJ(2.5),evtoJ(2.5),evtoJ(2.6),3*10^(-9),1*10^(-9),5*10^(-9),1.0,0.0,0.0,"Real") 

# ╔═╡ 813dcd90-6bad-4811-900d-8f8f3ea67fd0
md""" 
$$\tag{Figure 3.3}$$"""

# ╔═╡ af6c21bc-fd3f-46e7-830d-58fae6dab7d3
md"But when the energy of the electron is marginally smaller than the potentials $E = 2.4 eV$ the effect of the dampening is far more pronounced. This demonstrates the phenomena of quantum tunnelling for an electron, and how exaggerated the rate of decay is (especially in comparison to the classical case where the electron has greater energy) relative to the example $(3.3)$."

# ╔═╡ 39cbf241-c4af-4411-aaba-532b35eeac3a
waveFunction3(m,evtoJ(2.5),evtoJ(2.5),evtoJ(2.4),3*10^(-9),1*10^(-9),5*10^(-9),1.0,0.0,0.0,"Real") 

# ╔═╡ a9006f30-a9c7-415e-bfa5-e0212fbebd11
md""" 
$$\tag{Figure 3.4}$$"""

# ╔═╡ 7952d09c-586e-48db-a577-65375cdffc8a
md"The graph below follows the expected pattern with the decrease in the displacement of the troughs. The equations used in section $2$ are analytical in nature. It is good for simple systems however as the systems begins to be more complex it begins to falter, becoming less accurate though it is still a good approximation for the double barrier."

# ╔═╡ 80426867-5e2f-41da-851d-1aad5c9126a9
waveFunction3(m,evtoJ(2.5),evtoJ(2.5),evtoJ(2.6),3*10^(-9),1*10^(-9),5*10^(-9),1.0,0.0,0.0,"Tr") 

# ╔═╡ cfad5c07-c4a2-45b9-a4b6-89beaadc8cfd
md""" 
$$\tag{Figure 3.5}$$"""

# ╔═╡ 877bb084-6e7c-46f4-b2ef-c8ba682bcbd5
md"### 4. Variable potential (Arbitrary number of regions)"

# ╔═╡ 8d6ee29a-4a8f-461c-8722-527a5508efe1
md"The function below outputs a variable potential as opposed to the constant potential used in the previous two examples. The variable potential follows a Gaussian shape for the barrier. Using the equation, from $[1]$ - Figure $(116)$ :

$$V(x) = Ve^{-x^2} \tag{1}$$"

# ╔═╡ 3a3cc777-957c-480a-95e2-ff3003242847
function U(x::Float64,U::Float64)
	U * exp(-(x^2))
end

# ╔═╡ 269a6d76-66b5-4230-87fc-19c92785c07b
md"This function truncates a float variable to an arbitrary number of significant figures. I have struggled to find an alternative that doesn't decrease the accuracy. Potentially doing this in a coding language better suited to numerical solutions with maths in mind would let me sidestep this function. Similarly, a suitable package would be useful."

# ╔═╡ 3ce84e57-8be5-4134-a171-2539ad705c1b
function truncateToSf(x::Float64, sf::Int)
    if x != 0.0
        # Calculate the order of magnitude
        order_magnitude = floor(log10(abs(x)))
        # Calculate the scaling factor based on significant figures
        scale_factor = 10^(sf - 1 - order_magnitude)
        # Truncate the scaled number
        truncated_number = trunc(x * scale_factor) / scale_factor
        return truncated_number
    else
        return x
    end
end

# ╔═╡ b2e27d7c-437a-46d1-accf-d6ed1fa6df5f
md"Wavefunction $3$ here takes the coefficients of the final region $N$, $(A_N$ and $B_N)$, which are used to compute the coefficients of the different regions. The $x_1$ and $x_{n-1}$ are the first and last boundary values. $N$, the number of regions. The lowercase $v$ is the value for an arbitrarily chosen region which will output the wavefunction at the boundary of regions $v-1$ and $v$. The uppercase $V$ is to select the coefficient $V$ for the variable potential defined in the equation for the Gaussian barrier above."

# ╔═╡ e9154db9-cf7b-49eb-8b7f-4fd1dc0de6b1
md"The if statement below ensures the selected Region $v$ is not beyond the total number of Regions $N$ selected. And as mentioned previously when selecting an arbitrary region $v$ the wavefunction is plotted at the boundary of $v-1$ and $v$. There it must be 2 or greater as v-1 has to be greater than or equal to 1. Similarly, N has to be at least $2$ for it to be a non-trivial case."

# ╔═╡ c72d20b4-5617-4b7e-b464-9cefab50cd86
#=
	if v > N 
		return md"#### The specified region 'v' to create the wavefunction is more than the total number of regions"
	elseif v == 1
		return md"#### The specified region 'v' to create the wavefunction has to be greater than 1 to compare it to a region one less"
	end
=#

# ╔═╡ 2b4be712-5bdb-45d2-b2b3-7d44c09a9ab6
md"Below I defined the region boundaries where deltax is the constant width of each region. Due to the nature of Julia, there is some inaccurate rounding which is aided using the earlier defined truncating function. The arbitrary boundary values between regions are defined by $BR_i$. Each region has a set number of $x$ data points defined by $y$. $R_i$ is an $N*y$ matrix where each array gives the domain for each region. Shown by the for loop."

# ╔═╡ f677a1e8-1b0d-404e-9ee2-2aac2b84047a
#=
	deltax = (xnsub1 - x1)/(N-2)
	deltax = truncateToSf(deltax, 3)
	
	#Region boundaries
	BRi = (x1 - deltax):deltax:(xnsub1 + deltax)
	
	#Creation of regions
	dx = 0.01*10^(-9)
	y = Int(round(((BRi[2] - BRi[1])/dx)) + 1)
	Ri = ones(Float64, N , y)
	for i in 1:(N)
		Ri[i, :] = BRi[i]:dx:(BRi[i+1] + 0.01*dx)	
	end
=#

# ╔═╡ 485a8dab-575d-4ba0-8cc6-3b3d3d37ee44
md"Following that various wavevectors, coefficients and transfer matrices that are initialised for later use."

# ╔═╡ 77af237b-73c2-46ae-85c5-e5bed1dac5ed
#=
	#Defining the wavevectors at regions '1' 'v' 'v+1' and 'n' by taking the mean of the wavevectors at each point in each region
		kn = k(m,E,U(mean(Ri[N, :]),V)) 
		k1 = k(m,E,U(mean(Ri[1, :]),V))
		kv = k(m,E,U(mean(Ri[v, :]),V))
		kvminus1 = k(m,E,U(mean(Ri[v-1]),V))
		
		AV = 1.0
		BV = 1.0
		AVminus1 = 1.0
		BVminus1 = 1.0
		ANBN = [An ; Bn]
		T1N = [1.0 0.0; 0.0 1.0]
		TVN = [1.0 0.0; 0.0 1.0]
		TVminus1N = [1.0 0.0; 0.0 1.0]
		TE = ones(N-1)
=#

# ╔═╡ 4f3c0619-03de-46af-9791-6582ce745f79
md"In the following piece of code the iterative wavevector $k_i$ for each region is formed using the middle of the domain for each region (by taking the mean) to calculate the average potential for each region. Using $T_{i,i+1}$ and the following if statements I define the transfer matrices for $T_{V,N}$ and $T_{V-1,N}$ to initialise the coefficients $A_V$ and $B_V$ as well as $A_{V-1}$ and $B_{V-1}$ respectively. Following, which I define the wavevectors $\psi_v$ and $\psi_{v-1}$ with all the defined variables." 

# ╔═╡ 8a28ed21-93d6-43a7-9829-9981614aadb1
#=
	#Forms the T1N (transfer matrix) from region 1 to the nth region 
		for i in 1:(N-1)
			ki = k(m,E,U(mean(Ri[i, :]),V)) 
			kiplus1 = k(m,E,U(mean(Ri[i+1, :]),V)) 
			Tiplus1 = transferMatrix(ki,kiplus1,BRi[i+1])
			T1N = T1N * Tiplus1
			
			#Creates a transfer matrix for regions n to v, to help create a wavefunction at an arbitrary region v and the same for region v+1
			if i >= (v)
				TVN = TVN * Tiplus1
				#println(TVN)
			end
			
			if i >= (v-1)
				TVminus1N = TVminus1N * Tiplus1
				#println(TVminus1N)
			end
			
			if i == (N-1)
				AVBV = TVN * ANBN
				AV = AVBV[1]
				BV = AVBV[2]
				
				AVminus1BVminus1 = TVminus1N * ANBN
				AVminus1 = AVminus1BVminus1[1]
				BVminus1 = AVminus1BVminus1[2]
			end
		end
=#

# ╔═╡ 98596d44-87cf-451b-968d-450d9c0c9511
md"Here I define the transmission and reflection probabilities using the equations I previously defined. I used secondary equations for both $(Tr1$ and $R1)$ which rely on the elements of the $T_{1,N}$ matrix. I use them as a failsafe to ensure they are approximately the same ensuring there are no errors. They were not used here, I used them to test the accuracy of the probabilities."

# ╔═╡ 9e5e3db4-484b-4bab-ba30-0ee0a6768e73
#=
Transmission Probability 
Tr = (An/A1)^2 * (kn/k1)
Tr1 = (1/(T1N[1][1])^2)*(kn/k1)
Reflection Probability
R = (B1/A1)^2
R1 = (T1N[2][1]/T1N[1][1])^2
=#

# ╔═╡ 7c4ac936-cf83-40f1-9fa9-06654e154f38
md"Finally I defined the transmission as a function of energy, essentially identical to that used for the finite barrier."

# ╔═╡ 48830300-825d-4491-a3eb-746c6ba5ecfd
function waveFunction4(m::Float64,E::Float64,AN::Float64,BN::Float64,x1::Float64, xnsub1::Float64,N::Int64,v::Int64,V::Float64,RorIorTr::String)
	
	if v > N 
		return md"#### The specified region 'v' to create the wavefunction is more than the total number of regions"
	elseif v == 1
		return md"#### The specified region 'v' to create the wavefunction has to be greater than 1 to compare it to a region one less"
	end
	
	if N < 2
		println("You need to input at least two regions")
	else
		#Defining the different widths of the regions
		xnsub1 = xnsub1
		x1 = x1 
		
		deltax = (xnsub1 - x1)/(N-2)
		deltax = truncateToSf(deltax, 3)
		
		#Region boundaries
		BRi = (x1 - deltax):deltax:(xnsub1 + deltax)
		
		#Creation of regions
		dx = 0.01
		y = Int(round(((BRi[2] - BRi[1])/dx)) + 1)
		Ri = ones(Float64, N , y)
		
		for i in 1:(N)
			Ri[i, :] = BRi[i]:dx:(BRi[i+1] + 0.01*dx)	
		end
		
		#Defining the wavevectors at regions '1' 'v' 'v+1' and 'n' by taking the mean of the wavevectors at each point in each region
		kn = k(m,E,U(mean(Ri[N, :]),V)) 
		k1 = k(m,E,U(mean(Ri[1, :]),V))
		kv = k(m,E,U(mean(Ri[v, :]),V))
		kvminus1 = k(m,E,U(mean(Ri[v-1]),V))
		
		AV = 1.0
		BV = 1.0
		AVminus1 = 1.0
		BVminus1 = 1.0
		ANBN = [AN ; BN]
		T1N = [1.0 0.0; 0.0 1.0]
		TVN = [1.0 0.0; 0.0 1.0]
		TVminus1N = [1.0 0.0; 0.0 1.0]
		TE = ones(N-1)
		
		#Forms the T1N (transfer matrix) from region 1 to the nth region 
		for i in 1:(N-1)
			ki = k(m,E,U(mean(Ri[i, :]),V)) 
			kiplus1 = k(m,E,U(mean(Ri[i+1, :]),V)) 
			Tiplus1 = transferMatrix(ki,kiplus1,BRi[i+1])
			T1N = T1N * Tiplus1

			#Creates a transfer matrix for regions n to v, to help create a wavefunction at an arbitrary region v and the same for region v+1
			if i >= (v)
				TVN = TVN * Tiplus1
			end
			
			if i >= (v-1)
				TVminus1N = TVminus1N * Tiplus1
			end
			
			if i == (N-1)
				AVBV = TVN * ANBN
				AV = AVBV[1]
				BV = AVBV[2]
				
				AVminus1BVminus1 = TVminus1N * ANBN
				AVminus1 = AVminus1BVminus1[1]
				BVminus1 = AVminus1BVminus1[2]
			end
		end
		
		#Test to see if the array is defined without any 'NaN' and 'Inf' values
		#=
		if isArrayDefined(xv)
    		println("All values in the array are defined and not NaN or Inf.")
		else
    		println("The array contains NaN or Inf values.")
		end
		=#
		
		psivminus1 = AVminus1.*exp.(i*kvminus1.*Ri[v-1, :]) .+ BVminus1.*exp.(-i*kvminus1.*Ri[v-1, :])
		psiv = AV.*exp.(i*kv.*Ri[v, :]) .+ BV.*exp.(-i*kv.*Ri[v, :])
		
		A1B1 = T1N * ANBN
		A1 = A1B1[1]
		B1 = A1B1[2]
		
		#Transmission Probability 
		Tr = (AN/A1)^2 * (kn/k1)
		#Tr1 = (1/(T1N[1][1])^2)*(kn/k1)
		#Reflection Probability
		R = (B1/A1)^2
		#R1 = (T1N[2][1]/T1N[1][1])^2

		#To create the x (Energy(E)) and the y (Transmission(TE))
		dE = evtoJ(0.01)
		E1 = dE:dE:evtoJ(10.0)
		a = Int(E1[end]/dE)
		TE = ones(a)
		for i in 1:a
			kappa = k(m,E1[i],U(mean(Ri[2, :]), V))
			kprime = k(m,E1[i],U(mean(Ri[2, :]), V))
			ka = k(m,E1[i],0.0)
			
			if E1[i] < U(mean(Ri[2, :]), V)
				TE[i] = 1/(1 + (0.25)*(kappa/ka + ka/kappa)^(2)*(sinh(kappa*deltax))^2)
			elseif E1[i] == U(mean(Ri[2, :]), V)
				TE[i] = 1/(1 + (ka*deltax/2)^2)
			elseif E1[i] > U(mean(Ri[2, :]), V)
				TE[i] = 1/(1 + (0.25)*(kprime/ka - ka/kprime)^(2)*(sin(kprime*deltax))^2)
			end
		end

		if RorIorTr == "Real"
			#real
			plot(Ri[v-1, :] ,real(psivminus1) ,xlabel= "x",ylabel = "psi", label= "psivminus1" )
			plot!(Ri[v, :] ,real(psiv) ,xlabel= "x",ylabel = "psi", label= "psiv")
			plot!(title="Wave function (real) in regions $(v-1) & $(v)")
			
		elseif RorIorTr == "Imag"
			#imaginary
			plot!(Ri[v-1, :] ,imag(psivminus1) ,xlabel= "x",ylabel = "psi", label= "psivminus1" )
			plot(Ri[v, :] ,imag(psiv) ,xlabel= "x",ylabel = "psi", label= "psiv")
			plot!(title="Wave function (imaginary) in region $(v-1) & $(v)")

		elseif RorIorTr == "Tr"
			#Transmission probability against the energy
			plot(E1 ,TE ,xlabel= "Energy (J)",ylabel= "Tr", label= "T(E)")
			plot!(title="Transmission Probability")
			
		elseif RorIorTr != "Real" || "Imag" || "Tr"
			md"##### Input either the string 'Real' or 'Imag' to plot either the real or imaginary components of the wavefunction. And 'Tr' to plot the transmission as a function of energy"
		end
	end
end

# ╔═╡ 7afc1576-efe1-4650-9e6d-4779e26dc2a9
md"Syntax to call the function: $waveFunction4(m,E,A_N,B_N,x_1,x_{N-1},N,v,V,RorIorTr)$"

# ╔═╡ 11f11e8a-96fe-4c53-93bd-92cbba78b89b
md"Below $(4.1)$ I have plotted the five-piece approximation for a smooth repelling Gaussian potential. Refer to $[1]$ - Figure $14$ (Page $48$). In comparison to the literature, the plot does not have a slight dip after it levels off. The plot shows a steep slope immediately then levels off and does not follow the expected increasing then decreasing slope with a point of inflection at the value of $V = 5 eV$. This seems to follow the upper bound of the transmission probability for the finite barrier $(2.2)$ as opposed to the expected lower bound. "

# ╔═╡ db77f352-104a-4ac8-b61e-0f1088d8f51f
waveFunction4(m,evtoJ(5.5),1.0,0.0,-2.4,2.4,5,3,evtoJ(5.0) ,"Tr")

# ╔═╡ bca70dfb-7904-4996-89c1-34bf35d8adff
md"""
$$\tag{Figure 4.1}$$"""

# ╔═╡ e21ea900-95ff-43cc-a356-5d40b2cddc16
md"Here I have plotted the eleven-piece approximation for the same system. Refer to $[1]$ - Figure $15$ (Page $49$). This shouldn't have the dip as expected however similarly it seems to follow the upper bound as opposed to the lower bound of the finite barrier once again as expected."

# ╔═╡ 3997fabb-1d7d-4e65-a9a7-24066077b53c
waveFunction4(m,evtoJ(5.5),1.0,0.0,-2.4,2.4,11,3,evtoJ(5.0) ,"Tr")

# ╔═╡ d2eacf0b-48cb-4c7c-a18d-0bcd1c3c13f7
md"""
$$\tag{Figure 4.2}$$"""

# ╔═╡ 1d3343e9-bb3f-4091-9acc-132eaeb2907f
md"The five-piece approximation of an attracting Gaussian potential for when $V=-5eV$. Refer to $[1]$ - Figure $16$ (Page $50$). This shape follows the literature far more closely, yet it still does not have the expected dip. As it is an attracting barrier, the electron's energy is not relevant even if the energy of the electron is less than that of the potential barrier."

# ╔═╡ 6ac2c085-6bf3-49e7-a767-c723a531b076
waveFunction4(m,evtoJ(5.5),1.0,0.0,-2.4,2.4,5,3,evtoJ(-5.0) ,"Tr")

# ╔═╡ 4845d8e0-8abe-4cc2-ab5a-178988b64f05
md"""
$$\tag{Figure 4.3}$$"""

# ╔═╡ df6ba833-e5ca-476c-9c79-9139810996f4
md"This is the eleven-piece approximation for the same system. This follows the exact shape as expected. Without the dip as shown in the literature, refer to $[1]$ Figure $17$ (Page $50$)."

# ╔═╡ 9bf50566-f60c-440b-8168-08d7bb70525b
waveFunction4(m,evtoJ(5.1),1.0,0.0,-2.4,2.4,11,3,evtoJ(-5.0) ,"Tr")

# ╔═╡ de69dcd6-aa60-4589-a0a0-1309718c599f
md"""
$$\tag{Figure 4.4}$$"""

# ╔═╡ 1bec5464-ed80-44a1-ae34-c8d165833ad4
md"### Conclusions"

# ╔═╡ 07752783-ee50-4d82-b2c4-77c3a5ec394c
md"To conclude, the biggest inaccuracies were in the lack of a dip after levelling off for the five-piece approximations $(4.1$ & $4.3)$. Moreover, the two plots for the repelling barriers $(4.1$ & $4.2)$ did not pass $Tr = 0.5$ when $E = V$ (which in this case, $V = 5eV$). This should be the case as for a repelling barrier its transmission can only be greater than its reflection when the energy of the particle is greater than the potential barrier."

# ╔═╡ e6384b83-bc4c-4651-92f8-60b9bf106730
md"However, for the eleven-piece attracting Gaussian potential the plot is modelled exactly how the electron would behave for a range of energies. There were issues with the rounding and accuracy of values in Julia which made it difficult to define the boundary values and the regions. I also chose to take the variable potentials from the centre of each region due to how it was modelled in reference $[1]$."

# ╔═╡ e21a8b0b-43c5-4b47-8e85-ca4d05495f43
md"The lack of accuracy in the transmission probability graphs for the arbitrary number of regions (Section $4$) comes from the analytical nature of the equations (shown in Section $2$). There are many limitations, analytical solutions tend not to work well for more complex systems. Take the double barrier example (Section $3$). The equations I defined in Section $2$ assume that it is a single rectangular potential barrier, where the potential on the left of the barrier is equal to the potential on the right $V_L = V_R$. If $U_1$ and $U_2$ are different it would be a poor approximation."

# ╔═╡ 25ae9634-d33a-46af-87be-c9122559877a
md"This is further exasperated for the Gaussian potential (Section $4$) when using rectangular barriers as an approximation (five and eleven rectangle approximation). If I were to do this again I would tailor the solutions for each system. I would also select approximations which more closely follow the variable potential (i.e. a far greater number of regions which slowly get closer to the shape of the potential."

# ╔═╡ a233fd2a-7832-4261-b503-d9f186e3d0e9
md"When attempting to plot the wavefunctions for the arbitrary number of regions in a Gaussian shape [Section $4$], numerous issues arose. As was shown in Section $3$ the width of the barrier can control how much decay occurs. For a Gaussian barrier, where the
width of the potential gets smaller as you go up, it was difficult to make sense of the output. Moreover, the equation given did not adhere well as it was plotted using approximations, meaning changing the parameters did not have a clear meaning to what was being changed within the system. Ideally, a better method would go down the analytical route as opposed to the numerical method obtained from $[1]$. Naturally, it is not a guarantee that an analytical solution would be possible, or feasible, hence a better numerical method may be useful. "

# ╔═╡ b91b2959-79fb-45fb-8445-74098925823c
md"### References"

# ╔═╡ ec095b61-dd94-4c06-8e06-4816c4c710c8
md"##### $[1]$ = Robert Gilmore (2013) [*Elementary Quantum Mechanics in One Dimension*](http://www.physics.drexel.edu/~bob/NewBook.pdf)"

# ╔═╡ f2a57229-4164-49ef-9d94-9eca8dbb45a1
md"##### $[2]$ = Wikipedia Contributers (2023) [*Resonances in scattering from potentials*](https://en.wikipedia.org/wiki/Resonances_in_scattering_from_potentials)"

# ╔═╡ 135d540e-9f7a-4a17-8d73-a031c56631a7
md"##### $[3]$ = (2013) [*4.4 The Potential Step*](https://www.hansrajcollege.ac.in/hCPanel/uploads/elearning/elearning_document/study_material_from_quantum_mechanics_by_N_zettilli.pdf)"

# ╔═╡ dd0e23ec-6e07-4dda-aea1-68f04d67c62d
md"##### $[4]$ = Nikola Opacak, Vitomir Milanovic and Jelena Radovanovic (2017) [*Transmission singularities in resonant electron tunnelling through double complex potential barrier*](https://arxiv.org/ftp/arxiv/papers/1704/1704.08890.pdf)"

# ╔═╡ dca9b00e-787e-40d0-bbc5-795033965c53
md"##### $[5]$ = Mark Greenaway (2023) [*The transfer matrix method*](https://learn.lboro.ac.uk/mod/resource/view.php?id=1440719)"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Plots = "~1.39.0"
PlutoUI = "~0.7.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "27e6ac469491e741848d1f8b2e233b959955f0a0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "cd67fc487743b2f0fd4380d4cbd3a24660d0eec8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.3"

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
deps = ["UUIDs"]
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "5372dbbf8f0bdb8c700db5367132925c0771ef7e"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

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
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

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
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

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

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

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
git-tree-sha1 = "5eab648309e2e060198b45820af1a37182de3cce"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "9fb0b890adab1c0a4a475d4210d51f228bfc250d"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.6"

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
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

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
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

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
version = "2022.10.11"

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
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

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
version = "1.9.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

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
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "7c29f0e8c575428bd84dc3c72ece5178caa67336"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.2+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
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
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

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
git-tree-sha1 = "5165dfb9fd131cf0c6957a3a7605dede376e7b63"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

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
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

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
git-tree-sha1 = "49cbf7c74fafaed4c529d47d48c8f7da6a19eb75"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.1"
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
git-tree-sha1 = "a72d22c7e13fe2de562feda8645aa134712a87ee"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.17.0"

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
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "24b81b59bd35b3c42ab84fa589086e19be919916"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.11.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cf2c7de82431ca6f39250d2fc4aacd0daa1675c0"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.4+0"

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
version = "1.2.13+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "47cf33e62e138b920039e8ff9f9841aafe1b733e"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.35.1+0"

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
version = "5.8.0+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

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
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

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
# ╟─56c20ce8-e208-46db-80ae-ecbe7df7ad9e
# ╟─cead5407-18d3-4ca6-b98e-e852308e1937
# ╟─92fac372-efd4-4dc9-80f0-fded906ec865
# ╟─9bf9ca86-e6ca-43b8-bac4-bea979cf728b
# ╟─d50eff2a-55f2-4681-95c6-31d076f1c709
# ╟─a841a6ea-ff03-474c-a942-4c3cb74257f5
# ╟─2a2068d1-6419-44ef-89e4-2dc192991b9e
# ╟─e25552ec-811f-48b9-a1f6-9615ae05e554
# ╟─0ba047a4-b6f9-4ca3-9597-c3e2f7418bfb
# ╟─4aeaeaea-10ae-4d29-aba6-3bb47c5ebdf0
# ╟─ff4f1c98-5b3d-4823-b900-d9ba5dde07ef
# ╟─896d5022-0e2d-4610-89f2-fda508b6d41c
# ╟─9f4f6842-89a3-4cc9-955a-60b484e1a8c3
# ╟─ac810460-6acf-4d2c-bd2e-6d97cebefed7
# ╟─c85bad5c-0491-429f-b02e-998c59da6b23
# ╟─332f42be-7674-4717-b2c3-e522ac0970f3
# ╟─fe63d788-3a0a-4748-b577-00e67e031576
# ╠═554b9b10-6f27-11ee-2678-f10719029989
# ╟─8af54e36-f964-426b-b790-2ad1a0f1ac94
# ╠═82789e60-7dc4-4452-b328-57c9e7cdab22
# ╟─c264e0ac-f50c-4fce-acf2-421bf9d8f01e
# ╠═26e91094-9a06-4a9b-a58e-d9510c07af8c
# ╟─9195cd37-9de5-4797-96ba-0e2be3002dad
# ╠═68cc8a11-f3d6-4ac7-baf0-b55cec5e3ee7
# ╟─0941fa2c-3249-4bc2-9885-899effc0e3fe
# ╠═c207319c-dee6-4a87-b968-7ac8c86f047d
# ╟─49ae3345-2d58-4424-877f-6f5f1de44546
# ╠═dc480e3f-828a-4398-a944-68064c708443
# ╟─81284779-be86-46d4-9e68-817225790be6
# ╟─309e640a-0fef-4209-893f-3c6068d0bfd2
# ╟─bb9d49af-1772-4c85-ac90-29e19f533c86
# ╟─839ad551-f112-4899-a0dd-ba63fdfe5788
# ╟─da9904d7-3ab6-4cb5-b66a-f2bc94c0585b
# ╠═0a135be0-38ce-45d7-80d5-e10ec09f81ba
# ╟─d9bb8263-5b0d-485a-9e6a-bdb02f3c8443
# ╟─e2e9f666-f840-47d8-9a73-04157395d822
# ╟─3a4dba24-705b-4329-a8b4-6a31e2d5eb87
# ╟─136e351f-f373-43a9-a799-3bb4e820d49e
# ╟─5728e695-b021-4ba8-b2b0-d1c1e4592de8
# ╟─d9977075-05ea-4625-9002-7f2a767544c7
# ╟─8ac53ff1-1b91-44a8-ace2-81134f67743b
# ╟─6c102db4-74e7-42fc-b310-03ec11d53642
# ╟─11b9d64a-f4f1-4571-ad14-2a8959a81ff0
# ╟─64618b84-5fa2-45f1-a697-57cc3fa8378c
# ╟─0a509173-c9ff-431c-a5ff-2f50393465da
# ╟─960ce184-befd-41f0-8f13-84ea7c3cd007
# ╟─9d59944a-5a39-4157-beb0-d2a4d6c5dbfd
# ╟─4aae4592-19c1-41dd-aadf-2c96e68df96f
# ╟─583a7e06-2a9e-481e-b21d-923f4012bf9a
# ╟─c8c81211-f8c5-412b-a297-14ae8f12359e
# ╟─22edd46c-a68b-4e6b-a419-7940ea18efc0
# ╟─c5918e70-99f8-488f-996e-7153dc6c12e2
# ╠═e25e8620-8713-401e-83de-083fc269763c
# ╠═f936738a-42d9-4e9b-bcbb-86a252998c79
# ╟─71c797bf-abd2-4791-a367-a52b7c083b52
# ╟─1b6d2638-4a48-4eb3-a95e-5ecdab18bc73
# ╟─8757956d-6797-40ba-ae3a-d087fba12999
# ╟─71f4c76b-fc1e-4828-86e7-db49fbd311ee
# ╟─ba34511b-367f-49d6-9cee-3da45070ad75
# ╟─6cda3ee2-eebf-4330-a713-d1ad8f6e6b7f
# ╟─4007531f-5c53-47bc-8400-5d2c8fd7786c
# ╟─6e0fc9bf-8abb-444b-9653-af7a97e61f2b
# ╟─3dca1e0b-fbe5-42e6-87e2-81d283491829
# ╟─2726ed0f-d546-43a0-b9f2-3c1728054126
# ╟─39e540d3-d039-4f09-8469-252055162ae9
# ╠═059713dc-c14b-4261-822f-84650f595775
# ╟─6d9824c2-a1d3-4168-933b-d8a8b7a5a7f9
# ╟─87154415-4a49-474e-bb21-a5424d4de711
# ╠═0e676326-bc3e-4f73-9d62-99d3e02169c7
# ╟─269cf590-06cc-402d-8c07-77c85055fbdf
# ╟─5dcd0bf5-b13f-4ea3-92c7-9c66f2fac085
# ╟─d393922b-a0e9-4524-b71f-0310369a0e4d
# ╟─8ac27cfd-9c9b-4f81-b6ba-adcae7e2e13e
# ╟─a03ef348-17b8-4058-9fea-5e36c3f9a1ec
# ╟─9335ffb6-098e-44db-8ac4-554b3a3495dc
# ╟─4794b945-9ac0-465c-ad65-713331acf94e
# ╟─9be1351c-8fd5-412f-b852-6226c14b5f7d
# ╟─49dd843d-9791-427b-842b-2fb1d3642486
# ╟─5674456b-a714-4b33-aa17-ec5656328e30
# ╟─813dcd90-6bad-4811-900d-8f8f3ea67fd0
# ╟─af6c21bc-fd3f-46e7-830d-58fae6dab7d3
# ╟─39cbf241-c4af-4411-aaba-532b35eeac3a
# ╟─a9006f30-a9c7-415e-bfa5-e0212fbebd11
# ╟─7952d09c-586e-48db-a577-65375cdffc8a
# ╟─80426867-5e2f-41da-851d-1aad5c9126a9
# ╟─cfad5c07-c4a2-45b9-a4b6-89beaadc8cfd
# ╟─877bb084-6e7c-46f4-b2ef-c8ba682bcbd5
# ╟─8d6ee29a-4a8f-461c-8722-527a5508efe1
# ╠═3a3cc777-957c-480a-95e2-ff3003242847
# ╟─269a6d76-66b5-4230-87fc-19c92785c07b
# ╠═3ce84e57-8be5-4134-a171-2539ad705c1b
# ╟─b2e27d7c-437a-46d1-accf-d6ed1fa6df5f
# ╟─e9154db9-cf7b-49eb-8b7f-4fd1dc0de6b1
# ╠═c72d20b4-5617-4b7e-b464-9cefab50cd86
# ╟─2b4be712-5bdb-45d2-b2b3-7d44c09a9ab6
# ╠═f677a1e8-1b0d-404e-9ee2-2aac2b84047a
# ╟─485a8dab-575d-4ba0-8cc6-3b3d3d37ee44
# ╠═77af237b-73c2-46ae-85c5-e5bed1dac5ed
# ╟─4f3c0619-03de-46af-9791-6582ce745f79
# ╠═8a28ed21-93d6-43a7-9829-9981614aadb1
# ╟─98596d44-87cf-451b-968d-450d9c0c9511
# ╠═9e5e3db4-484b-4bab-ba30-0ee0a6768e73
# ╟─7c4ac936-cf83-40f1-9fa9-06654e154f38
# ╠═48830300-825d-4491-a3eb-746c6ba5ecfd
# ╟─7afc1576-efe1-4650-9e6d-4779e26dc2a9
# ╟─11f11e8a-96fe-4c53-93bd-92cbba78b89b
# ╟─db77f352-104a-4ac8-b61e-0f1088d8f51f
# ╟─bca70dfb-7904-4996-89c1-34bf35d8adff
# ╟─e21ea900-95ff-43cc-a356-5d40b2cddc16
# ╟─3997fabb-1d7d-4e65-a9a7-24066077b53c
# ╟─d2eacf0b-48cb-4c7c-a18d-0bcd1c3c13f7
# ╟─1d3343e9-bb3f-4091-9acc-132eaeb2907f
# ╟─6ac2c085-6bf3-49e7-a767-c723a531b076
# ╟─4845d8e0-8abe-4cc2-ab5a-178988b64f05
# ╟─df6ba833-e5ca-476c-9c79-9139810996f4
# ╟─9bf50566-f60c-440b-8168-08d7bb70525b
# ╟─de69dcd6-aa60-4589-a0a0-1309718c599f
# ╟─1bec5464-ed80-44a1-ae34-c8d165833ad4
# ╟─07752783-ee50-4d82-b2c4-77c3a5ec394c
# ╟─e6384b83-bc4c-4651-92f8-60b9bf106730
# ╟─e21a8b0b-43c5-4b47-8e85-ca4d05495f43
# ╟─25ae9634-d33a-46af-87be-c9122559877a
# ╟─a233fd2a-7832-4261-b503-d9f186e3d0e9
# ╟─b91b2959-79fb-45fb-8445-74098925823c
# ╟─ec095b61-dd94-4c06-8e06-4816c4c710c8
# ╟─f2a57229-4164-49ef-9d94-9eca8dbb45a1
# ╟─135d540e-9f7a-4a17-8d73-a031c56631a7
# ╟─dd0e23ec-6e07-4dda-aea1-68f04d67c62d
# ╟─dca9b00e-787e-40d0-bbc5-795033965c53
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
