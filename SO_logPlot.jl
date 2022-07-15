cd(@__DIR__)

using PyPlot
using CSV, DataFrames
using Polynomials
using LaTeXStrings
using Printf

function loglogplot(DF, size, title)
	fig, ax = plt.subplots(figsize = (7, 7))

	ax.set_xscale("log");	# Set logarithmic scale
	ax.set_yscale("log");
	ax.set_xlabel("N", fontsize = 14)
	ax.set_ylabel("Time [sec]", fontsize = 14)
	ax.set_title("Execution time as function of system size - "*title, fontweight ="bold")
	ax.set_xlim([90,size+10])
	ax.tick_params(labelsize = 12)

	# Learning, direct
	ax.scatter(DF[!, :N], DF[!, :L_direct], label="Learning, direct", s=10, marker="o",c="firebrick")

	logX = log.(DF[!, :N])
	logY = log.(DF[!, :L_direct])
	f = fit(logX, logY, 1)
	intercept, slope = f.coeffs
	y_fit = exp(intercept) .* (DF[!, :N]).^ slope

	n = @sprintf("%.3E", exp(intercept))
	m = @sprintf("%.3f", slope)
	label_txt = raw"y(x) = " * n * raw"*" * raw"$x^{" * m * raw"}$"
	plot(DF[!, :N], y_fit, label=label_txt, color="firebrick")

	# Learning, on-the-fly
	ax.scatter(DF[!, :N], DF[!, :L_onthefly], label="Learning, on-the-fly", s=10, marker="s",c="darkblue")

	logX = log.(DF[!, :N])
	logY = log.(DF[!, :L_onthefly])
	f = fit(logX, logY, 1) 							# log(y) = mlog(x) + n
	intercept, slope = f.coeffs 					# n, m
	y_fit = exp(intercept) .* (DF[!, :N]).^ slope	# y(x)=exp(n)*x^m

	n = @sprintf("%.3E", exp(intercept))
	m = @sprintf("%.3f", slope)
	label_txt = raw"y(x) = " * n * raw"*" * raw"$x^{" * m * raw"}$"
	plot(DF[!, :N], y_fit, label=label_txt, color="darkblue")


	# No Learning
	ax.scatter(DF[!, :N], DF[!, :NL], label="No Learning", s=10, marker="^",c="grey")

	logX = log.(DF[!, :N])
	logY = log.(DF[!, :NL])
	f = fit(logX, logY, 1) 							# log(y) = mlog(x) + n
	intercept, slope = f.coeffs 					# n, m
	y_fit = exp(intercept) .* (DF[!, :N]).^ slope	# y(x)=exp(n)*x^m

	n = @sprintf("%.3E", exp(intercept))
	m = @sprintf("%.3f", slope)
	label_txt = raw"y(x) = " * n * raw"*" * raw"$x^{" * m * raw"}$"
	plot(DF[!, :N], y_fit, label=label_txt, color="grey")

	legend()
	ax.legend(loc="best", fontsize=12)
end

path = "./times/"

filename1 = "times_array_10resets.csv"
df1 = CSV.read(path*filename1, DataFrame, header=1)
loglogplot(df1, 10000, "10 resets")

fig_name = "logPlot_N10000"
savefig(path * fig_name * ".png")
# savefig(path * fig_name * ".svg")

filename2 = "times_array_1000resets.csv"
df2 = CSV.read(path*filename2, DataFrame, header=1)
loglogplot(df2, 2000, "1000 resets")

fig_name = "logPlot_N2000"
savefig(path * fig_name * ".png")
# savefig(path * fig_name * ".svg")

##. Separate fits by range - 10 resets

fig, ax = plt.subplots(figsize = (7, 7))

ax.set_xscale("log");	# Set logarithmic scale
ax.set_yscale("log");
ax.set_xlabel("N", fontsize = 14)
ax.set_ylabel("Time [sec]", fontsize = 14)
ax.set_title("Execution time as function of system size - 10 resets", fontweight ="bold")
ax.set_xlim([90,10010])
ax.tick_params(labelsize = 12)

# Learning, direct
ax.scatter(df1[!, :N], df1[!, :L_direct], label="Learning, direct", s=10, marker="o",c="firebrick")

N_max = 1100

logX = log.(df1[df1.N .<= N_max , :N])
logY = log.(df1[df1.N .<= N_max , :L_direct])

f = fit(logX, logY, 1)
intercept, slope = f.coeffs
y_fit = exp(intercept) .* (df1[df1.N .<= N_max, :N]).^ slope

n = @sprintf("%.3E", exp(intercept))
m = @sprintf("%.3f", slope)
label_txt = raw"y(x) = " * n * raw"*" * raw"$x^{" * m * raw"}$"
plot(df1[df1.N .<= N_max, :N], y_fit, label=label_txt)

N_min = 1100
N_max = 10000

logX = log.(df1[N_min .< df1.N .<= N_max , :N])
logY = log.(df1[N_min .< df1.N .<= N_max , :L_direct])

f = fit(logX, logY, 1)
intercept, slope = f.coeffs
y_fit = exp(intercept) .* (df1[N_min .< df1.N .<= N_max, :N]).^ slope

n = @sprintf("%.3E", exp(intercept))
m = @sprintf("%.3f", slope)
label_txt = raw"y(x) = " * n * raw"*" * raw"$x^{" * m * raw"}$"
plot(df1[N_min .< df1.N .<= N_max, :N], y_fit, label=label_txt)

# Learning, on-the-fly
ax.scatter(df1[!, :N], df1[!, :L_onthefly], label="Learning, on-the-fly", s=10, marker="s",c="darkblue")

logX = log.(df1[!, :N])
logY = log.(df1[!, :L_onthefly])
f = fit(logX, logY, 1) 							# log(y) = mlog(x) + n
intercept, slope = f.coeffs 					# n, m
y_fit = exp(intercept) .* (df1[!, :N]).^ slope	# y(x)=exp(n)*x^m

n = @sprintf("%.3E", exp(intercept))
m = @sprintf("%.3f", slope)
label_txt = raw"y(x) = " * n * raw"*" * raw"$x^{" * m * raw"}$"
plot(df1[!, :N], y_fit, label=label_txt, color="darkblue")


# No Learning
ax.scatter(df1[!, :N], df1[!, :NL], label="No Learning", s=10, marker="^",c="grey")

logX = log.(df1[!, :N])
logY = log.(df1[!, :NL])
f = fit(logX, logY, 1) 							# log(y) = mlog(x) + n
intercept, slope = f.coeffs 					# n, m
y_fit = exp(intercept) .* (df1[!, :N]).^ slope	# y(x)=exp(n)*x^m

n = @sprintf("%.3E", exp(intercept))
m = @sprintf("%.3f", slope)
label_txt = raw"y(x) = " * n * raw"*" * raw"$x^{" * m * raw"}$"
plot(df1[!, :N], y_fit, label=label_txt, color="grey")

legend()
ax.legend(loc="best", fontsize=12)
fig_name = "logPlot_N10000_fits"
savefig(path * fig_name * ".png")
savefig(path * fig_name * ".svg")

##. Separate fits by range - 1000 resets

fig, ax = plt.subplots(figsize = (7, 7))

ax.set_xscale("log");	# Set logarithmic scale
ax.set_yscale("log");
ax.set_xlabel("N", fontsize = 14)
ax.set_ylabel("Time [sec]", fontsize = 14)
ax.set_title("Execution time as function of system size - 1000 resets", fontweight ="bold")
ax.set_xlim([90,2010])
ax.tick_params(labelsize = 12)

# Learning, direct
ax.scatter(df2[!, :N], df2[!, :L_direct], label="Learning, direct", s=10, marker="o",c="firebrick")

N_max = 1200

logX = log.(df2[df2.N .<= N_max , :N])
logY = log.(df2[df2.N .<= N_max , :L_direct])

f = fit(logX, logY, 1)
intercept, slope = f.coeffs
y_fit = exp(intercept) .* (df2[df2.N .<= N_max, :N]).^ slope

n = @sprintf("%.3E", exp(intercept))
m = @sprintf("%.3f", slope)
label_txt = raw"y(x) = " * n * raw"*" * raw"$x^{" * m * raw"}$"
plot(df2[df2.N .<= N_max, :N], y_fit, label=label_txt)

N_min = 1200
N_max = 2000

logX = log.(df2[N_min .< df2.N .<= N_max , :N])
logY = log.(df2[N_min .< df2.N .<= N_max , :L_direct])

f = fit(logX, logY, 1)
intercept, slope = f.coeffs
y_fit = exp(intercept) .* (df2[N_min .< df2.N .<= N_max, :N]).^ slope

n = @sprintf("%.3E", exp(intercept))
m = @sprintf("%.3f", slope)
label_txt = raw"y(x) = " * n * raw"*" * raw"$x^{" * m * raw"}$"
plot(df2[N_min .< df2.N .<= N_max, :N], y_fit, label=label_txt)

# Learning, on-the-fly
ax.scatter(df2[!, :N], df2[!, :L_onthefly], label="Learning, on-the-fly", s=10, marker="s",c="darkblue")

logX = log.(df2[!, :N])
logY = log.(df2[!, :L_onthefly])
f = fit(logX, logY, 1) 							# log(y) = mlog(x) + n
intercept, slope = f.coeffs 					# n, m
y_fit = exp(intercept) .* (df2[!, :N]).^ slope	# y(x)=exp(n)*x^m

n = @sprintf("%.3E", exp(intercept))
m = @sprintf("%.3f", slope)
label_txt = raw"y(x) = " * n * raw"*" * raw"$x^{" * m * raw"}$"
plot(df2[!, :N], y_fit, label=label_txt, color="darkblue")


# No Learning
ax.scatter(df2[!, :N], df2[!, :NL], label="No Learning", s=10, marker="^",c="grey")

logX = log.(df2[!, :N])
logY = log.(df2[!, :NL])
f = fit(logX, logY, 1) 							# log(y) = mlog(x) + n
intercept, slope = f.coeffs 					# n, m
y_fit = exp(intercept) .* (df2[!, :N]).^ slope	# y(x)=exp(n)*x^m

n = @sprintf("%.3E", exp(intercept))
m = @sprintf("%.3f", slope)
label_txt = raw"y(x) = " * n * raw"*" * raw"$x^{" * m * raw"}$"
plot(df2[!, :N], y_fit, label=label_txt, color="grey")

legend()
ax.legend(loc="best", fontsize=12)
fig_name = "logPlot_N2000_fits"
savefig(path * fig_name * ".png")
savefig(path * fig_name * ".svg")
