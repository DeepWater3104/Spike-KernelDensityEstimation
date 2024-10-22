using KernelDensity
using Printf

SpikeTime_typed = convert(Array{Float64,1}, SpikeTime)
U = kde(SpikeTime_typed,bandwidth=10)

num_spikes = length(SpikeTime)
R_b = U.density * num_spikes / N

p1 = scatter(SpikeTime, SpikeNeuron, xlims=(0,1000), label="", ylabel="Neuron", ms=3)
p2 = histogram(SpikeTime, bins=(range(0,1000,step=10)), xlims=(0,1000), label="", ylabel="number of spikes")
p3 = plot(U.x, R_b, xlims=(0, 1000), label="", ylabel="IPBR[Hz]", xlabel="Time[ms]")
plot(p1, p2, p3, layout=(3,1), size=(800, 600))
savefig("figure/result.png")

# integration
dt = U.x[2] - U.x[1]
area = 0
for i=1:length(U.x)
    global area = area + U.density[i]*dt
end
@printf("area:%f", area)
