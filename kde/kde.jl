using KernelDensity
using Printf

SpikeTime_typed = convert(Array{Float64,1}, SpikeTime)
U = kde(SpikeTime_typed,bandwidth=10)

p1 = scatter(SpikeTime, SpikeNeuron, xlims=(0,1000), label="", title="Raster Plot", ylabel="Neuron")
p2 = histogram(SpikeTime, bins=(range(0,1000,step=10)), xlims=(0,1000), label="", title="number of spikes", ylabel="number of spikes")
p3 = plot(U.x, U.density, xlims=(0, 1000), label="", title="kernel density", xlabel="Time[ms]")
plot(p1, p2, p3, layout=(3,1), size=(900, 600))
savefig("figure/result.png")

# integration
dt = U.x[2] - U.x[1]
area = 0
for i=1:length(U.x)
    global area = area + U.density[i]*dt
end
@printf("area:%f", area)
