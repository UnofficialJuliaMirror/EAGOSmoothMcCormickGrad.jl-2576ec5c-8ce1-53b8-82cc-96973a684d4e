#workspace()

using EAGOSmoothMcCormickGrad
using IntervalArithmetic
using StaticArrays
#using DataFrames
#using CSV
using Plots
using PyPlot

# Plots the differentiable McCormick relaxations and
# original function x*y over [0,200] by [0,400], respectively
EAGOSmoothMcCormickGrad.set_diff_relax(1)

xl = 100
xu = 200
yl = -400
yu = -100
intv = 10
nx = (xu-xl)/intv+1
ny = (yu-yl)/intv+1
Vals = zeros(Float64,nx*ny,9)
x_grid = zeros(Float64,nx,ny)
y_grid = zeros(Float64,nx,ny)
values = zeros(Float64,nx,ny)
ccval = zeros(Float64,nx,ny)
cvval = zeros(Float64,nx,ny)
ccgrad1 = zeros(Float64,nx,ny)
ccgrad2 = zeros(Float64,nx,ny)
cvgrad1 = zeros(Float64,nx,ny)
cvgrad2 = zeros(Float64,nx,ny)

count1 = 1
for i=yl:intv:yu
    count2 = 1
    for j=xl:intv:xu
        temp1 = SMCg{2,Float64}(Float64(j),Float64(j),seed_g(Float64,1,2),
                                 seed_g(Float64,1,2),Interval(xl,xu),
                                 false, 1.0*[Interval(xl,xu),Interval(yl,yu)],[(xu-xl)/2.0,(yu-yl)/2.0])
        temp2 = SMCg{2,Float64}(Float64(i),Float64(i),seed_g(Float64,2,2),
                                 seed_g(Float64,2,2),Interval(yl,yu),
                                 false, 1.0*[Interval(xl,xu),Interval(yl,yu)],[(xu-xl)/2.0,(yu-yl)/2.0])
        x_grid[count2,count1] = j # x
        y_grid[count2,count1] = i # y
        values[count2,count1] = (j)*(i) #x*y
        temp_MC = temp1*temp2 # relax x*y

        ccval[count2,count1] = temp_MC.cc
        cvval[count2,count1] = temp_MC.cv
        ccgrad1[count2,count1] = temp_MC.cc_grad[1]
        ccgrad2[count2,count1] = temp_MC.cc_grad[2]
        cvgrad1[count2,count1] = temp_MC.cv_grad[1]
        cvgrad2[count2,count1] = temp_MC.cv_grad[2]
        count2 += 1
    end
    count1 += 1
end

#plot(layout = 4)
surface(x_grid,y_grid,cvval, xlab = "x", ylab = "y", zlab = "z")
p1 = surface(x_grid,y_grid,cvgrad1, xlab = "x", ylab = "y", zlab = "z")
#p2 = surface(x_grid,y_grid,cvgrad2, xlab = "x", ylab = "y", zlab = "z")
#pout = plot(p1,p2,layout = 2)
#pout = plot(p1,layout = 1)
#show(pout)
gui()
#=
subplot(221)
surface(x_grid,y_grid,ccval)
subplot(222)
surface(x_grid,y_grid,ccval)
subplot(223)
surface(x_grid,y_grid,ccval)
subplot(224)
surface(x_grid,y_grid,ccval)
suptitle("2x2 Subplot")
=#
#gui()
#DF = convert(DataFrame,Vals)
#CSV.write("/home/mewilhel/Desktop/MultiplicationPlot3.csv", DF)


seed1 = seed_g(Float64,1,2)
seed2 = seed_g(Float64,2,2)
x1 = SMCg{2,Float64}(0.0,0.0,seed1,seed1,Interval(-200.0,200.0),false, [Interval(-200.0,200.0),Interval(0.0,400.0)],[0.0,200.0])
y1 = SMCg{2,Float64}(200.0,200.0,seed2,seed2,Interval(0.0,400.0),false, [Interval(-200.0,200.0),Interval(0.0,400.0)],[0.0,200.0])
z1 = x1*y1
println("z1.cc: $(z1.cc)")
println("z1.cv: $(z1.cv)")

x2 = SMCg{2,Float64}(170.0,170.0,seed1,seed1,Interval(100.0,240.0),false, [Interval(100.0,240.0),Interval(100.0,400.0)],[170.0,250.0])
y2 = SMCg{2,Float64}(250.0,250.0,seed2,seed2,Interval(100.0,400.0),false, [Interval(100.0,240.0),Interval(100.0,400.0)],[170.0,250.0])
z2 = x2*y2
println("z2.cc: $(z2.cc)")
println("z2.cv: $(z2.cv)")

x3 = SMCg{2,Float64}(-200.0,-200.0,seed1,seed1,Interval(-300.0,-100.0),false, [Interval(-300.0,-100.0),Interval(-400.0,-200.0)],[-200.0,-300.0])
y3 = SMCg{2,Float64}(-300.0,-300.0,seed2,seed2,Interval(-400.0,-200.0),false, [Interval(-300.0,-100.0),Interval(-400.0,-200.0)],[-200.0,-300.0])
z3 = x3*y3
println("z3.cc: $(z3.cc)")
println("z3.cv: $(z3.cv)")

x4 = SMCg{2,Float64}(150.0,150.0,seed1,seed1,Interval(100.0,200.0),false, [Interval(100.0,200.0),Interval(-500.0,-100.0)],[150.0,-300.0])
y4 = SMCg{2,Float64}(-250.0,-250.0,seed2,seed2,Interval(-500.0,-100.0),false, [Interval(100.0,200.0),Interval(-500.0,-100.0)],[150.0,-300.0])
z4 = x4*y4
println("z4.cc: $(z4.cc)")
println("z4.cv: $(z4.cv)")

x5 = SMCg{2,Float64}(-150.0,-150.0,seed1,seed1,Interval(-200.0,-100.0),false, [Interval(-200.0,-100.0),Interval(200.0,400.0)],[-150.0,300.0])
y5 = SMCg{2,Float64}(300.0,300.0,seed2,seed2,Interval(200.0,400.0),false, [Interval(-200.0,-100.0),Interval(200.0,400.0)],[-150.0,300.0])
z5 = x5*y5
println("z5.cc: $(z5.cc)")
println("z5.cv: $(z5.cv)")
