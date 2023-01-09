# function to simulate on replicate of the process
simulate.resistance = function(generations,resident,mutation,r,s,u,threshold)
{
  
  for(j in 1:generations)
  {
    res.exp = resident*(1-r)
    mut.exp = mutation*(1+s)
    
    resident = rpois(1,res.exp)
    mutation = rpois(1,mut.exp)
    
    new.mutations = rbinom(1,resident,u)
    
    resident = resident - new.mutations
    mutation = mutation + new.mutations
    
    if(resident == 0 & mutation == 0)
    {
      return(0)
    }
    
    if(mutation > threshold)
    {
      return(1)
    }
  }
  
  return(NA)
}

# function to simulate average of multiple replicate of the process
get.resc.prop = function(gen,wt,mut,stress,sel.coefficient,mut.rate,rep)
{
  p.rescue.all = numeric(rep)
  for(i in 1:rep)
  {
    p.rescue.all[i] = simulate.resistance(generations = gen,resident = wt,
                                          mutation = mut,r = stress,s = sel.coefficient,
                                          u = mut.rate,threshold = 1/sel.coefficient)
  }
  return(mean(p.rescue.all,na.rm=T))
}

# simulate rep replications to get the rescue probability
rep = 100
# number of wild type individuals at beginning
wt = 1000
# number of mutations initially 
mut = 0

# set parameter values 
mut.rate = 0.0001
stress = 0.2
sel.coefficient = 0.1

# set maximum number of generations to simulate
gen = 1000

# estimate resistance probability for one set of parameters:
p = get.resc.prop(gen  = gen,wt = wt,mut = mut,stress = stress,
                  sel.coefficient = sel.coefficient,mut.rate = mut.rate,
                  rep = rep)


#-----
# A example


r.vec = c(0.0001,0.1,0.2,0.3,0.4,0.5)
n = length(r.vec)
p.resc = numeric(n)
rep = 1000

for(i in 1:n)
{
  r = r.vec[i]
  p.resc[i] = get.resc.prop(gen  = gen,wt = wt,mut = mut,stress = r,
                                sel.coefficient = sel.coefficient,
                                mut.rate = mut.rate,
                                rep = rep)
}

plot(r.vec,p.resc,pch=16,
     xlab = "env. stress (r)",ylab = "Prob. resistance")

r.vals = seq(0,1,by=0.01)
theory = 1 - exp(- 2 * wt *mut.rate* sel.coefficient/r.vals )
lines(r.vals,theory)

legend("topright",c("simulations","theory"),col="black",lwd=c(1,1)
       ,pch=c(16,-1),bty="n")

#-----
#2b question : time-dependent r
#Intial parameters settings, feel free to modify them at your needs.
time = c(0:50)
n = length(time)
r.vec = numeric(n)
p.resc = numeric(n)
pop.size = numeric(n+1)
pop.size[1]= 1000

#We are going to alter r in a range between 0 to 0.2 in these simulations
#for the sake of seeing differences between the situations

#1. What happens if stress is increasing over time
for(i in 1:n) {
  r = 0.01 + 0.01*(time[i]/4)
  r.vec[i] = r
  pop.size[i+1] = pop.size[i]*exp(-r)
  p.resc[i] = get.resc.prop(gen  = gen,wt = wt,mut = mut,stress = r,
                            sel.coefficient = sel.coefficient,
                            mut.rate = mut.rate,
                            rep = rep)
}

par(mfrow = c(3,1))

#Plotting 3 figures
#Environmental stress over time
plot(time,r.vec,pch = 15,
     xlab = "generation", ylab = "env.stress (r)")
#Probability of resistance/rescue over time
plot(time,p.resc,pch=17,
     xlab = "generation", ylab = "Prob. resistance")
#Population of agents over time
plot(time, pop.size[1:n],
     xlab = "generation", ylab = "Individuals")

#2. What happens if stress is decreasing over time
for(i in 1:n) {
  r = 0.2 - (0.01*time[i]/4)
  r.vec[i] = r
  pop.size[i+1] = pop.size[i]*exp(-r)
  p.resc[i] = get.resc.prop(gen  = gen,wt = wt,mut = mut,stress = r,
                            sel.coefficient = sel.coefficient,
                            mut.rate = mut.rate,
                            rep = rep)
}

par(mfrow = c(3,1))

#Plotting 3 figures
#Environmental stress over time
plot(time,r.vec,pch = 15,
     xlab = "generation", ylab = "env.stress (r)")
#Probability of resistance/rescue over time
plot(time,p.resc,pch=17,
     xlab = "generation", ylab = "Prob. resistance")
#Population of agents over time
plot(time, pop.size[1:n],
     xlab = "generation", ylab = "Individuals")

#3. What happens if stress is cycling, for example changing in a periodic way
for(i in 1:n) {
  r = 0.1 - (sin(0.1*time[i])/10)
  r.vec[i] = r
  pop.size[i+1] = pop.size[i]*exp(-r)
  p.resc[i] = get.resc.prop(gen  = gen,wt = wt,mut = mut,stress = r,
                            sel.coefficient = sel.coefficient,
                            mut.rate = mut.rate,
                            rep = rep)
}
par(mfrow = c(3,1))

#Plotting 3 figures
#Environmental stress over time
plot(time,r.vec,pch = 15,
     xlab = "generation", ylab = "env.stress (r)")
#Probability of resistance/rescue over time
plot(time,p.resc,pch=17,
     xlab = "generation", ylab = "Prob. resistance")
#Population of agents over time
plot(time, pop.size[1:n],
     xlab = "generation", ylab = "Individuals")
