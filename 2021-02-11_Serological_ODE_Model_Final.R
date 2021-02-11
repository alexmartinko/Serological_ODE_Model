#Script used to model split Luciferase serological assay system in Elledge and Zhou et al NBT 2021
#Created by Alex Martinko 2020/07/08
#Comments updated on 2021/02/11

#Installation of the required deSolve package
library(deSolve)

#Function that simulates the described thermodynamic scheme over time (t), given starting concentrations (state), and Kf/Kr values for each rate (rates)
ODE_sim <- function(t, state, rates) {
  with(as.list(c(state, rates)), {
    #Ordinary differential equations (ODEs) to represent the serological assay as described in Extended Data Fig. 3A
    dA = -k1f*C*A - k1f*D*A - k1f*E*A + k1r*D + k1r*G + k1r*H
    dB = -k1f*C*B - k1f*E*B - k1f*D*B + k1r*E + k1r*I + k1r*H
    dC = -k1f*C*A - k1f*C*B + k1r*D + k1r*E
    dD = -k1r*D - k1f*D*A - k1f*D*B + k1f*C*A + k1r*G + k1r*H
    dE = -k1r*E - k1f*E*A - k1f*E*B + k1f*C*B + k1r*H + k1r*I
    dG = -k1r*G + k1f*D*A
    dH = -k1r*H - k1r*H + k1f*D*B  + k1f*E*A
    dI = -k1r*I + k1f*E*B
    list(c(dA, dB, dC,dD,dE,dG,dH,dI))
  })
}

########## Code below was used to generate Extended Data Fig. 3B
### A loop to iterate simulations over a range of Antibody concentrations (C)

# Input variables for the function below. Can be customized.
C_vals = c(30,10,3.33,1.111111,0.37037,0.123457,0.041152,0.013717,0.004572,0.001524,0.000508) #Range of Ab concentrations (C) tested in Extended Data Fig. 3B. These numbers were inspired by experimental data. Units are arbitrary, but were assumed as nM for this study. 
eq_matrix <- matrix( data=NA,nrow=length(C_vals), ncol=8) #Matrix is built to the dimensions of C_vals and the number species in the model.
rownames(eq_matrix) <- c(30,10,3.33,1.111111,0.37037,0.123457,0.041152,0.013717,0.004572,0.001524,0.000508) #Row names are filled with concentration (C) values from above. 
colnames(eq_matrix) <- c("A","B","C","D","E","G","H","I") #Col names are filled with the names of different species in ODE model (F is excluded as a letter choice b/c F=False in R)
times <- seq(0, 1000, by = .1) #This is the time range for the ODE calculator. This was determined to be enough time to reach equilibrium for the scenerios computed in the paper.
rates = c(k1f=1,k1r=0.001) #This is where the affinity of the Ab is defined. For the simulations in Extended Data Fig. 3B, k1r ranged from 0.001 to 100. Units for (k1r/k1f) are arbitrary, but were assumed as nM for this study

#Function to compute and output concentrations at equilibrium across a range of concentrations of C.
conc_matrix = function(x){ #x in the case of Extended Data Fig. 3B was the variable "C_vals" as defined above.
  counter <- 0 #Counter is used to calulcate percent completion.
  for (i in x){ #Loops through range of concentrations defined in x
    state <- c(A=1,B=1,C=i, D = 0, E=0, G=0, H=0, I=0) #Use this line to define/edit T=0 concentrations.
    out <- ode(y = state, times = times, func = ODE_sim, parms = rates) #Calls deSolve and simulates equilibration.
    eq_matrix[paste(i, sep = ''),]  <<-  #Finds the appropriate position in eq_matrix to populate.
      out[length(out[1,]),2:9]           #Appends the final row of "out" (the final concs at equilibration) to eq_matrix.
    counter <- counter + 1
    print(paste(round(counter/(length(x))*100, 2), "% completed", sep = '')) #Prints percent completion.
  } 
  write.csv(eq_matrix, file = "YYYY-MM-DD_Customize_your_file_name_eq_matrix_1pM_k1r.csv") #CSV file is generated with filled in matrix.
}

# To generate Extended Data Fig. 3B, the above function was run (i.e. conc_matrix(C_vals)) for k1r values ranging from 0.001 to 100. The values computed in column "I" in each case were plotted against the values in C_vals.



########## Code below was used to generate Extended Data Fig. 3C
### A loop to iterate simulations over a range of Antibody affinities


# Input variables for the function below. Can be customized.
K_vals = c(1:10 %o% 10^(-2:2)) #Range of values to scan across different Ab affinities (k1r/k1f). Units for (k1r/k1f) are arbitrary, but were assumed as nM for this study.
K_vals <- K_vals[-seq(11,51, 10)] #Removes redundant affinity vals. This is specific to the current range of vals.
kd_matrix <- matrix( data=NA,nrow=length(K_vals), ncol=8) #Matrix is built to the dimensions of K_vals and the number species in the model.
rownames(kd_matrix) <- c(K_vals) #Row names are filled with Kd inputs for K_vals.
colnames(kd_matrix) <- c("A","B","C","D","E","G","H","I") #Col names are filled with the names of different species in ODE model (F is excluded as a letter choice b/c F=False in R).
times <- seq(0, 1000, by = .1) #This is the time range for the ODE calculator. This was determined to be enough time to reach equilibrium for the scenerios computed in the paper.
state <- c(A=1,B=1,C=0.1, D = 0, E=0, G=0, H=0, I=0) #Use this line to edit T=0 concentrations. For Extended Data Fig. 3C, A/B were varried from 1 to 100, C was set to 0.1, and all others were set to 0. Units are arbitrary, but were assumed as nM for this study..

#Function to compute and output concentrations at equilibrium across a range of Antibody affinities (k1r)
kd_func = function(x){ #x in the case of Extended Data Fig. 3C was the variable "K_vals" as defined above.
  counter <- 0 #Counter is used to calulcate percent completion.
  for (i in x){ #Loops through range of k1r values defined in x
    rates = c(k1f=1,k1r=i)
    out <- ode(y = state, times = times, func = ODE_sim, parms = rates) #Calls deSolve and simulates equilibration.
    kd_matrix[paste(i, sep = ''),]  <<-   #Finds the appropriate position in kd_matrix to populate.
      out[length(out[1,]),2:9]           #Appends the final row of "out" (the final concs at equilibration) to kd_matrix,
    counter <- counter + 1
    print(paste(round(counter/(length(x))*100, 2), "% completed", sep = '')) #Prints percent completion.
  } 
  write.csv(kd_matrix, file = "YYYY-MM-DD_Customize_your_file_name_kq_matrix_1nM_AandB.csv") #CSV file is generated with filled in matrix.
}

# To generate Extended Data Fig. 3C, the above function was run (i.e. kd_func(K_vals)) for A and B concentration values ranging from 1 to 100. The values computed in column "I" for each case were plotted against the values in K_vals.


