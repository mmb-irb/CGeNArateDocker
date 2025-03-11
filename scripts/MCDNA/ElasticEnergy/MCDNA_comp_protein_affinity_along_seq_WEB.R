#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


# Get stiffness and ground state values from sequence #

#library("IRanges")


if (length(args)<5) {
  stop("Five arguments must be supplied: seq,pdbCode,position,helparmsPath,stiffPath,stiffTetraPath", call.=FALSE)
		     }

seq = as.character(args[1])
pdb = as.character(args[2])
pos = as.numeric(args[3])
helprms = as.character(args[4])
stif_file = as.character(args[5])
stif_tetra_file = as.character(args[6])
#in_folder = as.character(args[2])
#out_folder = as.character(args[4])

#inp_prot = read.table(sprintf("%s/mcdna_prot.config",in_folder))
#class(inp_prot[,1]) = "numeric"
#class(inp_prot[,2]) = "character"
#class(inp_prot[,3]) = "numeric"
#prot_files = as.character(inp_prot[1,2])
#helprms = Sys.getenv("MCDNA_PROT_HELPRMS")
prot_files = sprintf("%s/%s.helprms",helprms,pdb)
prot_conf = read.table(prot_files)

# To test you can use and out_folder as below
#seq_file="/orozco/homes/pluto/jwalther/Programs/Outputs_MCDNA/atomistic_protdna/seq_90bp.dat"
#in_folder="/orozco/homes/pluto/jwalther/Programs/Outputs_MCDNA/atomistic_protdna"
#out_folder="/orozco/homes/pluto/jwalther/Programs/Outputs_MCDNA/atomistic_protdna"
###



#seq = scan(sprintf("%s",seq_file), what="string", sep=NULL)
stif = read.table(file = stif_file, fill=TRUE)
stif_tetra = read.table(file = stif_tetra_file, fill=TRUE)


seq_indsteps = sapply(1:(nchar(seq)), function(i) substr(seq, start=i, stop=(i)))
seq_bpsteps = sapply(1:(nchar(seq)-1), function(i) substr(seq, start=i, stop=(i+1)))
seq_tetsteps = sapply(1:(nchar(seq)-3), function(i) substr(seq, start=i, stop=(i+3)))



get_stif = function(i)
{

if(i %in% c(1,(nchar(seq[1])-1))){
for(j in seq(from=1, to=length(stif[,1]), by=8)){		
	if (as.character(stif[j,1]) == seq_bpsteps[i]){ #print("di") 
							#print(seq_bpsteps[i])
							out = stif[(j+1):(j+7),]
							}
						}
	    } else{
for(j in seq(from=1, to=length(stif_tetra[,1]), by=8)){		
	if (as.character(stif_tetra[j,1]) == seq_tetsteps[i-1]){ #print("tet") 
								 #print(seq_tetsteps[i-1])
	   							 out = stif_tetra[(j+1):(j+7),]
						   		}
							}
	    }


return(out)
}

stif_all = lapply(1:(nchar(seq[1])-1), get_stif)
gs_all = do.call("rbind", lapply(stif_all, function(x)c(as.numeric(as.character(x[7,1])),unlist(x[7,2:6]))))




#Function
iterate_inp = function(inp){

gs_red = gs_all[inp:(inp+length(prot_conf[,1])-1),]
stif_red = stif_all[inp:(inp+length(prot_conf[,1])-1)]

diff = prot_conf - gs_red

calc_ener = function(i,diff){

# Function to obtain individual energy for a given bp-step in each snapshot of the trajectory

        vec1 = as.numeric(diff[i,])     # select i snapshot
        mat = as.matrix(stif_red[[i]][1:6,])
	class(mat)="numeric"
        ener_bp = vec1 %*% mat %*% vec1 # individual energy per bp-step (for first bp-step)
	return(ener_bp)
                         }

energy = unlist(lapply(1:length(stif_red),calc_ener, diff))

el_ener = energy/2

return(el_ener)
}

results = lapply(c(1:(length(gs_all[,1])-length(prot_conf[,1]))),iterate_inp)
en_penalty = unlist(lapply(results, function(x) sum(x)))
out_energy_penalty = cbind(1:length(en_penalty),en_penalty)
colnames(out_energy_penalty) = c("Number of base pair","Energy (kcal/mol)")
out_char = sprintf("The energy penalty at the initial position you indicated (%d) is %.2f. The minimum energy penalty of %.2f is found at the initial position %d.",pos,round(en_penalty[pos],2),round(min(en_penalty),2),which(en_penalty == min(en_penalty)))
# Output values
#write.csv(out_energy_penalty, file=sprintf("%s/protein_energy_penalty.csv",out_folder))
write.csv(out_energy_penalty)
#write(out_char, file=sprintf("%s/protein_energy_penalty_initial_position.dat",out_folder))

