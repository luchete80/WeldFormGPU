add_library (cuda 
							#Functions.cu
							#Interaction.cu
							#NbSearch.cu
							#Particle_d.cu
							#Interaction.cu
							#Interaction_global.cu
							Domain_d.cu
              Geom.cu
							#tensor.cu)
							
#set_property(TARGET cuda PROPERTY CUDA_SEPARABLE_COMPILATION ON)
