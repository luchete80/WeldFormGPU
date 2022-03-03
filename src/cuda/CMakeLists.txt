add_library (cuda 
							#Functions.cu
							#Interaction.cu
							#NbSearch.cu
							#Particle_d.cu
							#Interaction.cu
							Domain_d.cu
							tensor.cu)
							
#set_property(TARGET cuda PROPERTY CUDA_SEPARABLE_COMPILATION ON)
