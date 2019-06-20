project(cuda_library LANGUAGES CUDA)

add_library(mylib /home/fabio/uavlab/dune/src/MultiAgent/UAVagent/kernel.cu)

target_link_libraries(dune-core mylib)

