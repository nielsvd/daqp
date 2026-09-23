set(CMAKE_C_STANDARD_LIBRARIES "${CMAKE_C_STANDARD_LIBRARIES}") # math for C api
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -Wextra")

add_executable(bug1 bug1.c)
target_link_libraries(bug1 daqpstat)

