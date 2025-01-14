add_library(math
        goldbach/goldbach.cpp
        utils/prime_utils.cpp
)

target_include_directories(math
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(math
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/goldbach
            ${CMAKE_CURRENT_SOURCE_DIR}/utils
)

add_executable(goldbach_verifier main.cpp)
target_link_libraries(goldbach_verifier PRIVATE math)
