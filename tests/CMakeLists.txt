include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG main
)
FetchContent_MakeAvailable(googletest)

add_executable(goldbach_test
    goldbach/goldbach_test.cpp
)

target_link_libraries(goldbach_test
    PRIVATE
        math
        gtest
        gtest_main
)

add_test(NAME goldbach_test COMMAND goldbach_test)
