cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_C_FLAGS_DEBUG "-g -O0")
#set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -g -O1 -fno-omit-frame-pointer -fno-optimize-sibling-calls")



set(PROJECT_AND_EXE_NAME "AstroGen2")
project(${PROJECT_AND_EXE_NAME})


# Setzen Sie den Ausgabeordner für Bibliotheken
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Setzen Sie den Zwischenspeicherordner für Objektdateien
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj)

# Alle .cpp Dateien rekursiv finden
file(GLOB_RECURSE USER_SOURCE_FILES
    src/*.cpp
)

# Funktion, um rekursiv alle Verzeichnisse zu finden
function(add_subdirectories_recursive dir)
    file(GLOB children RELATIVE ${dir} ${dir}/*)
    foreach(child ${children})
        if(IS_DIRECTORY ${dir}/${child})
            include_directories(${dir}/${child})
        endif()
    endforeach()
endfunction()


# Hauptverzeichnis und alle Unterverzeichnisse hinzufügen
include_directories(${PROJECT_SOURCE_DIR}/src)
add_subdirectories_recursive(${PROJECT_SOURCE_DIR}/src)

find_package(OpenMP REQUIRED)

add_executable(${PROJECT_AND_EXE_NAME} ${USER_SOURCE_FILES})

if(OpenMP_FOUND)
    target_link_libraries(AstroGen2 PUBLIC OpenMP::OpenMP_CXX)
endif()

# Compiler-Optionen und Link-Bibliotheken
target_compile_options(${PROJECT_AND_EXE_NAME} PRIVATE -Wall)

#target_link_libraries(${PROJECT_AND_EXE_NAME} ws2_32 gdi32)

# Define custom target to clean build directory
add_custom_target(clean-build
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
)

# Define custom target to build using an external terminal
add_custom_target(build-in-external-terminal
    COMMAND start cmd /c "${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --config Debug"
    DEPENDS ${PROJECT_AND_EXE_NAME}
)
