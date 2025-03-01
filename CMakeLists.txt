cmake_minimum_required(VERSION 3.10)
project(c_ml C)

set(SRC src/main.c src/module1.c src/module2.c)
include_directories(include)

add_executable(my_project ${SRC})
