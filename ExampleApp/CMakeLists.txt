PROJECT(ExampleApp)

FILE(GLOB headers ../include/MTL/*.h)
FILE(GLOB math_headers ../include/MTL/Math/*.h)
FILE(GLOB tools_headers ../include/MTL/Tools/*.h)

ADD_EXECUTABLE(Example Example.cpp ${headers} ${math_headers} ${tools_headers})

if (NOT MSVC)
  TARGET_LINK_LIBRARIES(Example stdc++fs)
endif()
