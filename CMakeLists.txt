cmake_minimum_required(VERSION 2.8)
project(ample NONE)

install(
    DIRECTORY
        ample
    DESTINATION
        lib/py2
)

install(
    DIRECTORY
        docs examples i2 include scripts testfiles
    DESTINATION
        share/ample
)

if (WIN32)
  install(
      PROGRAMS 
          bin/ample.bat 
      DESTINATION 
          bin
  )
else()
  install(
      PROGRAMS
          bin/ample
      DESTINATION
          bin
)
endif()
