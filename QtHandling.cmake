# Path to Qt5.
# Both WITH_QT5 and CMAKE_PREFIX_PATH can be changed in CMake.vars.
# This is just to handle the case WITH_QT5 is and CMAKE_PREFIX_PATH is not.
IF(WITH_QT5)
  IF(WIN32)
    IF(NOT(CMAKE_PREFIX_PATH))
      SET(CMAKE_PREFIX_PATH "C:\\Qt\\Qt5.0.2\\5.0.2\\msvc2010_opengl\\")
    ENDIF()
  ENDIF()
ENDIF()

# Add QT.
INCLUDE(${CMAKE_AGROS_DIRECTORY}/CMakeQt.cmake)
IF(DEFINED QT_USE_FILE)
  INCLUDE(${QT_USE_FILE})
ENDIF()
ADD_DEFINITIONS(${QT_DEFINITIONS})

# Link OpenGL and Zlib when needed.
IF(WITH_QT5)
  IF(WIN32)
    FIND_PACKAGE(OpenGL REQUIRED)
    FIND_PACKAGE(ZLIB REQUIRED)
  ENDIF()
ENDIF(WITH_QT5)