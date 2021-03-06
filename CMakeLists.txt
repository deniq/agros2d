# Set project's name
PROJECT( agros2d )
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# For Win64
if(${CMAKE_CL_64})
  set(WIN64 YES)
else(${CMAKE_CL_64})
  set(WIN64 NO)
endif(${CMAKE_CL_64})

# Debug yes / no.
SET(AGROS_DEBUG NO)
SET(WITH_QT5 NO)
SET(CMAKE_AGROS_DIRECTORY "${CMAKE_HOME_DIRECTORY}/")

set(WITH_OPENMP YES)
IF(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(WITH_OPENMP NO)
ENDIF()

# Allow to override the default values in CMake.vars:
INCLUDE(CMake.vars OPTIONAL)

# Handle Qt.
include(QtHandling.cmake)

# versions
SET(VERSION_MAJOR 3)
SET(VERSION_MINOR 2)
SET(VERSION_SUB 0)

# git revision
#execute_process(COMMAND git rev-list HEAD --count OUTPUT_VARIABLE NUMBER_OF_COMMITS)
#STRING(REGEX REPLACE "(\r?\n)+$" "" VERSION_GIT "${NUMBER_OF_COMMITS}")

# date
MACRO (TODAY RESULT)
    IF (WIN32)
        # EXECUTE_PROCESS(COMMAND "date" "/T" OUTPUT_VARIABLE ${RESULT})
        EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1" ${RESULT} ${${RESULT}})
    ELSEIF (UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1" ${RESULT} ${${RESULT}})
    ELSE (WIN32)
        MESSAGE(SEND_ERROR "date not implemented")
        SET(${RESULT} 000000)
    ENDIF (WIN32)
ENDMACRO (TODAY)

TODAY(RESULT)
string(SUBSTRING ${RESULT} 0 2 YEAR)
string(SUBSTRING ${RESULT} 2 2 MONTH)
string(SUBSTRING ${RESULT} 4 2 DAY)
SET(VERSION_YEAR 20${YEAR})
SET(VERSION_MONTH ${MONTH})
SET(VERSION_DAY ${DAY})



# Naming of projects built.
SET(AGROS_UTIL agros2d_util)
SET(PYTHONLAB_LIBRARY agros2d_pythonlab_library)
SET(PYTHONLAB_BINARY pythonlab_binary)
SET(AGROS_LIBRARY agros2d_library)
SET(AGROS_SOLVER agros2d_solver)
SET(AGROS_GENERATOR agros2d_generator)
SET(AGROS_BINARY agros2d_binary)
SET(AGROS_REMOTE agros2d_remote)
SET(PARALUTION_LIBRARY agros2d_3dparty_paralution)
SET(BSON_LIBRARY agros2d_3dparty_bson)
SET(CTEMPLATE_LIBRARY agros2d_3dparty_ctemplate)
SET(DXFLIB_LIBRARY agros2d_3dparty_dxflib)
SET(POLY2TRI_LIBRARY agros2d_3dparty_poly2tri)
SET(QCUSTOMPLOT_LIBRARY agros2d_3dparty_qcustomplot)
SET(QUAZIP_LIBRARY agros2d_3dparty_quazip)
SET(STB_TRUETYPE_LIBRARY agros2d_3dparty_stb_truetype)
SET(QTSINGLEAPPLICATION_LIBRARY agros2d_3dparty_qtsingleapplication)
SET(MATIO_LIBRARY agros2d_3dparty_matio)
SET(TRIANGLE_LIBRARY agros2d_3dparty_triangle)

SET(AGROS_SOLVER_EXTERNAL solver_external)

# Hermes and Hermes common
IF(MSVC)
  IF(WIN64)
    SET(HERMES_COMMON_LIBRARY hermes_common_64)
    SET(HERMES_LIBRARY hermes2d_64)
  ELSE(WIN64)
    SET(HERMES_COMMON_LIBRARY hermes_common)
    SET(HERMES_LIBRARY hermes2d)
  ENDIF(WIN64)
ELSE(MSVC)
  IF(AGROS_DEBUG)
    SET(HERMES_COMMON_LIBRARY hermes_common-debug)
    SET(HERMES_LIBRARY hermes2d-debug)
  ElSE(AGROS_DEBUG)
    SET(HERMES_COMMON_LIBRARY hermes_common)
    SET(HERMES_LIBRARY hermes2d)
  ENDIF(AGROS_DEBUG)
ENDIF()

# Output paths.
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_AGROS_DIRECTORY})
SET(LIBRARY_OUTPUT_PATH ${CMAKE_AGROS_DIRECTORY}/libs)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_AGROS_DIRECTORY})
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_AGROS_DIRECTORY}/libs)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_AGROS_DIRECTORY}/libs)

IF(WIN32)
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_AGROS_DIRECTORY})
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_AGROS_DIRECTORY})
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_AGROS_DIRECTORY}/libs)
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_AGROS_DIRECTORY}/libs)
  SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_AGROS_DIRECTORY}/libs)
  SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_AGROS_DIRECTORY}/libs)
ENDIF(WIN32)

# Set global compiler parameters.
IF (WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()
IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated")
  INCLUDE_DIRECTORIES(/usr/include/google)  
ENDIF()

IF(MSVC)
  INCLUDE_DIRECTORIES(c:/hpfem/hermes/dependencies/include)
  INCLUDE_DIRECTORIES(d:/hpfem/hermes/dependencies/include)
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /NODEFAULTLIB:libcmtd /NODEFAULTLIB:libcmt")
ENDIF(MSVC)

# Python
FIND_PACKAGE(PythonLibs 2.7 REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})

# Include current dir
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

# Turn on Unicode
ADD_DEFINITIONS(-DUNICODE -D_UNICODE)

# Enable debugging symbols for all files all the time:
IF(AGROS_DEBUG)
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Disable all warnings and turn on only important ones:
    # SET(CMAKE_CXX_FLAGS "-w ${CMAKE_CXX_FLAGS} -fPIC")
    # SET(CMAKE_CXX_FLAGS "-Wuninitialized -Wvla -Wsign-compare ${CMAKE_CXX_FLAGS}")
    SET(DEBUG_FLAGS   "-g")
    SET(CMAKE_BUILD_TYPE Debug)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEBUG_FLAGS}")
  ENDIF()
ElSE(AGROS_DEBUG)
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    SET(CMAKE_BUILD_TYPE Release)
    SET(RELEASE_FLAGS "-DNDEBUG -g -O3 -Ofast")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${RELEASE_FLAGS}")
  ENDIF()
ENDIF(AGROS_DEBUG)

# This overrides CXX flags for MSVC
IF(MSVC)
  # On MSVC, we will always use Debug, otherwise the plugins do not work.
  SET(CMAKE_BUILD_TYPE Debug)
  SET(MSVC_DEFINES "/DWIN32 /D_WINDOWS /Dpopen=_popen /Dpclose=_pclose /D__value=_value /Dfinite=_finite /Dhypot=_hypot /Disatty=_isatty /Dfileno=_fileno /D_CRT_SECURE_NO_WARNINGS /MP")
  SET(CMAKE_CXX_FLAGS_DEBUG_INIT          "/D_DEBUG /Od /Ob2 /MDd /Zi ${MSVC_DEFINES}")
  SET(CMAKE_CXX_FLAGS_RELEASE_INIT        "/DNDEBUG /O2 /Ob2 /MD ${MSVC_DEFINES}")
  SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SAFESEH:NO")
  SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /SAFESEH:NO")
  SET (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /SAFESEH:NO")
ENDIF(MSVC)

# Create file with preprocessor definitions exposing the build settings to the source code.
configure_file(
  util/config.h.cmake
  util/config.h)
configure_file(
  resources_source/resources.rc.cmake
  resources_source/resources.rc)

# Include OUR header files location
INCLUDE(IncludeSubdirs.cmake)

# Subdirectories.
ADD_SUBDIRECTORY(3rdparty)
# Also looks for Xerces, UMFPACK, etc.
ADD_SUBDIRECTORY(hermes)
# Include external directories.
INCLUDE_DIRECTORIES(${UMFPACK_INCLUDE_DIR})

message("CMake: Start Agros2D")
message("=====================")
ADD_SUBDIRECTORY(util)
ADD_SUBDIRECTORY(pythonlab-library)
ADD_SUBDIRECTORY(pythonlab-binary)
ADD_SUBDIRECTORY(agros2d-library)
ADD_SUBDIRECTORY(agros2d-binary)
ADD_SUBDIRECTORY(agros2d-solver)
ADD_SUBDIRECTORY(agros2d-generator)
ADD_SUBDIRECTORY(agros2d-remote)
ADD_SUBDIRECTORY(solver-external)
#ADD_SUBDIRECTORY(agros2d-forms)

message("\nAgros2D Configuration results")
message("---------------------")
message(" Version: ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_SUB}")

IF(AGROS_DEBUG)
  message(" Build version: Debug")
ELSE(AGROS_DEBUG)
    message(" Build version: Release")
ENDIF(AGROS_DEBUG)
message("---------Qt---------")
IF(WITH_QT5)
  message(" Build with Qt5")
ELSE(WITH_QT5)
  message(" Build with Qt4")
  message("  Qt use file: ${QT_USE_FILE}")
  message("  Qt definitions: ${QT_DEFINITIONS}")
  message("  Qt libraries: ${QT_LIBRARIES}")
ENDIF(WITH_QT5)
message("---------------------")
message("\n")

# pixmap
FILE(GLOB PIXMAP_FILES "resources/images/agros2d.xpm")
# desktop
FILE(GLOB DESKTOP_FILES "agros2d.desktop")
# resources
SET(RESOURCES_FILES "resources")

# share
INSTALL(FILES ${PIXMAP_FILES} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/pixmaps)
INSTALL(FILES ${DESKTOP_FILES} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/applications)
INSTALL(DIRECTORY ${RESOURCES_FILES} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/agros2d)

# other files
MACRO(ADD_INFO_FILES)
   FOREACH(_it ${ARGN})
     IF(NOT IS_DIRECTORY ${_it})
       GET_FILENAME_COMPONENT(name ${_it} NAME)
       IF(NOT ${_it} MATCHES "^/\\\\..*$;~$")
         SET_SOURCE_FILES_PROPERTIES(${_it} PROPERTIES HEADER_FILE_ONLY TRUE)
       ENDIF()
     ENDIF()
   ENDFOREACH()
ENDMACRO()

set(info resources/modules/acoustic.xml
         resources/modules/current.xml
         resources/modules/elasticity.xml
         resources/modules/electrostatic.xml
         resources/modules/flow.xml
         resources/modules/heat.xml
         resources/modules/magnetic.xml
         resources/modules/rf_te.xml
         resources/modules/rf_tm.xml
         resources/modules/richards.xml
         resources/couplings/current-heat.xml
         resources/couplings/heat-elasticity.xml
         resources/couplings/magnetic-heat.xml
         resources/python/functions_agros2d.py
         resources/python/functions_pythonlab.py
         resources/templates/generator/coupling_CMakeLists_txt.tpl
         resources/templates/generator/coupling_interface_cpp.tpl
         resources/templates/generator/coupling_interface_h.tpl
         resources/templates/generator/filter_cpp.tpl
         resources/templates/generator/filter_h.tpl
         resources/templates/generator/errorcalculator_cpp.tpl
         resources/templates/generator/errorcalculator_h.tpl
         resources/templates/generator/force_cpp.tpl
         resources/templates/generator/force_h.tpl
         resources/templates/generator/localvalue_h.tpl
         resources/templates/generator/localvalue_cpp.tpl
         resources/templates/generator/surfaceintegral_cpp.tpl
         resources/templates/generator/surfaceintegral_h.tpl
         resources/templates/generator/volumeintegral_cpp.tpl
         resources/templates/generator/volumeintegral_h.tpl
         resources/templates/generator/interface_h.tpl
         resources/templates/generator/module_CMakeLists_txt.tpl
         resources/templates/generator/weakform_cpp.tpl
         resources/templates/generator/weakform_h.tpl
         resources/templates/generator/coupling_weakform_cpp.tpl
         resources/templates/generator/coupling_weakform_h.tpl
         resources/templates/generator/plugins_CMakeLists_txt.tpl
         resources/templates/generator/interface_cpp.tpl
         resources/templates/generator/doc_index.tpl
         resources/templates/generator/equations.tpl
         resources/templates/generator/equations.tex
         resources/templates/panels/empty.tpl
         resources/templates/panels/example_form.tpl
         resources/templates/panels/example_problem.tpl
         resources/templates/panels/example_python.tpl
         resources/templates/panels/integrals.tpl
         resources/templates/panels/local_point_values.tpl
         resources/templates/panels/material.tpl
         resources/templates/panels/problem.tpl
         resources/templates/panels/welcome.tpl
         resources/templates/panels/style_common.css
         resources/templates/panels/style_results.css
         resources_source/python/agros2d.pyx
         resources_source/python/pyfield.pxi
         resources_source/python/pygeometry.pxi
         resources_source/python/pyparticletracing.pxi
         resources_source/python/pyproblem.pxi
         resources_source/python/pythonlab.pyx
         resources_source/python/pyview.pxi
)
add_info_files(${info})
add_custom_target(other SOURCES ${info})
