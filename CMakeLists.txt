cmake_minimum_required( VERSION 2.8 )

project( VTK_STRAIN_GLYPH )

find_package( VTK REQUIRED )
include( ${VTK_USE_FILE} )

enable_testing()
add_subdirectory( test )
