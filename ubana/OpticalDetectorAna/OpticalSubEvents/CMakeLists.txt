# find_ups_product is meant to be called from the top CMakeLists.txt file

#add_subdirectory( cfdiscriminator_algo )
#add_subdirectory( subevent_algo )

#SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -D__BUILD_ROOT_DICT__" )

art_make( LIBRARY_NAME OpticalSubEvents
	  SUBDIRS cfdiscriminator_algo subevent_algo	  
	  EXCLUDE SubEventBuilder_module.cc
	  LIB_LIBRARIES
	       ubobj_Optical
	       ${ART_FRAMEWORK_CORE} 
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas
               art_Persistency_Provenance canvas
               art_Utilities canvas
	       ${ROOT_BASIC_LIB_LIST}	       
	)
		

simple_plugin( SubEventBuilder "module"
	       OpticalSubEvents
	       ubobj_Optical
	       lardataobj_RecoBase
	       larcorealg_Geometry
	       larcore_Geometry_Geometry_service
	       lardata_Utilities
               ${ART_FRAMEWORK_CORE}
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas
               art_Persistency_Provenance canvas
               art_Utilities canvas
	       ${ROOT_BASIC_LIB_LIST}	       
               BASENAME_ONLY
             )


install_headers()
install_source()

install_fhicl()
FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )
