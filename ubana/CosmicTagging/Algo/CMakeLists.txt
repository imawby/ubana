# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMinuit -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make( BASENAME_ONLY
    	  LIBRARY_NAME     ubana_CosmicTagging_CosmicTaggingAlgorithms
	  LIB_LIBRARIES
                	  larsim_PhotonPropagation
			  larsim_PhotonPropagation_PhotonVisibilityService_service
			  ubana_BasicTool_GeoAlgo
			  ubreco_SelectionTool_OpT0FinderBase
                          UBOpReadoutMap_service
  		          ${ART_FRAMEWORK_CORE}
   		          ${ART_FRAMEWORK_PRINCIPAL}
		          ${ART_FRAMEWORK_BASIC}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
   		          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
     		          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
     		          ${ART_PERSISTENCY_COMMON}
                          ${ART_PERSISTENCY_PROVENANCE}
                          ${ART_UTILITIES}
                          ${ROOT_BASIC_LIB_LIST}
                          ${FHICLCPP}
                          cetlib cetlib_except
	)

install_headers()
#install_fhicl()
install_source()
