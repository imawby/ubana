
art_make( BASENAME_ONLY
	  LIBRARY_NAME  ubana_OpticalDetectorAna
	  LIB_LIBRARIES
			larana_OpticalDetector
			larana_OpticalDetector_OpHitFinder
                        ${ROOT_BASIC_LIB_LIST}

	  MODULE_LIBRARIES larcorealg_Geometry
	  	        larcore_Geometry_Geometry_service
		        lardata_Utilities
	                ubcore_Geometry
			UBOpReadoutMap_service
                        lardataobj_OpticalDetectorData
			larana_OpticalDetector
                        lardataobj_RawData
			lardataobj_RecoBase
			larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			nusimdata_SimulationBase
			larsim_PhotonPropagation
			ubana_OpticalDetectorAna
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
			nusimdata_SimulationBase
        )							   

install_headers()
install_source()
install_fhicl()
add_subdirectory(job)
add_subdirectory(OpticalSubEvents)
add_subdirectory(DLPMTPreCuts)
