# CMakeLists are the devil and anybody who tells you they understand them is a liar
# This one is stolen from the AnalysisTree directory, which is a good bet for 
# a first try.

art_make( 
  MODULE_LIBRARIES 
    larcorealg_Geometry
    larcore_Geometry_Geometry_service
    larsim_Simulation
    nutools_ParticleNavigation
    lardataobj_Simulation
    lardata_Utilities
    larevt_Filters
    lardataobj_RawData
    lardataobj_RecoBase 
    lardata_RecoObjects 
    larreco_Deprecated
    larreco_RecoAlg
    larpandora_LArPandoraInterface
    SignalShapingServiceMicroBooNE_service
    nusimdata_SimulationBase
    ubraw_RawData
    ubobj_Trigger
    ubobj_Optical
    ${LARRECO_LIB}
    ${LARDATA_LIB}
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    art_Persistency_Common canvas
    art_Persistency_Provenance canvas
    art_Utilities canvas
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    cetlib 
    cetlib_except
    ${ROOT_GEOM}
    ${ROOT_XMLIO}
    ${ROOT_GDML}
    ${ROOT_BASIC_LIB_LIST}
  )

# install_headers()
install_fhicl()
install_source()
