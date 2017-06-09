# Make sure we get the drawing tools
add_subdirectory(Drawers2D)

# basic source code CMakeLists.txt
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )

FILE( GLOB src_files *.cxx )

SET(CET_DICT_FUNCTIONS TRUE)


art_make( BASENAME_ONLY
          LIBRARY_NAME  EventDisplayLibs
          LIB_LIBRARIES lareventdisplay_EventDisplay_ColorDrawingOptions_service
	                    lareventdisplay_EventDisplay_InfoTransfer_service
                        lardataobj_Simulation
                        lardataobj_RawData
                        lardataobj_RecoBase
                        lardata_RecoObjects
                        lardata_Utilities
                        lardataobj_AnalysisBase
                        larevt_CalibrationDBI_IOVData
                        larevt_CalibrationDBI_Providers
                        ${SIMULATIONBASE}
                        ${EVENTDISPLAYBASE}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_IO}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common
                        art_Persistency_Provenance
                        art_Utilities
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_EVE_LIB_LIST}
                        ${ROOT_EVE}
                        ${ROOT_EG}
                        ${ROOT_RGL}
                        ${ROOT_X3d}
                        ${ROOT_BASIC_LIB_LIST}
      SERVICE_LIBRARIES EventDisplayLibs
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO}
                        ${ART_FRAMEWORK_IO_Root_RootInput_source}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common
                        art_Persistency_Provenance
                        art_Persistency_Common
                        art_Persistency_Provenance
                        art_Utilities
                        art_Framework_IO_Root_RootInput_source
                        art_Framework_IO_Root
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_EVE}
                        ${ROOT_BASIC_LIB_LIST}
	   MODULE_LIBRARIES EventDisplayLibs
	                    lareventdisplay_EventDisplay_InfoTransfer_service
                        lardataobj_Simulation
                        lardataobj_RawData
                        lardataobj_RecoBase
                        lardata_RecoObjects
                        lardata_Utilities
                        lardataobj_AnalysisBase
	                    ${SIMULATIONBASE}
                        ${EVENTDISPLAYBASE}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO}
                        ${ART_FRAMEWORK_IO_Root_RootInput_source}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common
                        art_Persistency_Provenance
                        art_Utilities
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_GUI}
                        ${ROOT_GEOM}
                        ${ROOT_EVE}
                        ${ROOT_EG}
                        ${ROOT_RGL}
                        ${ROOT_BASIC_LIB_LIST}
                )

install_headers()
install_source()
install_fhicl()
