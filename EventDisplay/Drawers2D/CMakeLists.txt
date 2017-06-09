# basic source code CMakeLists.txt
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )

cet_enable_asserts()

art_make( BASENAME_ONLY
          TOOL_LIBRARIES  larcore_Geometry
                          lardataobj_RecoBase
                          larcore_Geometry_Geometry_service
                          lardata_Utilities
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          art_Persistency_Common
                          art_Persistency_Provenance
                          ${ART_FRAMEWORK_SERVICES_BASIC}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                          ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                          art_Utilities
                          canvas
                          ${MF_MESSAGELOGGER}
                          ${MF_UTILITIES}
                          ${FHICLCPP}
                          ${CETLIB}
                          ${Boost_FILESYSTEM_LIBRARY}
                          ${Boost_SYSTEM_LIBRARY}
                          ${CLHEP}
                          ${ROOT_EVE_LIB_LIST}
                          ${ROOT_EVE}
                          ${ROOT_EG}
                          ${ROOT_RGL}
                          ${ROOT_X3d}
                          ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_source()
install_fhicl()
