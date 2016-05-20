find_package(SmartFactory)
find_package(HistAsymmErrors)

if (SmartFactory_FOUND)
include_directories(
	${SMARTFACTORY_INCLUDE_DIR}
)

link_directories(
	${SMARTFACTORY_LIBRARY_DIR}
)

add_executable(example_write example_write.cpp)
target_link_libraries(example_write ${ROOT_LIBRARIES} ${SMARTFACTORY_LIBRARIES} HistAsymmErrors)

add_executable(example_read example_read.cpp)
target_link_libraries(example_read ${ROOT_LIBRARIES} ${SMARTFACTORY_LIBRARIES} HistAsymmErrors)

add_executable(example_multi_read example_multi_read.cpp)
target_link_libraries(example_multi_read ${ROOT_LIBRARIES} ${SMARTFACTORY_LIBRARIES} HistAsymmErrors)

add_executable(example_display example_display.cpp)
target_link_libraries(example_display ${ROOT_LIBRARIES} HistAsymmErrors)

endif(SmartFactory_FOUND)

add_executable(example_display2 example_display2.cpp)
target_link_libraries(example_display2 ${ROOT_LIBRARIES} HistAsymmErrors)
