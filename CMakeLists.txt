# Specifica la versione minima di CMake richiesta per questo progetto. Se l’utente ha una versione più vecchia, CMake interrompe la configurazione con un errore
cmake_minimum_required(VERSION 3.9)

# Definisce il nome del progetto (Project), i linguaggi usati (C++ e C) e la versione del progetto (1.0.0).
project(Project LANGUAGES CXX C VERSION 1.0.0)


################################################################################
# Set C++ Standard: imposta lo standard C++ da usare a C++20, e forza l'uso di quello standard (non versioni precedenti)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_minimum_required(VERSION 3.9)

project(Project LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Quando si compila in modalità Debug, abilita una serie di warning e controlli molto severi:
# -Wall: attiva molti warning utili.
# -Wextra: attiva ancora più warning.
# -pedantic-errors: segnala errori se il codice non è conforme allo standard.

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

#Su Windows, aggiunge un flag (-Wa,-mbig-obj) necessario quando si compila un file con molti simboli (come quelli generati da librerie template pesanti tipo Eigen).
# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

# Add dependencies
################################################################################
## Eigen3: Cerca la libreria Eigen3 installata nel sistema e la rende obbligatoria (REQUIRED). Aggiunge la libreria a quelle da linkare.
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND Project_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Threads: Cerca il supporto multi-threading (standard std::thread) e lo aggiunge alle librerie da linkare.
find_package(Threads REQUIRED)
list(APPEND Project_LINKED_LIBRARIES PRIVATE Threads::Threads)

## GTest; Trova la libreria Google Test e aggiunge tutti i componenti utili (gtest, gtest_main, gmock, gmock_main) al progetto per eseguire i test.
find_package(GTest REQUIRED)

list(APPEND Project_LINKED_LIBRARIES PRIVATE GTest::gtest)
list(APPEND Project_LINKED_LIBRARIES PRIVATE GTest::gtest_main)
list(APPEND Project_LINKED_LIBRARIES PRIVATE GTest::gmock)
list(APPEND Project_LINKED_LIBRARIES PRIVATE GTest::gmock_main)

# Insert Sources: Include la directory src, che presumibilmente ha un proprio CMakeLists.txt. Appende le variabili dichiarate lì (src_sources, ecc.) alle variabili principali del progetto. (successivamente anche per test)
################################################################################
add_subdirectory(src)

list(APPEND Project_SOURCES ${src_sources})
list(APPEND Project_HEADERS ${src_headers})
list(APPEND Project_INCLUDE ${src_includes})

add_subdirectory(src_test)

list(APPEND Project_TEST_SOURCES ${src_test_sources})
list(APPEND Project_TEST_HEADERS ${src_test_headers})
list(APPEND Project_TEST_INCLUDE ${src_test_includes})

# Altra cartella del progetto, probabilmente relativa all’esportazione di dati per Paraview. Appende i file relativi alla lista delle sorgenti principali.
add_subdirectory(ExportParaview)

list(APPEND Project_SOURCES ${src_paraview_sources})
list(APPEND Project_HEADERS ${src_paraview_headers})
list(APPEND Project_INCLUDE ${src_paraview_includes})


# Create executable
################################################################################
#Eseguibile principale; Crea l’eseguibile principale (chiamato Project) con main.cpp e tutti i file sorgenti e header.
add_executable(${PROJECT_NAME}
      main.cpp
      ${Project_SOURCES}
      ${Project_HEADERS})

# Configura: Le librerie da linkare, Le cartelle d’inclusione, L’opzione -fPIC, necessaria per creare codice riutilizzabile in librerie dinamiche (specialmente su Linux).
target_link_libraries(${PROJECT_NAME} ${Project_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${Project_INCLUDE})
target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

# Crea l’eseguibile di test (es. Project_TEST) con i file main_test.cpp e tutte le sorgenti (inclusi i test).
add_executable(${PROJECT_NAME}_TEST
      main_test.cpp
      ${Project_SOURCES}
      ${Project_HEADERS}
      ${src_test_sources}
      ${src_test_headers})

#Come sopra: collega librerie, include cartelle, abilita compilazione con -fPIC
target_link_libraries(${PROJECT_NAME}_TEST ${Project_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME}_TEST PRIVATE ${Project_INCLUDE} 
${src_test_includes})
target_compile_options(${PROJECT_NAME}_TEST PUBLIC -fPIC)

############################################################################
#riassunto: Questo CMakeLists.txt:

# - Richiede CMake >= 3.9.
# - Usa C++20.
# - Abilita warning e flag specifici in Debug e su Windows.
# - Include librerie come Eigen3, Threads, e GTest.
# - Organizza il progetto in tre sottodirectory: src, src_test, ExportParaview.
# - Crea due eseguibili: uno per il programma principale, uno per i test.
