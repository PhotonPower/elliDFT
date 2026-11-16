find_package(Catch2 3 REQUIRED)
include(Catch)

# Unit tests: single functions, no reference data needed.
add_executable(unit_tests
  unit/test_core.cpp
  unit/test_density_fitting.cpp
  unit/test_xyz.cpp
  unit/test_basis_library.cpp
  unit/test_basis_set.cpp
  unit/test_boys.cpp
  unit/test_hermite.cpp
  unit/test_one_electron.cpp
  unit/test_dispersion_fit.cpp
  unit/test_ellipsometry_fit.cpp
  unit/test_ellipsometry_materials2.cpp
  unit/test_ellipsometry_model.cpp
  unit/test_ellipsometry_stack2.cpp
  unit/test_least_squares.cpp
  unit/test_optics.cpp
  unit/test_range_separation.cpp
  unit/test_nonlinear.cpp
  unit/test_realtime.cpp
  unit/test_response_operator.cpp
  unit/test_response_properties.cpp
  unit/test_two_electron.cpp
  unit/test_xc.cpp
  unit/test_xc_kernel.cpp
  unit/test_jk.cpp
  unit/test_kohn_sham.cpp
  unit/test_basis_on_grid.cpp
  unit/test_lebedev.cpp
  unit/test_molecular_grid.cpp
  unit/test_rhf.cpp
  unit/test_scf_options.cpp
  unit/test_normalization.cpp
  unit/test_solid_harmonics.cpp
)
target_link_libraries(unit_tests PRIVATE
  ellidft Catch2::Catch2WithMain nlohmann_json::nlohmann_json ellidft_warnings)
target_compile_definitions(unit_tests PRIVATE ELLIDFT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_include_directories(unit_tests PRIVATE support)
catch_discover_tests(unit_tests TEST_PREFIX "unit:" ADD_TAGS_AS_LABELS)

# Regression tests: comparison with the stored reference data
# (tests/reference/data, produced by the PySCF scripts in tests/reference).
add_executable(regression_tests
  regression/test_m0_basis.cpp
  regression/test_m1_reference_data.cpp
  regression/test_m1_one_electron.cpp
  regression/test_m1_two_electron.cpp
  regression/test_m1_hartree_fock.cpp
  regression/test_m2_grid.cpp
  regression/test_m2_basis_values.cpp
  regression/test_m2_kohn_sham.cpp
  regression/test_m3_density_fitting.cpp
  regression/test_m4_hybrids.cpp
  regression/test_m5_response.cpp
)
target_link_libraries(regression_tests PRIVATE
  ellidft Catch2::Catch2WithMain nlohmann_json::nlohmann_json ellidft_warnings)
target_compile_definitions(regression_tests PRIVATE ELLIDFT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_include_directories(regression_tests PRIVATE support)
# Tests tagged [slow] get the ctest label "slow"; scripts/build_and_test.sh
# leaves them out unless called with --slow.
catch_discover_tests(regression_tests TEST_PREFIX "regression:" ADD_TAGS_AS_LABELS
  PROPERTIES TIMEOUT 7200)

# Program ellidft-ellipsometry: the example and a test input must give their
# known answers (data computed with the program itself from known models).
add_test(NAME app:ellipsometry-example
  COMMAND ellidft-ellipsometry fit ${CMAKE_SOURCE_DIR}/examples/ellipsometry/fit.json)
set_tests_properties(app:ellipsometry-example PROPERTIES
  PASS_REGULAR_EXPRESSION "film\\.thickness_nm +13(1\\.99|2\\.00)" LABELS "ellipsometry")
add_test(NAME app:ellipsometry-density
  COMMAND ellidft-ellipsometry fit ${CMAKE_SOURCE_DIR}/tests/data/ellipsometry/density_fit.json)
set_tests_properties(app:ellipsometry-density PROPERTIES
  PASS_REGULAR_EXPRESSION "film\\.density_g_cm3 +1\\.047 " LABELS "ellipsometry")

# Program ellidft-realtime: with --compare it reports its deviation from the
# linear response equations; it must be below 1e-4.
add_test(NAME app:realtime-h2
  COMMAND ellidft-realtime ${CMAKE_SOURCE_DIR}/data/geometries/h2.xyz cc-pVDZ --functional HF
          --gamma 0.15 --time-step 0.2 --time 90 --directions z --frequencies 0,0.1,0.3 --compare)
set_tests_properties(app:realtime-h2 PROPERTIES
  PASS_REGULAR_EXPRESSION "agreement with linear response: [0-9.]+e-(0[5-9]|[1-9][0-9])" LABELS "realtime")

# Program ellidft-nonlinear on its smallest example. H2 with s functions
# responds only along the bond: the orientational average is one fifth of the
# value along the bond, which is the third entry of the list by direction.
add_test(NAME app:nonlinear-h2-clean
  COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_SOURCE_DIR}/examples/nonlinear/h2_kerr.json.partial.json)
add_test(NAME app:nonlinear-h2
  COMMAND ellidft-nonlinear ${CMAKE_SOURCE_DIR}/examples/nonlinear/h2_kerr.json)
# second run: the finished wavelength comes from the file of the first
add_test(NAME app:nonlinear-h2-continued
  COMMAND ellidft-nonlinear ${CMAKE_SOURCE_DIR}/examples/nonlinear/h2_kerr.json)
set_tests_properties(app:nonlinear-h2-clean PROPERTIES FIXTURES_SETUP nonlinear_clean LABELS "nonlinear")
set_tests_properties(app:nonlinear-h2 PROPERTIES FIXTURES_REQUIRED nonlinear_clean FIXTURES_SETUP nonlinear_first
  PASS_REGULAR_EXPRESSION "Kerr effect at 780 nm: gamma[^\n]*residual [0-9.e-]+\n" LABELS "nonlinear")
set_tests_properties(app:nonlinear-h2-continued PROPERTIES FIXTURES_REQUIRED nonlinear_first
  PASS_REGULAR_EXPRESSION "Kerr effect at 780 nm: gamma[^\n]*from the earlier run" LABELS "nonlinear")

# Main program: control file in, results and protocol out.
add_test(NAME app:ellidft-control-file
  COMMAND ellidft-run ${CMAKE_SOURCE_DIR}/tests/data/control/h2.json)
set_tests_properties(app:ellidft-control-file PROPERTIES
  PASS_REGULAR_EXPRESSION "total energy -1\\.11671432.*results written to" LABELS "control")
add_test(NAME app:ellidft-rejects-unknown-keys
  COMMAND ellidft-run ${CMAKE_SOURCE_DIR}/tests/data/control/unknown_key.json)
set_tests_properties(app:ellidft-rejects-unknown-keys PROPERTIES
  PASS_REGULAR_EXPRESSION "scf: unknown key \"tolerance\"" LABELS "control")

# Checkpoint: the first run starts from nothing, the second continues.
add_test(NAME app:ellidft-checkpoint-clean
  COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_SOURCE_DIR}/tests/data/control/output/restart_test.chk)
add_test(NAME app:ellidft-checkpoint-first-run
  COMMAND ellidft-run ${CMAKE_SOURCE_DIR}/tests/data/control/restart.json)
add_test(NAME app:ellidft-checkpoint-second-run
  COMMAND ellidft-run ${CMAKE_SOURCE_DIR}/tests/data/control/restart.json)
set_tests_properties(app:ellidft-checkpoint-clean PROPERTIES FIXTURES_SETUP checkpoint_clean LABELS "control")
set_tests_properties(app:ellidft-checkpoint-first-run PROPERTIES
  FIXTURES_REQUIRED checkpoint_clean FIXTURES_SETUP checkpoint_written LABELS "control"
  PASS_REGULAR_EXPRESSION "restart    not used: no checkpoint file.*results written to")
set_tests_properties(app:ellidft-checkpoint-second-run PROPERTIES
  FIXTURES_REQUIRED checkpoint_written LABELS "control"
  PASS_REGULAR_EXPRESSION "density of a converged SCF, 3 finished response.*converged after 2 iterations.*from the checkpoint.*from the checkpoint.*from the checkpoint")

# Program ellidft-kerr-film on constructed inputs: gamma grows by exactly 1500 a.u. per unit at 780 nm.
add_test(NAME app:kerr-film
  COMMAND ellidft-kerr-film --density 1.18 --molar-mass 100.12 --refractive-index 1.49
          ${CMAKE_SOURCE_DIR}/tests/data/kerr_film/units_1.json ${CMAKE_SOURCE_DIR}/tests/data/kerr_film/units_2.json
          ${CMAKE_SOURCE_DIR}/tests/data/kerr_film/units_3.json)
set_tests_properties(app:kerr-film PROPERTIES
  PASS_REGULAR_EXPRESSION "at 780 nm: gamma per unit 1500 a.u." LABELS "nonlinear")

# Uniaxial film with a rough surface (effective medium) on a glass plate with back face, named materials.
add_test(NAME app:ellipsometry-glass
  COMMAND ellidft-ellipsometry fit ${CMAKE_SOURCE_DIR}/tests/data/ellipsometry/glass_fit.json)
set_tests_properties(app:ellipsometry-glass PROPERTIES
  PASS_REGULAR_EXPRESSION "film\\.thickness_nm +(240 |239\\.9999|240\\.0000)" LABELS "ellipsometry")
