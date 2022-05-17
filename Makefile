##
## FEM benchmarks
##
bindir := bin
benchmarks = \
	poisson2d_cg \
	poisson3d \
	assembly \
	reorder_mesh \
	stream \
	gather \
	pointer_chase
libs = util linear_algebra linear_solver mesh fem variational_form solvers benchmarks_common
unittest = $(bindir)/unittest

clean-benchmarks = $(benchmarks:%=%-clean)
clean-libs = $(libs:%=%-clean)
clean-unittest = $(unittest:%=%-clean)

all: $(libs) $(benchmarks)

check: run-unit-tests

clean: $(clean-benchmarks) $(clean-libs) $(clean-unittest)
	rm -rf docs/html

html:
	doxygen docs/doxygen.config

.PHONY: all check clean html
.PHONY: $(libs) run-unit-tests

##
## Configuration
##

INCLUDES = -iquote src

CFLAGS += -g -Wall

ifdef ENABLE_ASAN
CFLAGS += -fsanitize=address
CXXFLAGS += -fsanitize=address
endif

ifdef ENABLE_UBSAN
CFLAGS += -fsanitize=undefined
CXXFLAGS += -fsanitize=undefined
endif

ifndef NO_OPENMP
CFLAGS += -fopenmp
CXXFLAGS += -fopenmp
endif

ifndef NO_MPI
MPI_LIBS += -L"${MPI_HOME}/lib" -lmpi
CFLAGS += -DHAVE_MPI
CXXFLAGS += -DHAVE_MPI
endif

ifndef NO_METIS
METIS_INCLUDES = -isystem "${METIS_ROOT}/include"
METIS_LIBS = -L"${METIS_ROOT}/lib" -lmetis
CFLAGS += -DHAVE_METIS
CXXFLAGS += -DHAVE_METIS
endif

ifndef NO_PETSC
PETSC_INCLUDES = -isystem "${PETSC_DIR}/include"
PETSC_LIBS = -L"${PETSC_DIR}/lib" -lpetsc
CFLAGS += -DHAVE_PETSC
CXXFLAGS += -DHAVE_PETSC
endif

ifndef NO_CUDA
ifndef NVCC
NVCC = nvcc
endif
CUDA_CFLAGS = $(foreach x,$(CFLAGS),--compiler-options $(x))
CUDA_LIBS = -lcudart -lstdc++
endif

ifndef NO_LIBPFM
CFLAGS += -DHAVE_LIBPFM
CXXFLAGS += -DHAVE_LIBPFM
LDFLAGS += -lpfm
endif

ifeq ($(NO_GTEST),)
ifdef GTEST_ROOT
GTEST_INCLUDES = -isystem "$(GTEST_ROOT)/include"
GTEST_LIBS = -Wl,-rpath,"$(GTEST_ROOT)/lib" -L"$(GTEST_ROOT)/lib"
endif
GTEST_LIBS += -lgtest -lgtest_main
endif


##
## Various utilities
##
util_a = src/util/util.a
util_c_sources = \
	src/util/array.c \
	src/util/combinations.c \
	src/util/cuthill_mckee.c \
	src/util/matrix_market.c \
	src/util/parse.c \
	src/util/partition.c \
	src/util/perf_events.c \
	src/util/perf_sessions.c \
	src/util/permutation.c \
	src/util/prefix_sum.c \
	src/util/print.c \
	src/util/radix_sort.c \
	src/util/sample_statistics.c \
	src/util/selections.c \
	src/util/sort.c \
	src/util/stream_compaction.c
util_c_headers = \
	src/util/array.h \
	src/util/binary.h \
	src/util/combinations.h \
	src/util/cuthill_mckee.h \
	src/util/matrix_market.h \
	src/util/parse.h \
	src/util/partition.h \
	src/util/perf_events.h \
	src/util/perf_sessions.h \
	src/util/permutation.h \
	src/util/power_of_two_radix.h \
	src/util/prefix_sum.h \
	src/util/print.h \
	src/util/sample_statistics.h \
	src/util/selections.h \
	src/util/sort.h \
	src/util/stream_compaction.h
util_c_objects := $(foreach x,$(util_c_sources),$(x:.c=.o))
$(util_c_objects): %.o: %.c $(util_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(util_a): $(util_c_objects)
	$(AR) $(ARFLAGS) $@ $^
util: $(util_a)
util-clean:
	rm -f $(util_c_objects) $(util_a)


##
## Linear algebra
##
linear_algebra_a = src/linear_algebra/linear_algebra.a
linear_algebra_c_sources = \
	src/linear_algebra/sparse/csr_matrix_int32.c \
	src/linear_algebra/sparse/csr_matrix_int32_int32.c \
	src/linear_algebra/sparse/ellpack_matrix_int32.c \
	src/linear_algebra/tensor.c \
	src/linear_algebra/tensor_mode.c
linear_algebra_c_headers = \
	src/linear_algebra/sparse/csr_matrix_int32.h \
	src/linear_algebra/sparse/csr_matrix_int32_int32.h \
	src/linear_algebra/sparse/ellpack_matrix_int32.h \
	src/linear_algebra/matrix/matrix12x4.h \
	src/linear_algebra/matrix/matrix12x4d256.h \
	src/linear_algebra/matrix/matrix16x4.h \
	src/linear_algebra/matrix/matrix16x4d256.h \
	src/linear_algebra/matrix/matrix16x8.h \
	src/linear_algebra/matrix/matrix16x8d512.h \
	src/linear_algebra/matrix/matrix2x2.h \
	src/linear_algebra/matrix/matrix3x3.h \
	src/linear_algebra/matrix/matrix3x3d256.h \
	src/linear_algebra/matrix/matrix3x3s.h \
	src/linear_algebra/matrix/matrix3x3sd256.h \
	src/linear_algebra/matrix/matrix3x3sx4.h \
	src/linear_algebra/matrix/matrix3x3sx4d256.h \
	src/linear_algebra/matrix/matrix3x3sx8.h \
	src/linear_algebra/matrix/matrix3x3sx8d512.h \
	src/linear_algebra/matrix/matrix3x3x4.h \
	src/linear_algebra/matrix/matrix3x3x4d256.h \
	src/linear_algebra/matrix/matrix3x3x8.h \
	src/linear_algebra/matrix/matrix3x3x8d512.h \
	src/linear_algebra/matrix/matrix4x8.h \
	src/linear_algebra/matrix/matrix4x12.h \
	src/linear_algebra/matrix/matrix4x12d256.h \
	src/linear_algebra/matrix/matrix4x16.h \
	src/linear_algebra/matrix/matrix4x16d256.h \
	src/linear_algebra/matrix/matrix4x4.h \
	src/linear_algebra/matrix/matrix4x4d256.h \
	src/linear_algebra/matrix/matrix8x12.h \
	src/linear_algebra/matrix/matrix8x12d512.h \
	src/linear_algebra/matrix/matrix8x16.h \
	src/linear_algebra/matrix/matrix8x16d512.h \
	src/linear_algebra/matrix/matrix8x4.h \
	src/linear_algebra/matrix/matrix8x8.h \
	src/linear_algebra/matrix/matrix8x8d512.h \
	src/linear_algebra/tensor.h \
	src/linear_algebra/tensor_mode.h \
	src/linear_algebra/vector/vector16.h \
	src/linear_algebra/vector/vector16d256.h \
	src/linear_algebra/vector/vector2.h \
	src/linear_algebra/vector/vector3.h \
	src/linear_algebra/vector/vector3d256.h \
	src/linear_algebra/vector/vector3x4.h \
	src/linear_algebra/vector/vector3x4d256.h \
	src/linear_algebra/vector/vector3x8.h \
	src/linear_algebra/vector/vector3x8d512.h \
	src/linear_algebra/vector/vector4.h \
	src/linear_algebra/vector/vector4d256.h \
	src/linear_algebra/vector/vector8.h \
	src/linear_algebra/vector/vector8d512.h
linear_algebra_c_objects := $(foreach x,$(linear_algebra_c_sources),$(x:.c=.o))
$(linear_algebra_c_objects): %.o: %.c $(linear_algebra_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(linear_algebra_a): $(linear_algebra_c_objects)
	$(AR) $(ARFLAGS) $@ $^
linear_algebra: $(linear_algebra_a)
linear_algebra-clean:
	rm -f $(linear_algebra_c_objects) $(linear_algebra_a)


##
## Linear solvers
##
linear_solver_a = src/linear_solver/linear_solver.a
linear_solver_c_sources = src/linear_solver/ksp_solver.c
linear_solver_c_headers = src/linear_solver/ksp_solver.h
linear_solver_c_objects := $(foreach x,$(linear_solver_c_sources),$(x:.c=.o))
$(linear_solver_c_objects): %.o: %.c $(linear_solver_c_headers)
	$(CC) -c $(CFLAGS) $(PETSC_INCLUDES) $(INCLUDES) $< -o $@
$(linear_solver_a): $(linear_solver_c_objects)
	$(AR) $(ARFLAGS) $@ $^
linear_solver: $(linear_solver_a)
linear_solver-clean:
	rm -f $(linear_solver_c_objects) $(linear_solver_a)


##
## Meshes
##
mesh_a = src/mesh/mesh.a
mesh_c_sources = \
	src/mesh/cells/empty.c \
	src/mesh/cells/simplex0.c \
	src/mesh/cells/simplex1.c \
	src/mesh/cells/simplex2.c \
	src/mesh/cells/simplex3.c \
	src/mesh/cells/simplex3x4.c \
	src/mesh/cells/simplex3x8.c \
	src/mesh/mesh.c \
	src/mesh/mesh_entities.c \
	src/mesh/mesh_entity_function.c \
	src/mesh/mesh_entity_map.c \
	src/mesh/mesh_entity_types.c \
	src/mesh/mesh_function.c \
	src/mesh/mesh_partitioning.c \
	src/mesh/mesh_reordering.c \
	src/mesh/simplex.c \
	src/mesh/singleton_mesh.c \
	src/mesh/submesh.c \
	src/mesh/tetgen_mesh.c \
	src/mesh/unit_mesh.c
mesh_c_headers = \
	src/mesh/cells/empty.h \
	src/mesh/cells/simplex0.h \
	src/mesh/cells/simplex1.h \
	src/mesh/cells/simplex2.h \
	src/mesh/cells/simplex3.h \
	src/mesh/cells/simplex3x4.h \
	src/mesh/cells/simplex3x8.h \
	src/mesh/mesh.h \
	src/mesh/mesh_entities.h \
	src/mesh/mesh_entity_function.h \
	src/mesh/mesh_entity_map.h \
	src/mesh/mesh_entity_types.h \
	src/mesh/mesh_function.h \
	src/mesh/mesh_partitioning.h \
	src/mesh/mesh_reordering.h \
	src/mesh/simplex.h \
	src/mesh/singleton_mesh.h \
	src/mesh/submesh.h \
	src/mesh/tetgen_mesh.h \
	src/mesh/unit_mesh.h
mesh_c_objects := $(foreach x,$(mesh_c_sources),$(x:.c=.o))
$(mesh_c_objects): %.o: %.c $(mesh_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(mesh_a): $(mesh_c_objects)
	$(AR) $(ARFLAGS) $@ $^
mesh: $(mesh_a)
mesh-clean:
	rm -f $(mesh_c_objects) $(mesh_a)


##
## Finite element methods
##
fem_a = src/fem/fem.a
fem_c_sources = \
	src/fem/finite_element_space.c \
	src/fem/finite_element_space_function.c \
	src/fem/finite_element_types.c \
	src/fem/lagrange.c
fem_c_headers = \
	src/fem/finite_element_space.h \
	src/fem/finite_element_space_function.h \
	src/fem/finite_element_types.h \
	src/fem/lagrange.h \
	src/fem/finite_elements/lagrange.h \
	src/fem/finite_elements/simplex2_p1.h \
	src/fem/finite_elements/simplex2_p1_simplex2_p1.h \
	src/fem/finite_elements/simplex2_p2.h \
	src/fem/finite_elements/simplex2_p2_simplex2_p2.h \
	src/fem/finite_elements/simplex3_p1.h \
	src/fem/finite_elements/simplex3_p1_simplex3_p1.h \
	src/fem/finite_elements/simplex3_p2.h \
	src/fem/finite_elements/simplex3_p2_simplex3_p2.h \
	src/fem/finite_elements/simplex3x4_p1.h \
	src/fem/finite_elements/simplex3x4_p1_simplex3x4_p1.h \
	src/fem/finite_elements/simplex3x4_p2.h \
	src/fem/finite_elements/simplex3x4_p2_simplex3x4_p2.h \
	src/fem/finite_elements/simplex3x8_p1.h \
	src/fem/finite_elements/simplex3x8_p1_simplex3x8_p1.h \
	src/fem/finite_elements/simplex3x8_p2.h \
	src/fem/finite_elements/simplex3x8_p2_simplex3x8_p2.h
fem_c_objects := $(foreach x,$(fem_c_sources),$(x:.c=.o))
$(fem_c_objects): %.o: %.c $(fem_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(fem_a): $(fem_c_objects)
	$(AR) $(ARFLAGS) $@ $^
fem: $(fem_a)
fem-clean:
	rm -f $(fem_c_objects) $(fem_a)


##
## Variational forms
##
variational_form_a = src/variational_form/variational_form.a
variational_form_c_sources = \
	src/variational_form/assembly/assembly_algorithm.c \
	src/variational_form/bilinear_form/laplacian.c \
	src/variational_form/bilinear_form/mass_matrix.c \
	src/variational_form/bilinear_form/scatter_to_global_matrix.c \
	src/variational_form/coefficient_function.c \
	src/variational_form/functional/coefficient_integral.c \
	src/variational_form/functional/coefficient_l2_norm.c \
	src/variational_form/functional/gather_vertex_coordinates.c \
	src/variational_form/functional/transform_to_reference_cell.c \
	src/variational_form/functional/volume.c \
	src/variational_form/linear_form/argument_integral.c \
	src/variational_form/linear_form/load_vector.c \
	src/variational_form/linear_form/pointwise_projection.c \
	src/variational_form/linear_form/scatter_to_global_vector.c \
	src/variational_form/variational_form.c \
	src/variational_form/variational_form_types.c
variational_form_c_headers = \
	src/variational_form/assembly/assemble_cellwise_atomic.h \
	src/variational_form/assembly/assemble_cellwise_sequential.h \
	src/variational_form/assembly/assemble_rowwise.h \
	src/variational_form/assembly/assembly_algorithm.h \
	src/variational_form/assembly/element_matrix_cellwise_sequential.h \
	src/variational_form/bilinear_form/bilinear_forms.h \
	src/variational_form/coefficient_function.h \
	src/variational_form/functional/functionals.h \
	src/variational_form/linear_form/linear_forms.h \
	src/variational_form/variational_form.h \
	src/variational_form/variational_form_types.h
variational_form_c_objects := $(foreach x,$(variational_form_c_sources),$(x:.c=.o))
$(variational_form_c_objects): %.o: %.c $(variational_form_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(variational_form_a): $(variational_form_c_objects)
	$(AR) $(ARFLAGS) $@ $^
variational_form: $(variational_form_a)
variational_form-clean:
	rm -f $(variational_form_c_objects) $(variational_form_a)


##
## Solvers for variational problems based on finite element methods
##
solvers_a = src/solvers/solvers.a
solvers_c_sources = \
	src/solvers/poisson.c \
	src/solvers/projection.c
solvers_c_headers = \
	src/solvers/poisson.h \
	src/solvers/projection.h
solvers_c_objects := $(foreach x,$(solvers_c_sources),$(x:.c=.o))
$(solvers_c_objects): %.o: %.c $(solvers_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(solvers_a): $(solvers_c_objects)
	$(AR) $(ARFLAGS) $@ $^
solvers: $(solvers_a)
solvers-clean:
	rm -f $(solvers_c_objects) $(solvers_a)


##
## Benchmarks
##
benchmark_binaries := $(addprefix $(bindir)/,$(benchmarks))
$(benchmark_binaries): | $(bindir)
$(bindir):
	mkdir $(bindir)
$(benchmarks) : $(benchmark_binaries)


##
## Benchmarks for finite element algorithms
##
benchmarks_common_a = src/benchmarks/common/benchmarks_common.a
benchmarks_common_c_sources = \
	src/benchmarks/common/benchmark_options.c \
	src/benchmarks/common/mesh_options.c
benchmarks_common_c_headers = \
	src/benchmarks/common/benchmark_options.h \
	src/benchmarks/common/mesh_options.h
benchmarks_common_c_objects := $(foreach x,$(benchmarks_common_c_sources),$(x:.c=.o))
$(benchmarks_common_c_objects): %.o: %.c $(benchmarks_common_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(benchmarks_common_a): $(benchmarks_common_c_objects)
	$(AR) $(ARFLAGS) $@ $^
benchmarks_common: $(benchmarks_common_a)
benchmarks_common-clean:
	rm -f $(benchmarks_common_c_objects) $(benchmarks_common_a)


# 2D Poisson solver
poisson2d_cg = $(bindir)/poisson2d_cg
poisson2d_cg_c_sources = \
	src/benchmarks/poisson2d/program_options.c \
	src/benchmarks/poisson2d/poisson2d_cg_seq.c
poisson2d_cg_c_headers = src/benchmarks/poisson2d/program_options.h
poisson2d_cg_c_objects := $(foreach x,$(poisson2d_cg_c_sources),$(x:.c=.o))
$(poisson2d_cg_c_objects): %.o: %.c $(poisson2d_cg_c_headers)
	$(CC) -c $(CFLAGS) $(PETSC_INCLUDES) $(INCLUDES) $< -o $@
$(poisson2d_cg): $(poisson2d_cg_c_objects) $(solvers_a) $(variational_form_a) $(fem_a) $(mesh_a) $(linear_solver_a) $(linear_algebra_a) $(util_a)
	$(CC) $(CFLAGS) $(INCLUDES) $^ $(LDFLAGS) $(PETSC_LIBS) $(MPI_LIBS) -lm -o $@
poisson2d_cg-clean:
	rm -f $(poisson2d_cg_c_objects) $(poisson2d_cg)


# 3D Poisson solver
poisson3d = $(bindir)/poisson3d
poisson3d_c_sources = \
	src/benchmarks/poisson3d/program_options.c \
	src/benchmarks/poisson3d/poisson3d.c
poisson3d_c_headers = src/benchmarks/poisson3d/program_options.h
poisson3d_c_objects := $(foreach x,$(poisson3d_c_sources),$(x:.c=.o))
$(poisson3d_c_objects): %.o: %.c $(poisson3d_c_headers)
	$(CC) -c $(CFLAGS) $(PETSC_INCLUDES) $(INCLUDES) $< -o $@
$(poisson3d): $(poisson3d_c_objects) $(solvers_a) $(variational_form_a) $(fem_a) $(mesh_a) $(linear_solver_a) $(linear_algebra_a) $(util_a)
	$(CC) $(CFLAGS) $(INCLUDES) $^ $(LDFLAGS) $(PETSC_LIBS) $(MPI_LIBS) -lm -o $@
poisson3d-clean:
	rm -f $(poisson3d_c_objects) $(poisson3d)

# Benchmark for assembly of finite element matrices
assembly = $(bindir)/assembly
assembly_c_sources = \
	src/benchmarks/assembly/program_options.c \
	src/benchmarks/assembly/assembly_kernel.c \
	src/benchmarks/assembly/assembly.c
assembly_c_headers = \
	src/benchmarks/assembly/program_options.h \
	src/benchmarks/assembly/assembly_kernel.h
assembly_c_objects := $(foreach x,$(assembly_c_sources),$(x:.c=.o))
$(assembly_c_objects): %.o: %.c $(assembly_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(assembly): $(assembly_c_objects) $(benchmarks_common_a) $(solvers_a) $(variational_form_a) $(fem_a) $(mesh_a) $(linear_solver_a) $(linear_algebra_a) $(util_a)
	$(CC) $(CFLAGS) $(INCLUDES) $^ $(LDFLAGS) $(PETSC_LIBS) $(MPI_LIBS) -lm -o $@
assembly-clean:
	rm -f $(assembly_c_objects) $(assembly)

# A STREAM-like benchmark for measuring memory throughput
stream = $(bindir)/stream
stream_c_sources = \
	src/benchmarks/stream/program_options.c \
	src/benchmarks/stream/stream_kernel.c \
	src/benchmarks/stream/stream.c
stream_c_headers = \
	src/benchmarks/stream/stream_kernel.h \
	src/benchmarks/stream/program_options.h
stream_c_objects := $(foreach x,$(stream_c_sources),$(x:.c=.o))
$(stream_c_objects): %.o: %.c $(stream_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(stream): $(stream_c_objects) $(benchmarks_common_a) $(util_a)
	$(CC) $(CFLAGS) $(INCLUDES) $^ $(LDFLAGS) $(PETSC_LIBS) $(MPI_LIBS) -lm -o $@
stream-clean:
	rm -f $(stream_c_objects) $(stream)

# Micro-benchmarks for gather- and scatter-type memory accesses
gather = $(bindir)/gather
gather_c_sources = \
	src/benchmarks/gather/program_options.c \
	src/benchmarks/gather/gather_kernel.c \
	src/benchmarks/gather/gather.c
gather_c_headers = \
	src/benchmarks/gather/gather_kernel.h \
	src/benchmarks/gather/program_options.h
gather_c_objects := $(foreach x,$(gather_c_sources),$(x:.c=.o))
$(gather_c_objects): %.o: %.c $(gather_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(gather): $(gather_c_objects) $(benchmarks_common_a) $(util_a)
	$(CC) $(CFLAGS) $(INCLUDES) $^ $(LDFLAGS) $(PETSC_LIBS) $(MPI_LIBS) -lm -o $@
gather-clean:
	rm -f $(gather_c_objects) $(gather)

# Micro-benchmarks for pointer-chasing memory accesses
pointer_chase = $(bindir)/pointer_chase
pointer_chase_c_sources = \
	src/benchmarks/pointer_chase/program_options.c \
	src/benchmarks/pointer_chase/pointer_chase_kernel.c \
	src/benchmarks/pointer_chase/pointer_chase.c
pointer_chase_c_headers = \
	src/benchmarks/pointer_chase/pointer_chase_kernel.h \
	src/benchmarks/pointer_chase/program_options.h
pointer_chase_c_objects := $(foreach x,$(pointer_chase_c_sources),$(x:.c=.o))
$(pointer_chase_c_objects): %.o: %.c $(pointer_chase_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(pointer_chase): $(pointer_chase_c_objects) $(benchmarks_common_a) $(util_a)
	$(CC) $(CFLAGS) $(INCLUDES) $^ $(LDFLAGS) $(PETSC_LIBS) $(MPI_LIBS) -lm -o $@
pointer_chase-clean:
	rm -f $(pointer_chase_c_objects) $(pointer_chase)


##
## Utility programs
##
reorder_mesh = $(bindir)/reorder_mesh
reorder_mesh_c_sources = \
	src/benchmarks/reorder_mesh/mesh_reordering_algorithm.c \
	src/benchmarks/reorder_mesh/program_options.c \
	src/benchmarks/reorder_mesh/reorder_mesh.c
reorder_mesh_c_headers = \
	src/benchmarks/reorder_mesh/mesh_reordering_algorithm.h \
	src/benchmarks/reorder_mesh/program_options.h
reorder_mesh_c_objects := $(foreach x,$(reorder_mesh_c_sources),$(x:.c=.o))
$(reorder_mesh_c_objects): %.o: %.c $(reorder_mesh_c_headers)
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(reorder_mesh): $(reorder_mesh_c_objects) $(benchmarks_common_a) $(solvers_a) $(variational_form_a) $(fem_a) $(mesh_a) $(linear_solver_a) $(linear_algebra_a) $(util_a)
	$(CC) $(CFLAGS) $(INCLUDES) $^ $(LDFLAGS) $(METIS_LIBS) -lm -o $@
reorder_mesh-clean:
	rm -f $(reorder_mesh_c_objects) $(reorder_mesh)


##
## Unit tests
##
.PHONY: run-unit-tests
run-unit-tests: $(unittest)
	OMP_NUM_THREADS=1 $(unittest)

unittest_sources = \
	test/fem/finite_element_spaces/test_empty.cpp \
	test/fem/finite_element_spaces/test_simplex0.cpp \
	test/fem/finite_element_spaces/test_simplex1.cpp \
	test/fem/finite_element_spaces/test_simplex2.cpp \
	test/fem/finite_element_spaces/test_simplex3.cpp \
	test/fem/test_finite_element_types.cpp \
	test/fem/test_lagrange.cpp \
	test/linear_algebra/matrix/test_matrix12x4.cpp \
	test/linear_algebra/matrix/test_matrix12x4d256.cpp \
	test/linear_algebra/matrix/test_matrix12x8.cpp \
	test/linear_algebra/matrix/test_matrix16x4.cpp \
	test/linear_algebra/matrix/test_matrix16x4d256.cpp \
	test/linear_algebra/matrix/test_matrix16x8.cpp \
	test/linear_algebra/matrix/test_matrix16x8d512.cpp \
	test/linear_algebra/matrix/test_matrix2x2.cpp \
	test/linear_algebra/matrix/test_matrix3x3.cpp \
	test/linear_algebra/matrix/test_matrix3x3x4.cpp \
	test/linear_algebra/matrix/test_matrix3x3x8.cpp \
	test/linear_algebra/matrix/test_matrix4x12.cpp \
	test/linear_algebra/matrix/test_matrix4x16.cpp \
	test/linear_algebra/matrix/test_matrix4x16d256.cpp \
	test/linear_algebra/matrix/test_matrix4x4.cpp \
	test/linear_algebra/matrix/test_matrix4x4d256.cpp \
	test/linear_algebra/matrix/test_matrix4x4d512.cpp \
	test/linear_algebra/matrix/test_matrix4x8.cpp \
	test/linear_algebra/matrix/test_matrix4x8d512.cpp \
	test/linear_algebra/matrix/test_matrix8x12.cpp \
	test/linear_algebra/matrix/test_matrix8x12d512.cpp \
	test/linear_algebra/matrix/test_matrix8x16.cpp \
	test/linear_algebra/matrix/test_matrix8x16d512.cpp \
	test/linear_algebra/matrix/test_matrix8x4.cpp \
	test/linear_algebra/matrix/test_matrix8x8.cpp \
	test/linear_algebra/sparse/test_csr_matrix_int32.cpp \
	test/linear_algebra/sparse/test_csr_matrix_int32_int32.cpp \
	test/linear_algebra/sparse/test_ellpack_matrix_int32.cpp \
	test/linear_algebra/test_tensor.cpp \
	test/linear_algebra/test_tensor_mode.cpp \
	test/linear_algebra/vector/test_vector2.cpp \
	test/linear_algebra/vector/test_vector3.cpp \
	test/linear_algebra/vector/test_vector3x4.cpp \
	test/main.cpp \
	test/mesh/cells/test_empty.cpp \
	test/mesh/cells/test_simplex0.cpp \
	test/mesh/cells/test_simplex1.cpp \
	test/mesh/cells/test_simplex2.cpp \
	test/mesh/cells/test_simplex3.cpp \
	test/mesh/cells/test_simplex3x4.cpp \
	test/mesh/cells/test_simplex3x8.cpp \
	test/mesh/test_mesh.cpp \
	test/mesh/test_mesh_entities.cpp \
	test/mesh/test_mesh_entity_function.cpp \
	test/mesh/test_mesh_entity_map.cpp \
	test/mesh/test_mesh_entity_types.cpp \
	test/mesh/test_mesh_function.cpp \
	test/mesh/test_mesh_partitioning.cpp \
	test/mesh/test_mesh_reordering.cpp \
	test/mesh/test_simplex.cpp \
	test/mesh/test_submesh.cpp \
	test/mesh/test_singleton_mesh.cpp \
	test/mesh/test_unit_mesh.cpp \
	test/solvers/test_projection.cpp \
	test/test_utils/common_finite_element_spaces.cpp \
	test/test_utils/common_meshes.cpp \
	test/util/test_binary.cpp \
	test/util/test_combinations.cpp \
	test/util/test_cuthill_mckee.cpp \
	test/util/test_matrix_market.cpp \
	test/util/test_parse.cpp \
	test/util/test_partition.cpp \
	test/util/test_perf_events.cpp \
	test/util/test_perf_sessions.cpp \
	test/util/test_permutation.cpp \
	test/util/test_power_of_two_radix.cpp \
	test/util/test_prefix_sum.cpp \
	test/util/test_sample_statistics.cpp \
	test/util/test_search.cpp \
	test/util/test_selections.cpp \
	test/util/test_sort.cpp \
	test/util/test_stream_compaction.cpp \
	test/variational_form/test_argument_integral.cpp \
	test/variational_form/test_bilinear_form.cpp \
	test/variational_form/test_coefficient_integral.cpp \
	test/variational_form/test_coefficient_l2_norm.cpp \
	test/variational_form/test_functional.cpp \
	test/variational_form/test_gather_vertex_coordinates.cpp \
	test/variational_form/test_transform_to_reference_cell.cpp \
	test/variational_form/test_laplacian.cpp \
	test/variational_form/test_linear_form.cpp \
	test/variational_form/test_load_vector.cpp \
	test/variational_form/test_mass_matrix.cpp \
	test/variational_form/test_pointwise_projection.cpp \
	test/variational_form/test_scatter_to_global_matrix.cpp \
	test/variational_form/test_scatter_to_global_vector.cpp \
	test/variational_form/test_variational_form_types.cpp \
	test/variational_form/test_volume.cpp
unittest_headers = \
	test/test_utils/common_finite_element_spaces.hpp \
	test/test_utils/common_meshes.hpp
unittest_objects := $(foreach x,$(unittest_sources),$(x:.cpp=.o))
$(unittest_objects): %.o: %.cpp $(unittest_headers)
	$(CXX) $(CPPFLAGS) -c $(CXXFLAGS) $(INCLUDES) -iquote test $(GTEST_INCLUDES) $< -o $@
$(unittest): $(unittest_objects) $(solvers_a) $(variational_form_a) $(fem_a) $(mesh_a) $(linear_solver_a) $(linear_algebra_a) $(util_a)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) $(GTEST_LIBS) $(PETSC_LIBS) $(MPI_LIBS) $(METIS_LIBS) -lm -o $@
$(unittest)-clean:
	rm -f $(unittest_objects) $(unittest)
