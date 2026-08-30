#!/bin/sh
# =====================================================================
# Builds cppDE's two optional C libraries from upstream sources into a
# per-user prefix:
#
#   suitesparse   SuiteSparse/KLU. Required by the sparse Jacobian path
#                 of cppODE() and cvode().
#   sundials      SUNDIALS/CVODES. Required by the cvode() backend.
#                 Implies suitesparse: SUNDIALS links its sparse solver
#                 wrappers against an external KLU at its own build time.
#
# Usage: build-libs.sh [suitesparse|sundials]      (default: sundials)
#
# Invoked by ./configure for CPPDE_BUILD_SUITESPARSE / CPPDE_BUILD_SUNDIALS,
# and by the R helper cppDE::install_libs(). Requires network access and
# runs for several minutes, so it is never started implicitly.
#
# Both components install into one prefix keyed by the pinned version
# pair. Each component is skipped when already present, so a suitesparse
# build can later be completed to a full sundials build in place.
#
# Prints the prefix on stdout; all logging goes to stderr.
# =====================================================================
set -e

component="${1:-sundials}"
case "$component" in
  suitesparse) build_suitesparse=1; build_sundials=0 ;;
  sundials)    build_suitesparse=1; build_sundials=1 ;;
  *) echo "build-libs: unknown component '$component'" >&2; exit 2 ;;
esac

# Upstream releases. SUNDIALS compiles its KLU wrappers against whichever
# SuiteSparse it is pointed at, so the two pins are chosen as a pair.
SUNDIALS_VERSION="${CPPDE_SUNDIALS_VERSION:-7.4.0}"
SUITESPARSE_VERSION="${CPPDE_SUITESPARSE_VERSION:-7.10.0}"

# Additional CMake options, appended to the respective configure step and
# therefore overriding the defaults set below. Word-split on whitespace.
# Example: CPPDE_SUNDIALS_CMAKE_ARGS="-DCMAKE_C_COMPILER=icx"
SUITESPARSE_CMAKE_ARGS="${CPPDE_SUITESPARSE_CMAKE_ARGS:-}"
SUNDIALS_CMAKE_ARGS="${CPPDE_SUNDIALS_CMAKE_ARGS:-}"

log() { echo "build-libs: $*" >&2; }
die() { log "ERROR: $*"; exit 1; }

# R CMD INSTALL exports the compiler and flag settings it uses for R
# packages (CC, CFLAGS, FC, FFLAGS, LDFLAGS, ...). They are inappropriate
# for building an unrelated C library and break CMake's own probes -- the
# Fortran name-mangling test SUNDIALS needs for its LAPACK interface fails
# under them. Clearing them makes this build behave identically whether it
# is started from ./configure or from a shell. Compiler choices can still
# be passed through CPPDE_*_CMAKE_ARGS.
unset CC CXX CPP FC F77 F90 F95 OBJC OBJCXX \
      CFLAGS CXXFLAGS CPPFLAGS FFLAGS FCFLAGS OBJCFLAGS \
      LDFLAGS LIBS AR RANLIB NM MAKEFLAGS MFLAGS

# tools::R_user_dir() is the per-user data location prescribed by CRAN
# policy. CPPDE_LIBS_CACHE overrides it for filesystems chosen by the
# user (HPC work directories, persistent container layers).
cache_root="${CPPDE_LIBS_CACHE:-}"
if [ -z "$cache_root" ]; then
  cache_root=$(R --slave --no-save -e 'cat(tools::R_user_dir("cppDE", "cache"))' 2>/dev/null || true)
fi
if [ -z "$cache_root" ]; then
  cache_root="${XDG_CACHE_HOME:-$HOME/.cache}/cppDE"
fi
PREFIX="$cache_root/deps-sundials-${SUNDIALS_VERSION}-suitesparse-${SUITESPARSE_VERSION}"

# Completion is recorded by a stamp file written after the component has
# been installed and verified. Individual headers or libraries are not a
# safe marker: an interrupted `cmake --install` can leave a subset of them
# behind, which would then be mistaken for a finished build.
have_suitesparse() { [ -f "$PREFIX/.suitesparse-complete" ]; }
have_sundials()    { [ -f "$PREFIX/.sundials-complete" ]; }

files_suitesparse() {
  [ -f "$PREFIX/include/suitesparse/klu.h" ] && ls "$PREFIX"/lib/libklu.* >/dev/null 2>&1
}
files_sundials() {
  [ -f "$PREFIX/include/cvodes/cvodes.h" ] \
    && [ -f "$PREFIX/include/sundials/sundials_config.h" ] \
    && ls "$PREFIX"/lib/libsundials_cvodes.* >/dev/null 2>&1
}

if have_suitesparse; then build_suitesparse=0; fi
if [ "$build_sundials" = "1" ] && have_sundials; then build_sundials=0; fi

if [ "$build_suitesparse" = "0" ] && [ "$build_sundials" = "0" ]; then
  log "$component already present in $PREFIX"
  log "delete that directory to force a rebuild"
  printf '%s\n' "$PREFIX"
  exit 0
fi

command -v cmake >/dev/null 2>&1 || die "cmake not found; install cmake (>= 3.18) and retry"

fetcher=""
if command -v curl >/dev/null 2>&1; then
  fetcher="curl"
elif command -v wget >/dev/null 2>&1; then
  fetcher="wget"
elif command -v git >/dev/null 2>&1; then
  fetcher="git"
else
  die "need one of curl, wget or git to obtain the sources"
fi

njobs=$( (command -v nproc >/dev/null 2>&1 && nproc) \
        || (command -v sysctl >/dev/null 2>&1 && sysctl -n hw.ncpu) \
        || echo 2)
# One core is left to the rest of the machine, which matters on shared
# login nodes.
if [ "$njobs" -gt 2 ]; then
  njobs=$((njobs - 1))
fi

work=$(mktemp -d 2>/dev/null || mktemp -d -t cppde-libs)

# On a normal exit only the scratch directory goes away. On abort the
# prefix is discarded as well: this run was mid-build, so whatever is in
# there is partial, and a partial prefix must not survive to confuse the
# next attempt. A handler without an explicit exit would let the script
# continue, so this one terminates. Signal set as in the autoconf
# boilerplate, spelled with names because R CMD check flags numeric
# signals as a possible bashism.
on_exit()  { rm -rf "$work"; }
on_abort() {
  trap - HUP INT PIPE TERM
  log "aborted; discarding the incomplete build in $PREFIX"
  rm -rf "$work" "$PREFIX"
  exit 130
}
trap on_exit EXIT
trap on_abort HUP INT PIPE TERM

# $1 tarball URL, $2 git URL, $3 git tag, $4 destination directory name
fetch_source() {
  url="$1"; git_url="$2"; git_tag="$3"; dest="$4"
  case "$fetcher" in
    curl|wget)
      log "downloading $url"
      if [ "$fetcher" = "curl" ]; then
        curl -fsSL "$url" -o "$work/$dest.tar.gz" || die "download failed: $url"
      else
        wget -q "$url" -O "$work/$dest.tar.gz" || die "download failed: $url"
      fi
      mkdir -p "$work/$dest"
      # Release tarballs wrap their contents in one versioned top-level
      # directory, which --strip-components removes.
      tar -xzf "$work/$dest.tar.gz" -C "$work/$dest" --strip-components=1 \
        || die "could not unpack $dest"
      ;;
    git)
      log "cloning $git_url at $git_tag"
      git clone --depth 1 -b "$git_tag" "$git_url" "$work/$dest" >&2 \
        || die "git clone failed: $git_url"
      ;;
  esac
}

log "prefix: $PREFIX"
log "cores: $njobs"
if [ -n "$SUITESPARSE_CMAKE_ARGS" ]; then
  log "extra SuiteSparse CMake options: $SUITESPARSE_CMAKE_ARGS"
fi
if [ -n "$SUNDIALS_CMAKE_ARGS" ]; then
  log "extra SUNDIALS CMake options: $SUNDIALS_CMAKE_ARGS"
fi

# ---------------------------------------------------------------------
# SuiteSparse
#
# Only the modules KLU requires are enabled; the full project also builds
# CHOLMOD, UMFPACK and GraphBLAS, which cppDE never calls.
#
# CMAKE_INSTALL_LIBDIR is pinned to lib because CMake otherwise selects
# lib64 on Fedora/RHEL/SUSE and lib elsewhere.
# ---------------------------------------------------------------------
if [ "$build_suitesparse" = "1" ]; then
  log "building SuiteSparse $SUITESPARSE_VERSION (klu, amd, colamd, btf)"
  fetch_source \
    "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v${SUITESPARSE_VERSION}.tar.gz" \
    "https://github.com/DrTimothyAldenDavis/SuiteSparse.git" \
    "v${SUITESPARSE_VERSION}" \
    suitesparse

  cmake -S "$work/suitesparse" -B "$work/suitesparse/build" \
        -DCMAKE_INSTALL_PREFIX="$PREFIX" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DBUILD_SHARED_LIBS=ON \
        -DBUILD_STATIC_LIBS=OFF \
        -DSUITESPARSE_ENABLE_PROJECTS="suitesparse_config;amd;colamd;btf;klu" \
        -DSUITESPARSE_USE_OPENMP=OFF \
        -DSUITESPARSE_USE_CUDA=OFF \
        $SUITESPARSE_CMAKE_ARGS \
        >&2 || die "SuiteSparse configure failed"
  cmake --build "$work/suitesparse/build" -j "$njobs" >&2 || die "SuiteSparse build failed"
  cmake --install "$work/suitesparse/build" >&2 || die "SuiteSparse install failed"
  files_suitesparse || die "SuiteSparse installed but $PREFIX is incomplete"
  : > "$PREFIX/.suitesparse-complete"
fi

# ---------------------------------------------------------------------
# SUNDIALS
#
# KLU_INCLUDE_DIR / KLU_LIBRARY_DIR point at the SuiteSparse in this same
# prefix, so the wrappers are compiled against the copy they are linked
# against when a model is built.
#
# MPI is disabled: a SUNDIALS built with SUNDIALS_MPI_ENABLED=1 makes
# <sundials/sundials_types.h> include <mpi.h>, which would then have to
# resolve on every machine that compiles a model.
# ---------------------------------------------------------------------
configure_sundials() {
  # "$@": extra options for this attempt. Passed through unsplit -- a
  # -DLAPACK_LIBRARIES holding several linker flags is one argument, and
  # letting the shell split it makes cmake abort with "Unknown argument".
  rm -rf "$work/sundials/build"
  cmake -S "$work/sundials" -B "$work/sundials/build" \
        -DCMAKE_INSTALL_PREFIX="$PREFIX" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DBUILD_SHARED_LIBS=ON \
        -DBUILD_STATIC_LIBS=OFF \
        -DBUILD_CVODES=ON \
        -DBUILD_ARKODE=OFF \
        -DBUILD_IDA=OFF \
        -DBUILD_IDAS=OFF \
        -DBUILD_KINSOL=OFF \
        -DENABLE_MPI=OFF \
        -DENABLE_OPENMP=OFF \
        -DEXAMPLES_ENABLE_C=OFF \
        -DEXAMPLES_ENABLE_CXX=OFF \
        -DEXAMPLES_INSTALL=OFF \
        -DENABLE_KLU=ON \
        -DKLU_INCLUDE_DIR="$PREFIX/include/suitesparse" \
        -DKLU_LIBRARY_DIR="$PREFIX/lib" \
        "$@" \
        $SUNDIALS_CMAKE_ARGS \
        >&2
}

# CMake's FindLAPACK discards a caller-supplied LAPACK_LIBRARIES -- it runs
# `set(LAPACK_LIBRARIES)` before searching -- and then picks whatever BLAS it
# finds on the system. On a cluster that is rarely the one R uses: R links
# MKL loaded through a module, while FindLAPACK either fails outright (making
# SUNDIALS' `find_package(LAPACK REQUIRED)` fatal) or picks up an unrelated
# reference LAPACK. Shadowing the module is the only way to pin the choice.
#
# $1: the linker flags to link LAPACK and BLAS with.
write_lapack_shim() {
  mkdir -p "$work/cmake-modules"
  # CMake wants a ;-separated list; R reports space-separated linker flags.
  _shim_libs=$(printf '%s' "$1" | tr -s ' ' ';')
  for _shim_mod in LAPACK BLAS; do
    cat > "$work/cmake-modules/Find$_shim_mod.cmake" <<EOF
# Generated by cppDE's build-libs.sh -- pins BLAS/LAPACK to R's own.
set(LAPACK_LINKER_FLAGS "")
set(BLAS_LINKER_FLAGS "")
set(LAPACK_LIBRARIES "$_shim_libs")
set(BLAS_LIBRARIES "$_shim_libs")
set(LAPACK_FOUND TRUE)
set(BLAS_FOUND TRUE)
foreach(_lib BLAS LAPACK)
  if(NOT TARGET \${_lib}::\${_lib})
    add_library(\${_lib}::\${_lib} INTERFACE IMPORTED)
    set_target_properties(\${_lib}::\${_lib} PROPERTIES
                          INTERFACE_LINK_LIBRARIES "$_shim_libs")
  endif()
endforeach()
EOF
  done
}

if [ "$build_sundials" = "1" ]; then
  log "building SUNDIALS $SUNDIALS_VERSION (cvodes, dense + KLU solvers)"
  fetch_source \
    "https://github.com/LLNL/sundials/releases/download/v${SUNDIALS_VERSION}/sundials-${SUNDIALS_VERSION}.tar.gz" \
    "https://github.com/LLNL/sundials.git" \
    "v${SUNDIALS_VERSION}" \
    sundials

  # SUNDIALS' LAPACK-backed dense solver is built against the same BLAS and
  # LAPACK that R uses, so a model links one implementation rather than two.
  # SUNDIALS requires a 32-bit index size for its LAPACK interface.
  #
  # SUNDIALS_LAPACK_CASE/UNDERSCORES state the Fortran name-mangling scheme
  # outright. Without them SUNDIALS infers it by compiling a Fortran program,
  # which needs a Fortran compiler this build has no other use for -- and on
  # a cluster where gfortran lives behind a module, that probe is what fails.
  # Every Fortran compiler R is built with mangles dgetrf to dgetrf_.
  #
  # This is attempted, not required: R may report no BLAS at all, the
  # libraries may be absent, or the try-compile may fail. The build then
  # proceeds without the LAPACK solver, and ./configure detects its absence
  # and leaves cvode(lapack = TRUE) reporting that it is unavailable.
  r_blas=$(R CMD config BLAS_LIBS 2>/dev/null || true)
  r_lapack=$(R CMD config LAPACK_LIBS 2>/dev/null || true)
  lapack_link=$(printf '%s %s' "$r_lapack" "$r_blas" | sed 's/^ *//; s/ *$//')

  sundials_lapack=0
  if [ -n "$lapack_link" ]; then
    log "trying LAPACK-backed dense solver against R's BLAS ($lapack_link)"
    write_lapack_shim "$lapack_link"
    if configure_sundials \
         -DENABLE_LAPACK=ON \
         -DSUNDIALS_INDEX_SIZE=32 \
         -DSUNDIALS_LAPACK_CASE=lower \
         -DSUNDIALS_LAPACK_UNDERSCORES=one \
         -DCMAKE_MODULE_PATH="$work/cmake-modules"; then
      sundials_lapack=1
    else
      log "LAPACK-enabled configure failed (see the cmake output above);"
      log "continuing without the LAPACK-backed dense solver"
    fi
  else
    log "R reports no BLAS library; building without the LAPACK solver"
  fi

  if [ "$sundials_lapack" = "0" ]; then
    configure_sundials || die "SUNDIALS configure failed"
  fi

  cmake --build "$work/sundials/build" -j "$njobs" >&2 || die "SUNDIALS build failed"
  cmake --install "$work/sundials/build" >&2 || die "SUNDIALS install failed"
  files_sundials || die "SUNDIALS installed but $PREFIX is incomplete"
  : > "$PREFIX/.sundials-complete"

  if ls "$PREFIX"/lib/libsundials_sunlinsollapackdense.* >/dev/null 2>&1; then
    log "LAPACK-backed dense solver available"
  else
    log "built without the LAPACK-backed dense solver"
  fi
fi

# Verify the files ./configure probes for rather than relying on the
# cmake exit status.
files_suitesparse || die "build finished but SuiteSparse/KLU is missing from $PREFIX"
if [ "$component" = "sundials" ]; then
  files_sundials || die "build finished but SUNDIALS is missing from $PREFIX"
  ls "$PREFIX"/lib/libsundials_sunlinsolklu.* >/dev/null 2>&1 \
    || die "SUNDIALS built without KLU support in $PREFIX"
fi

log "done: $PREFIX"
printf '%s\n' "$PREFIX"
