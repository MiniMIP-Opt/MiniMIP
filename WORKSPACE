workspace(name = "de_zib_minimip")

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository", "new_git_repository")

# Bazel Extensions
## Bazel Skylib rules.
git_repository(
    name = "bazel_skylib",
    tag = "1.5.0",
    remote = "https://github.com/bazelbuild/bazel-skylib.git",
)
load("@bazel_skylib//:workspace.bzl", "bazel_skylib_workspace")
bazel_skylib_workspace()

## Bazel rules.
git_repository(
    name = "platforms",
    tag = "0.0.9",
    remote = "https://github.com/bazelbuild/platforms.git",
)

git_repository(
    name = "rules_cc",
    tag = "0.0.9",
    remote = "https://github.com/bazelbuild/rules_cc.git",
)

git_repository(
    name = "rules_proto",
    tag = "5.3.0-21.7",
    remote = "https://github.com/bazelbuild/rules_proto.git",
)

# Support for building Python targets in Bazel.
git_repository(
    name = "rules_python",
    tag = "0.31.0",
    remote = "https://github.com/bazelbuild/rules_python.git",
)

# Dependencies
## ZLIB
new_git_repository(
    name = "zlib",
    build_file = "@com_google_protobuf//:third_party/zlib.BUILD",
    tag = "v1.2.13",
    remote = "https://github.com/madler/zlib.git",
)

## Abseil-cpp
git_repository(
    name = "com_google_absl",
    tag = "20240116.2",
    #patches = ["//patches:abseil-cpp-20240116.2.patch"],
    #patch_args = ["-p1"],
    remote = "https://github.com/abseil/abseil-cpp.git",
)

## Re2
git_repository(
    name = "com_google_re2",
    tag = "2024-04-01",
    remote = "https://github.com/google/re2.git",
    repo_mapping = {"@abseil-cpp": "@com_google_absl"},
)

## Protobuf
# proto_library, cc_proto_library, and java_proto_library rules implicitly
# depend on @com_google_protobuf for protoc and proto runtimes.
# This statement defines the @com_google_protobuf repo.
git_repository(
    name = "com_google_protobuf",
    #patches = ["//patches:protobuf-v26.1.patch"],
    #patch_args = ["-p1"],
    tag = "v26.1",
    #tag = "v25.3",
    remote = "https://github.com/protocolbuffers/protobuf.git",
)
# Load common dependencies.
load("@com_google_protobuf//:protobuf_deps.bzl", "protobuf_deps")
protobuf_deps()

## Python
load("@rules_python//python:repositories.bzl", "py_repositories")
py_repositories()

load("@com_google_protobuf//bazel:system_python.bzl", "system_python")
system_python(
    name = "system_python",
    minimum_python_version = "3.8",
)

## Testing
git_repository(
    name = "com_google_googletest",
    tag = "v1.14.0",
    remote = "https://github.com/google/googletest.git",
)

git_repository(
    name = "com_google_benchmark",
    tag = "v1.8.3",
    remote = "https://github.com/google/benchmark.git",
)

# Google OR-Tools. Contains many LP/MIP related utilities, including GLOP
# (LP solver) and SCIP (MIP solver).
#
# Note, OR-Tools includes some of the Google utilities not yet (fully) included
# in abseil (https://github.com/google/or-tools/tree/stable/ortools/base).
# In particular, OR-Tools includes logging routines, which are not part of
# abseil yet. Importantly, we do not depend on and use "glog", because it's old
# and depends on  "gflags" (now replaced by abseil/Flags).
git_repository(
    name = "com_google_ortools",
    #tag = "v9.9",
    tag = "v9.10",
    remote = "https://github.com/google/or-tools.git",
)

# bliss -- needed to compile SCIP with symmetry breaking.
http_archive(
    name = "bliss",
    build_file = "@com_google_ortools//bazel:bliss.BUILD.bazel",
    patches = ["@com_google_ortools//bazel:bliss-0.73.patch"],
    sha256 = "f57bf32804140cad58b1240b804e0dbd68f7e6bf67eba8e0c0fa3a62fd7f0f84",
    url = "http://www.tcs.hut.fi/Software/bliss/bliss-0.73.zip",
)

# OR-Tools v9.9 depends on scip v810
# OR-Tools v9.10 will depends on scip v900
new_git_repository(
    name = "scip",
    build_file = "@com_google_ortools//bazel:scip.BUILD.bazel",
    patches = ["@com_google_ortools//bazel:scip.patch"],
    #patches = ["@com_google_ortools//bazel:scip-v900.patch"],
    patch_args = ["-p1"],
    #tag = "v810",
    tag = "v900",
    remote = "https://github.com/scipopt/scip.git",
)

# Eigen is needed by ortools to compile the PDLP solver, which isn't used by
# MiniMIP but is included in targets we depend on.
http_archive(
    name = "eigen",
    url = "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip",
    strip_prefix = "eigen-3.4.0",
    sha256 = "1ccaabbfe870f60af3d6a519c53e09f3dcf630207321dffa553564a8e75c4fc8",
    build_file_content =
"""
cc_library(
    name = 'eigen3',
    srcs = [],
    includes = ['.'],
    hdrs = glob(['Eigen/**']),
    defines = ["EIGEN_MPL2_ONLY",],
    visibility = ['//visibility:public'],
)
""")

new_git_repository(
    name = "soplex",
    build_file = "//:soplex.BUILD",
    patches = ["//:soplex-bazel.patch"],
    patch_args = ["-p1"],
    tag = "release-710",
    remote = "https://github.com/scipopt/soplex.git",
)
